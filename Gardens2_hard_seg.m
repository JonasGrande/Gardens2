%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Gardens2 algorithm                                                  %
%     Jonás Grande Barreto                                                %
%     María Del Pilar Gómez Gil                                           %
%     https://doi.org/10.1007/s11517-020-02270-1                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear %all
close all
clc

% NIfTI_files package can be found at
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
 %addpath(('NIfTI_files'))

%% inputs
% atlas_csf : atlas csf template 
% atlas_gm  : atlas gm template 
% atlas_wm  : atlas wm template  
% brain_msk : binari mask that indicates the whole brain volume
% feat_rep  : feature representation of the mri brain volume

% Inputs example
%     atlas_csf = load_nii('Brain_X_DT\Subject_Y\csf_template.nii');
%     atlas_gm  = load_nii('Brain_X_DT\Subject_Y\gm_template.nii');
%     atlas_wm  = load_nii('Brain_X_DT\Subject_Y\wm_template.nii');  
%     GT = load_nii('Brain_X_DT\Subject_Y\brain_gt.nii');
%     brain_msk = load_nii('Brain_X_DT\Subject_Y\brainmask.nii');  
%     load('Brain_X_DT\Subject_Y\feats_Subject_Y.mat','feat_rep')

%% Outputs
%  OutPut : Brain tissue segmentation (CSF, GM, and WM)
%  Cin_o  : Centroids (CSF, GM, and WM)

%% Fixed parameters
% Tissues to segment (CSF, GM, & WM)
tissues = 3; 

% Range of scans where the brain is located a standard brain MRI volume 
% contains 200 scans. The brain is usually located in the middle of the 
% MRI volume
rt =[63,201];  %<---example 

% Number of subregions used to compute each centroid
subreg_cent = 10;   %<---- default

% Percent of inclution for each subregion
percent_inclu = 0.5;  %<---- default

% Colormap
map1 = [ 1   1   1      
        0.9 0.9  0    
        0.0 0.1 0.9    
        0.0 0.9 0.0     
         0.9  0  0.1      
         1   0   1];  

figure
for subject = 1 : 1
%% Load files
   atlas_csf = load_nii('csf_template.nii');
   atlas_gm  = load_nii('gm_template.nii');
   atlas_wm  = load_nii('wm_template.nii'); 
   brain_msk = load_nii('IBSR_01_ana_brainmask.nii');  
   load('feats_Subject_x_nbf','feat_all_ibsr')   
   
   %% Atlas' reference
   [Atlas_parcial, Atlas_crisp] = Atlas_ref(brain_msk,atlas_csf,atlas_gm,atlas_wm);
   
   %% Index to adjoining scans
  [col,row,dip] = size( brain_msk.img);     
   midd = round((rt(2)-rt(1))/2)+rt(1);
   slix = midd- 20 : midd + 29;
   point = rt(1) : rt(2);
   scan_length =length(slix);
   
   %% Main Algorithm
   OutPut = zeros(row,col,scan_length);
   zave = zeros(row,col,scan_length);
   CI = zeros(3,15,scan_length);
   % From scan x to y
   for scan = 1 : scan_length
       % Reshape feature representation to image format
       data = zeros(row*col,35); 
       sliceg = slix(scan); 
       xslice = find(point == sliceg);
%        Atlas_guia = imrotate(ag.img(:,:,sliceg),90);
       Atlas_crisp_n = imrotate(Atlas_crisp(:,:,sliceg),90);
       Atlas_parcial_n = imrotate(Atlas_parcial(:,:,sliceg),90);
       mask = double(imrotate(round(abs(brain_msk.img(:,:,sliceg))),90)); 
       xs = find(mask); 
       data(xs,:) = feat_all_ibsr(1:length(xs),:,xslice);
       lnnan = isnan(data);
       data(lnnan) = eps;
 
       ld = 16:35;  
       data(:,ld) = [];   
       BI = reshape(data(:,13),row,col);
       aux = normalize(data(xs,:)); 
       data(xs,:) = aux;
       imshow(BI,[],'InitialMagnification','fit')
       %% >>>>>>>>>>>>>>>>>> Oversegmentation  <<<<<<<<<<<<<<<%%        
        [J3,~] = imgradient(BI,'sobel');   
        L=watershed(J3); % Watershed transformation
        LIO=logical(L); %
        subregion=max(L(:)); % # de subregiones
        subRs=zeros(subregion,size(data,2)); 
        for k=1:subregion
            index=find(L==k);
            subRs(k,:)=mean(data(index,:),1);
        end

        % False color subregion
%         g0 = LIO==0;
%         AA = La;
%         AA(g0) = 0;
%         lbls = label2rgb(AA,'jet','k', 'shuffle');
%         lbls(g0) = 0;
%         figure,imshow(lbls,[],'InitialMagnification','fit')
         
         %remove background
         subRs(1,:) = []; 
         subregion = size(subRs,1);
         La= double(L);
         La = La - 1; 
         lnan = L == 1;
         ln0 = L == 0;
        % La(lnan) = -1;
         La(lnan) = 0;
         La(ln0) = 0;   
         g0 = LIO==0;
         AA = BI;
         AA(g0) = 0;
                 
         %% Centroid computation
        [idx, Ci, tsamp] = refsubRs(subRs,LIO,La,data,Atlas_crisp_n,...
        subreg_cent,tissues);               
        [suger,Sub_inx,Cin_o,Gardened] = suggestedsubregion(subRs,percent_inclu,Ci,...
        Atlas_parcial_n,La,tissues,data);     

        %% Clustering subregions
        Dis_n = zeros(subregion,tissues+2);
        RN = 1 : subregion;
        Cin_old = Cin_o;     
        for k = 1: subregion               
        while true
            r = randi([min(RN),max(RN)]);
            ri = find(RN==r);
            if ~isempty(ri)
                RN(ri) = [];
                break
            end
        end
        B = subRs(r,:);
        l2 = find(La == r); 
        A = Cin_old; 
        Dis_n(r,1:tissues) = sum(A.^2  - 2*A.*B,2)';  
        [~,lr] = min(Dis_n(r,1:tissues));   
        Dis_n(r,end) = lr;
        Sub_inx(r) =lr;   
        Gardened(l2) = lr;          
        end

      %% Clustering of isolated voxels (Watershed borders voxels)
      isolated_vox = mask & ~La;
      lf = find(isolated_vox);        
      A = Cin_o;
      for q = 1 : length(lf)
          B = data(lf(q),:);
          srd = sum(A.^2 - 2*A.*B,2);
          [~,srm] = min(srd);
          Gardened(lf(q)) = srm;
      end
   
   CI(:,:,scan) = Cin_o; 
   zave(:,:,scan) = Gardened;  
   OutPut(:,:,scan) = Gardened;
   
    subplot(1,2,1)
    imshow(AA,[],'InitialMagnification','fit')
    title ('Brain scan')
    subplot(1,2,2)
    imshow(uint8(Gardened),map1,'InitialMagnification','fit')
    title ('Segmentation')
    pause (0.01)
 
   end
   save ('Gardens2_hard_output_nbf.mat','OutPut','Cin_o','zave','CI')
end

   
   