%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Gardens2 algorithm                                                  %
%     Jonás Grande Barreto                                                %
%     María Del Pilar Gómez Gil                                           %
%     https://doi.org/10.1007/s11517-020-02270-1                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inputs
% atlas_csf : atlas csf template 
% atlas_gm  : atlas gm template 
% atlas_wm  : atlas wm template 
% brain_msk : binari mask that indicates the whole brain volume
% feat_rep  : feature representation of the mri brain volume
% CI        : centroids for each tissue class


%% outputs
% partial maps : csf_u, gm_u, wm_u
clear;clc
close all

% NIfTI_files package can be found at
% https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
 %addpath(('NIfTI_files'))

%% Fixed parameters
% Tissues to segment (CSF, GM, & WM)
tissues = 3; 

% Range of scans where the brain is located a standard brain MRI volume 
% contains 200 scans. The brain is usually located in the middle of the 
% MRI volume
rt =[63,201];

% Level of fuzziness
expo = 2;

figure

% Load input files
CSF = load_nii('csf_template.nii');
GM  = load_nii('gm_template.nii');
WM  = load_nii('wm_template.nii'); 
brain_msk = load_nii('IBSR_01_ana_brainmask.nii');  
load('feats_Subject_x_nbf.mat','feat_all_ibsr')
load('Gardens2_hard_output_nbf.mat','OutPut','Cin_o','zave','CI')

for subject = 1 : 1
    %% Index to adjoining scans
    [col,row,dip] = size(brain_msk.img);     
    midd = round((rt(2)-rt(1))/2)+rt(1);
    slix = midd- 20 : midd + 29;
    point = rt(1) : rt(2);
    scan_length =length(slix);
    c1 = zeros(scan_length,15);
    c2 = zeros(scan_length,15);
    c3 = zeros(scan_length,15);
    K = zeros(tissues,15);
    initial_pure = zeros(row,col,scan_length);
    csf_u = zeros(row,col,scan_length);
    gm_u = zeros(row,col,scan_length);
    wm_u = zeros(row,col,scan_length);
%% Compute general centroids for the whole MRI volume
for qp = 1 : scan_length
    c1(qp,:) = CI(1,:,qp);
    c2(qp,:) = CI(2,:,qp);
    c3(qp,:) = CI(3,:,qp);
end
    
    K(1,:) = mean(c1,1,'omitnan');
    K(2,:) = mean(c2,1,'omitnan');
    K(3,:) = mean(c3,1,'omitnan');

    for scan = 1 : scan_length
        sliceg = slix(scan); 
        xslice = find(point == sliceg);
        %%  Match the tissue templates’ with the featurere presentation
        temp = zeros(row,col,tissues);
        tempa = zeros(row,col,tissues);
        mask = double(imrotate(round(abs(brain_msk.img(:,:,sliceg))),90)); 
        xs = find(mask);       
        csf = imrotate((CSF.img(:,:,sliceg)),90); 
        gm = imrotate((GM.img(:,:,sliceg)),90); 
        wm = imrotate((WM.img(:,:,sliceg)),90);    
        Z = zave(:,:,scan);         
        data = zeros(row*col,35); 
        data(xs,:) = feat_all_ibsr(1:length(xs),:,xslice);
        lnnan = isnan(data);
        data(lnnan) = eps;
        ld = 16:35;  
        data(:,ld) = [];   
        aux = normalize(data(xs,:)); 
        BI = reshape(data(:,13),row,col);
        Ci1=K; 
        
        U_new3  = stepfcm2(aux,Ci1,csf,gm,wm,mask,tissues,expo);          
         
        for qp = 1 : 3
            temp4 = zeros(row,col);
            temp4(xs) = U_new3(qp,:);
            tempa(:,:,qp) = temp4;
        end
        
        csf_u(:,:,scan) = tempa(:,:,1);
        gm_u(:,:,scan) = tempa(:,:,2);
        wm_u(:,:,scan) = tempa(:,:,3);
        
        subplot(2,2,1)
        imshow(BI,[],'InitialMagnification','fit')
        title ('MRI')
        subplot(2,2,2)
        imshow(tempa(:,:,1),[],'InitialMagnification','fit')
        title ('CSF partial map')
        subplot(2,2,3)
        imshow(tempa(:,:,2),[],'InitialMagnification','fit')
        title ('GM partial map')
        subplot(2,2,4)
        imshow(tempa(:,:,3),[],'InitialMagnification','fit')
        title ('WM partial map')
        pause(0.01)
    end
   save ('G2_partial_maps_nbf.mat','csf_u','gm_u','wm_u')
end


  