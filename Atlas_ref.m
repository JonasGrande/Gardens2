


function [Atlas_parcial, Atlas_crisp] = Atlas_ref(brain_msk,atlas_csf,atlas_gm,atlas_wm)
%     atlas_csf = load_nii('Gardens2 pro\csf_template.nii');
%     atlas_gm  = load_nii('Gardens2 pro\gm_template.nii');
%     atlas_wm  = load_nii('Gardens2 pro\wm_template.nii'); 

    [row,col,dip] = size(brain_msk.img);
    Atlas_parcial = zeros(row,col,dip);
    Atlas_crisp = zeros(row,col,dip);
  
    
    for slice = 1 : dip
        CSF = atlas_csf.img(:,:,slice)*100;
        GM = atlas_gm.img(:,:,slice)*100;
        WM = atlas_wm.img(:,:,slice)*100;

        CSFt = CSF;
        GMt = GM;
        WMt = WM;
%         CSFD=CSFt>10;
%         GMD=GMt>10;
%         WMD=WMt>10;

%         D = CSFD+GMD+WMD;
        [rowa,cola]=size(CSF);
        mask_g = logical((CSF+WM+GM)*100);

        X = round(rowa/2);
        Y = round(cola/2);
        lim_x = round(X/3);   
        lim_y = round(Y/3);   
        [ra,ca] = find(CSF);

        for j=1:length(ra)
            if (abs((ra(j)-X)^2+(ca(j)-Y)^2)>lim_x^2)
                CSF(ra(j),ca(j))=0;
            end
        end

        [ra,ca] = find(WM);
        
        for j=1:length(ra)
            if  ((ra(j)-X)^2)/(lim_x-4)^2 + ((ca(j)-Y)^2)/(lim_y+2)^2 > 4
                WM(ra(j),ca(j))=0;
            end
        end  
 
        Atlas_parcial(:,:,slice) = parcial(CSFt,GMt,WMt,mask_g);
        Atlas_crisp(:,:,slice) = duro(CSFt,GMt,WMt,mask_g);
    end
%     if salvar == 1
%          nii2 = make_nii(B,[pixdim(1),pixdim(3),pixdim(2)]);   
%          nii3 = make_nii(C,[pixdim(1),pixdim(3),pixdim(2)]);     
%          save_nii(nii2,'Brain_X_DT\Subject_Y\Atlas_parcial.nii');
%          save_nii(nii3,'Brain_X_DT\Subject_Y\Atlas_crisp.nii');
%     end

end


function Atlas_guia = parcial(CSF,GM,WM,mask_g)
parcial_t = 20;
% umbral para máscaras fuzzy
brainlbl = find(mask_g);
    Atlas_guia = zeros(size(CSF));
    for pxls = 1 : length(brainlbl)
        lxd = brainlbl(pxls);                
        if  (CSF(lxd)>=parcial_t && GM(lxd)>=parcial_t) 
            Atlas_guia (lxd) = 4;
        elseif  (GM(lxd)>=parcial_t && WM(lxd)>=parcial_t) 
            Atlas_guia (lxd) = 5;
        elseif length(unique([CSF(lxd),GM(lxd),WM(lxd)]))>1
            [maxval,maxmem] = max([CSF(lxd),GM(lxd),WM(lxd)]);   
            if maxval > 20 
             Atlas_guia(lxd) = maxmem;
            end
        end
    end
end

function Atlas_guia_d = duro(CSF,GM,WM,mask_g)

brainlbl = find(mask_g);
    Atlas_guia_d = zeros(size(CSF));
    for pxls = 1 : length(brainlbl)
        lxd = brainlbl(pxls);                
        [maxval,maxmem] = max([CSF(lxd),GM(lxd),WM(lxd)]);   
        if maxval > 20 
            Atlas_guia_d(lxd) = maxmem;
        end
        
    end
    
    [rowa,cola] = size(CSF);
     X = round(rowa/2);
     Y = round(cola/2);
     lim_x = round(X/3);   
%      lim_y = round(Y/3);   

     [ra,ca] = find(Atlas_guia_d == 1);

        for j=1:length(ra)
            if (abs((ra(j)-X)^2+(ca(j)-Y)^2)>lim_x^2)
                Atlas_guia_d(ra(j),ca(j))=4;
            end
        end
    
    
    
end
