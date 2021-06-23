function [idx, Ci, tsamp] = refsubRs(subRs,LIO,La,~,Atlas_crisp_n,subreg_cent,tissues)
    tsamp = zeros(1,3);
    [row,col] = size(LIO);
    subRsa = subRs;
    Ci = zeros(tissues,size(subRsa,2));
    idx = zeros(subreg_cent+1,tissues+2);    
      
    for k = 1 : tissues
        if k == 1
            colorin = 'red';
            Tissue = 'CSF';
        elseif k==2
            colorin = 'blue';
            Tissue = 'GM';
        else
            colorin = 'green';
            Tissue = 'WM';
        end
        lin = Atlas_crisp_n == k;
        clues_aux =zeros(row,col);
        clues_aux(lin) = k;

        ax1 = La&clues_aux;
        lio = ax1 == 1;
        Inreg = zeros(row,col);
        Inreg(lio) = La(lio);
        zvar = unique(Inreg);
        zvar(1) = [];        

        sup = zeros(length(zvar),4);
        sup(1:length(zvar),1) = subRsa(zvar,15);
        area_si = zeros(length(zvar),1);
        for a = 1 : length(zvar)
            a_ax = find(Inreg==zvar(a));
            area_si(a) = length(a_ax);
        end
        area_M = max(area_si);
        sup(1:length(zvar),2) = area_si;
        sup(1:length(zvar),3) = area_si/area_M; 
        sup(1:length(zvar),4) = sup(1:length(zvar),1)...
            ./sup(1:length(zvar),3);
        while sum(isnan(sup(:,4))) ~= length(zvar)
               [az,b] = max(sup(:,3));
                if az > 0.5 
                    [y,x] = find(La==zvar(b));
                    xy = round(length(y)/2);
                    text(x(xy),y(xy),'X','FontSize',8,'Color',colorin)
%                     title(sub), % [x(1),y(1)]
                     Cqk = subRsa(zvar(b),:);                      
                     idx(1,k) = zvar(b); 
                     sup(b,4) = nan;
                     subRsa(zvar(b),:) = nan; 
                    break
                else
                      sup(b,4) = nan; 
                      sup(b,3) = 0;
                end
            end
                                                                  % |
        disponible = 1 : length(zvar);                            % |
        subreg_ok = 0;                                         % |
        for qk = 1 : subreg_cent                                        % |  
            while ~isempty(disponible) &&  sum(isnan(sup(:,4))) ~= length(zvar)                           % |  
               dis = sum(Cqk.^2 + subRsa(zvar,:).^2 - 2*Cqk.*subRsa(zvar,:),2);
               [~,b] = min(dis,[],'omitnan');                                  % |                  
                l1 = find(La==zvar(b));                           % |
                    dim = length(l1);                             % |  
                    ax = La(l1);                                  % |  
                    ay = clues_aux(l1);                           % |  
                    aw = ax & ay;                                 % |
                    az = sum(aw)/dim;                             % |
                    if az > 0.7 && isempty(find(idx == zvar(b), 1)) % 
                       [y,x] = find(La==zvar(b));
                       xy = round(length(y)/2);
                        text(x(xy),y(xy),'X','FontSize',8,'Color',colorin)
%                         title(sub), %[x(1),y(1)]
                        idx(qk+1,k) = zvar(b);
                        sup(b,4) = nan;
                        subRsa(zvar(b),:) = nan;
                        subreg_ok = subreg_ok+1;                        
                        break
                    else
                        subRsa(zvar(b),:) = nan;
                    end
                    disponible(1)=[];
            end
        end
        
        disp([' Subregions for ', Tissue, ': ',num2str(subreg_ok),'/',num2str(subreg_cent)])

        tsamp(k) = subreg_ok/subreg_cent;
        lnx = idx(:,k);
        lnx0 = lnx == 0;
        lnx(lnx0) = [];
        Ci(k,:) = mean(subRs(lnx,:),1);
    end
end