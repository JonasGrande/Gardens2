
function [suger,Sub_inx,Cin_o,Gardened] = suggestedsubregion(subRs,percent_inclu,Ci,Atlas_parcial_n,La,tissues,data)      
        suger = zeros(length(subRs),3);
        Dis = zeros(length(subRs),3);
        Sub_inx = zeros(length(subRs),1);
        [row,col] = size(La);
        old = Ci;
        err3 = 0.00001;
        distan = 0.5;
%         subRsa = subRs;
        Gardened = zeros(row,col);
        for yu = 1 : 100
        for tk = 1 : tissues
            lin = Atlas_parcial_n == tk;
            clues_aux =zeros(row,col);
            clues_aux(lin) = tk;
            ax1 = La&clues_aux;
            lio = ax1 == 1;
            Inreg = zeros(row,col);
            Inreg(lio) = La(lio);
            zvar = unique(Inreg); 
            zvar(1) = [];                       
            tsc = 0;
            for ts = 1 : length(zvar)            % |  
                l1 = find(La==zvar(ts));          % |  
                dim = length(l1);                % |  
                ax = La(l1);                     % |  
                ay = clues_aux(l1);              % |  
                aw = ax & ay;                    % |  
                az = sum(aw)/dim;                % |  
                A = Ci(tk,:);
                B = subRs(zvar(ts),:);
                dis = sum(A.^2 + B.^2 - 2*A.*B,2);
                if az > percent_inclu &&  isempty(find(suger == zvar(ts), 1))...
                        && abs(dis)<=distan
                     tsc = tsc +1 ;   
                     suger(tsc,tk) = zvar(ts); 
                     Dis(tsc,tk) = dis; 
                     Sub_inx(zvar(ts)) = tk; 
                     Gardened(l1) = tk;
                     Cin_o = centroid_ci(data,Gardened);
                end
            end
        end
         dif_0 = sum(abs(Cin_o-old),2)';
         dif_t = dif_0 < err3;
         if (isequal(dif_t,ones(1,3))) % ||  yu > 50    %
             break
         else
         old = Cin_o;          
         end
        end
end

function Ci = centroid_ci(data,Gardened)
    Ci = zeros(3,size(data,2));
    for u = 1 : 3
         lt = find(Gardened == u);
         Ci(u,:) = sum(data(lt,:))./length(lt);
    end
end
