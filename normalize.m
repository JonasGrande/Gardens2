function shortfeat = normalize(allfeat) 
    shortfeat = zeros(size(allfeat,1),size(allfeat,2));
    a =allfeat==(inf);
    b =allfeat==(-inf);
    c = isnan(allfeat);
    allfeat(c) = eps;
    allfeat(a,1) =  1;
    allfeat(b,1) =  -1;
    [minv,~] = min(allfeat,[],'omitnan');
    [maxv,~] = max(allfeat,[],'omitnan');
    % uv=mean(allfeat,1);
    for k = 1 : size(allfeat,2)
         shortfeat(:,k) =(allfeat(:,k)-minv(k))/(maxv(k)-minv(k));
    end
%     shortfeat = normalize(allfeat',2)
end
