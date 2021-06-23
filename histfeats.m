function [feats]=histfeats(Iwindow);
    
    feats.Mean=0;
    feats.Variance=0;
    feats.Energy=0;
    feats.Entropy=0;
    feats.Skewness=0;
    feats.Kurtosis=0;
    feats.MPI=0;
    [pixelCounts, GLs] = imhist(uint8(Iwindow(:)));
    NM = sum(pixelCounts);     
    Prob  = pixelCounts./NM;
    feats.Mean     = sum(Prob.*GLs);
    feats.Variance=sum(Prob.*(GLs-feats.Mean).^2);
   
    term1    = Prob.*(GLs-feats.Mean).^3;
    term2    = sqrt(feats.Variance);
    feats.Skewness = term2^(-3)*sum(term1);
    
    term1    = Prob.*(GLs-feats.Mean).^4;
    term2    = sqrt(feats.Variance);
    feats.Kurtosis = term2^(-4)*sum(term1);
    
    feats.Energy   = sum(Prob.^2);
    feats.Entropy  = -nansum(Prob.*log2(Prob));   
    [~,feats.MPI]=max(Prob);
end
%     %energy
%     energy = sum((pixelCounts / NM) .^ 2);
%     %entropy 
%     pI = pixelCounts / NM;
%     entropy1 = -sum(pI(pI~=0) .* log2(pI(pI~=0)));
%     entropy1 =  -nansum(pI(sub).*log2(pI(sub)));%OK 9 
    %entropy2 = -sum(pI .* log2(pI+eps))
    
    
    