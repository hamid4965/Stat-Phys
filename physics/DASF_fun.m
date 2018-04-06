function [DASF,R2]=DASF_fun(w,BRF)
% w: leaf spectrum (700-900nm)
%BRF: canopy spectra (700-900nm)

ratio_dasf=BRF./w;
p=fitlm(BRF,ratio_dasf);
intercept=p.Coefficients.Estimate(1);
slope = p.Coefficients.Estimate(2);
DASF = intercept/(1-slope);

R2=p.Rsquared.Ordinary;
end