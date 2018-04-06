function [CB,CS,CBS,S,flag] =CRB(p,leaf,soil)
% This script is based on the equation 13a-14c in:
%"Silván-Cárdenas, J.L., Corona-Romero, N., 2017. Radiation budget of 
%vegetation canopies with reflective surface: A generalization using the
%Markovian approach. Remote Sens. Environ. 189, 118–131. 
%doi:https://doi.org/10.1016/j.rse.2016.11.019"
% inputs : 
% p:probabilities calculated from msa2si.m
% LS = leaf reflectance
% SR = Soil reflectance
%outputs : 
% CB : canopy radiation budget
% CS : canopy radiation budget (S problem) 
% CBS = canopy radiation budget (S problem) 

n=length(p);
if n<4
    error('p must be 4;[PLL PLS PSL i0]')
end
PLL = p(1);
PLS = p(2);
PSL = p(3);
i0 = p(4);
t0 = 1-i0;

qL = 1-PLL-PLS;
qS = 1-PSL;

% solution to BS problem
rBS = (qL*leaf*i0)./(1-PLL*leaf);
aBS = (1-leaf)*i0./ (1-PLL*leaf);
tBS = t0+(PLS*leaf*i0./(1-PLL*leaf));

% solution to the S problem
rS = (tBS-t0)*PSL./i0;
aS = (aBS*PSL)./i0;
tS = (rBS*PSL)./i0+qS;

%soil contribution 
S1 = ((soil.*tBS).*(tS-qS))./(1-soil.*rS);
S2 = ((soil.*tBS).*aS)./(1-soil.*rS);
S3 = ((soil.*tBS).*(rS-PSL))./(1-soil.*rS);
% simulated canopy spectra
rC = rBS + S1;
aC = aBS + S2;
tC = tBS + S3;

CB = [rC aC tC];  % This is canopy spectra
CS = [rS aS tS];   % solution to S problem
CBS = [rBS aBS tBS];  %Solution to BS problem
S = [S1 S2 S3];       % Soil contribution to CRB
tmp1=isnan(CB);
tmp2=isnan(CS);
tmp3=isnan(CBS);
tmp4=isnan(S);
if (sum(sum(tmp1+tmp2+tmp3+tmp4))==0 & sum(CB(:,1)<0)==0)

flag=1;
else
    flag=0;
%return
end
end