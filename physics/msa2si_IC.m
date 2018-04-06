function [p] = msa2si_IC(LS,SR,R,N,IC)
% MSA2SI - Two endmember spectral invariant test
%
% p = MSA2SI(LS,SR,R)
% p = MSA2SI(LS,SR,R,[Ns N])
% 
% Inputs:
%   LS - Leaf scattering signature
%   SR - soil reflectance signature
%   R - above canopy reflectance
%   Ns - number of wavelength samples used for parameter estimations
%   (default 30) 
%   N - number of estimations for each paramter of the MSA2 model (default
%   1000) 
%
% Outputs:
%  p = N-by-4 containing the estimated parameters [PLL PLS PSL i0];
%


if nargin < 4,
  N = 2000; Ns = 90; 
else
 Ns = N(1);IC=IC';
 if length(N) > 1,
   N = N(2);
 else
   N = 999;
 end
end
 LS = LS(:);
 SR = SR(:);
 R = R(:);
 M = length(LS);
 A = [R-SR 1-SR (1-R).*SR SR.*(1-LS)./LS];
 b = R./LS-SR;
 opt = optimset(optimset('fmincon'),'Algorithm','active-set','Display','off');
 %
 for n = 1:N,
    k = round(rand(1,Ns)*(M-1))+1;
    %c = max(min((A(k,:)'*A(k,:))\(A(k,:)'*b(k,1)),1),0);
    chi = [A(k,:) b(k,1)];
    c = fmincon(@(x) f(x,chi),IC,[],[],[],[],zeros(4,1),ones(4,1),[],opt);
    %pSL = max(roots([1-c(1)-c(2) c(2)-c(3)-(1-c(4))*(1-c(1)) (1-c(4))*c(3)]));
    %pLS = c(3)/pSL;
    %i0 = c(2)/(1-c(1)-pLS);
    %i0 = min(roots([1-c(1)+c(3) -(c(2)-c(3)+(1-c(4))*(1-c(1))) (1-c(4))*c(2)]));
    %roots([1-c(1)+c(3) -(c(2)-c(3)+(1-c(4))*(1-c(1))) (1-c(4))*c(2)])
    %pSL = 1-c(4)/(1-i0);
    %pLS = 1-c(1)-c(2)/i0;
    %p(n,:) = [c(1) pLS pSL i0];
    p(n,:) = c';
   
 end
 

function err = f(x,chi)
err = sqrt(mean((x(1)*chi(:,1)+(1-x(1)-x(2))*x(4)*chi(:,2)-x(2)*x(3)*chi(:,3)+(1-x(3))*(1-x(4))*chi(:,4)-chi(:,5)).^2,1));
