% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2P5B_noise(P,LS,SR,R)
 PM(1,1)=P(1);PM(1,2)=P(2);PM(2,1)=P(3);PM(2,2)=0;
 if (P(1)+P(2)>1 | 1-P(3)<0)  
 chi2=NaN;
 else
 w=[LS SR];
    a=[P(4) 1-P(4)]';
    j=1;
    [R_tmp,C_tmp,Ic] = msa(w,PM,a,j);
   noiseSigma_R = 0.04 * R_tmp;
   noise = noiseSigma_R .* (randn(1, length(R_tmp))');
   R_est = R_tmp + noise;
   
   noiseSigma_C = 0.04 * C_tmp;
   noise = noiseSigma_C .* (randn(1, length(C_tmp(:,1)))');
   C = C_tmp + noise;
   
   [CB,CS,CBS,S,flag] =CRB(P,LS,SR);
   
    if (any(C(:,1)<=0)|S(:,1)<0)
     chi2=NaN;
    else
chi2=norm(R-R_est);
    end
 end

