% _______________________________________________________________________
%
% chi2P5B.m
% merit function
% _______________________________________________________________________

function chi2=chi2P5B(P,LS,SR,R)
 PM(1,1)=P(1);PM(1,2)=P(2);PM(2,1)=P(3);PM(2,2)=0;
 if (P(1)+P(2)>1 | 1-P(3)<0)  
 chi2=NaN;
 else
 w=[LS SR];
    a=[P(4) 1-P(4)]';
    j=1;
    [R_est,C,Ic] = msa(w,PM,a,j);
     
    if (any(C(:,1)<=0))
     chi2=NaN;
    else
chi2=norm(R-R_est);
    end
 end

