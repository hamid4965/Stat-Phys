function IC_opt=myinv(LB,UB,LS,SR,R)

%rmse=zeros(1000,1);

for k=1:100
p0=rand(100,4);
sol2(k,:)=fminsearchbnd('chi2P5B',p0(k,:),LB,UB,[],LS,SR,R);

 PM(1,1)=sol2(k,1);PM(1,2)=sol2(k,2);PM(2,1)=sol2(k,3);PM(2,2)=0;
 w=[LS SR];
 a=[sol2(k,4) 1-sol2(k,4)]';
rmse(k,1)=norm(R-msa(w,PM,a,1));
    
end
%boxplot(sol2)
%hold on
 a=(1:1:size(rmse,1))';
 rmset= [rmse a];
 rmse_sort=sortrows(rmset,1);
 index = rmse_sort(1,2);
 IC_opt = sol2(index,:);
 %IC_opt=mean(p_sort);
end