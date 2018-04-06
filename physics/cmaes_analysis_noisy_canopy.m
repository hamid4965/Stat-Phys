% Analyze the results of cmaes with noisy objective function for the
% canopy scale.... the loaded cmaes_noise_caopy. mat is the results of the optim 
clear all
close all
clc
load cmaes_noise_canopy.mat
%load cmaes_canopy_all.mat
T=readtable('Spectrum.xlsx');
colName = T.Properties.VariableNames;
wave_tmp=colName(18:end);
for i=1:size(wave_tmp,2)
   a= wave_tmp(i);
   canopy_wave(i)=str2num(erase(a{1},'x'));
end
Canopy_ref= table2array(T(:,18:end))';
load RT_sim % This is sage leaf endmember send from Dar
load RT_sim_bt % This is bitterbrush leaf endmember send from Dar
load LRT_DASF.mat  % This is simulated by prospect model
I=find(RT(:,1)<2400);
w_leaf = RT(I,1);
LS_sage=RT(I,2)+RT(I,3);
LS_bt=RT_BT(I,2)+RT_BT(I,3);
load('soilspectra.mat')
load('soil_wave.txt')
soil_ref=table2array(soilspectra(:,2:end))';
soil_resampled = interp1(soil_wave,soil_ref,w_leaf,'spline');
canopy_resampled = interp1(canopy_wave,Canopy_ref,w_leaf,'spline');
LS_knya = interp1(LRT(:,1),LRT,w_leaf,'spline');

soil_id=T.Soil_id;
veg_type=T.VegType;
ID=T.ID;
C=cmaes_noise.C;
R_est=cmaes_noise.R_est;
P=cmaes_noise.P;
CB=cmaes_noise.CB;
flag=cmaes_noise.flag';
S=cmaes_noise.S;

for i=1:size(canopy_resampled,2)
  rms(i)=norm(canopy_resampled(:,i)-R_est(:,i));
  rmse(i)=rms(i)/sqrt(size(canopy_resampled,1))
end 
mean_rmse=mean(rmse);
hist(rmse)
% plot(R_est(:,95),'r')
% hold on
% plot(canopy_resampled(:,95),'b')
% hist(rmse)
% xlabel('Root mean square error')
% ylabel('Frequency')
% set(gcf,'color','w','Position', [250, 150,700,600]);
%  set(gca,'FontSize',16)
ID_sel = ID(logical(flag));
C_sel=C(:,:,ID_sel);
S_sel=S(:,:,ID_sel);
P_sel=P(:,ID_sel);
N_sel=T.PercentN(ID_sel);
ShrubName_sel=T.ShrubName(ID_sel);
LAI_Sel=T.LAI(ID_sel);
veg_type_sel=veg_type(ID_sel);
Gap_sel=T.Gap(ID_sel);
Height_sel=T.Height(ID_sel);
rmse_sel=rmse(ID_sel);
R_est_sel=R_est(:,ID_sel);
total_soil=R_est(:,ID_sel) - squeeze(C_sel(:,1,:));
canopy_resampled_sel=canopy_resampled(:,ID_sel);


%% 


conv=sum(flag);
boxplot(P_sel','Labels',{'PLL','PLS','PSL','i0'})
 set(gcf,'color','w','Position', [250, 150,700,600]);
 set(gca,'FontSize',16,'linew',1.5)
%  title(['Estimated ; ' 'converged=' num2str(conv) ' ' 'not converged='...
%      num2str(151-conv)])
 set(findobj(gca,'type','line'),'linew',1.5)
 ylabel('Probability');
 xlabel('Invariant parameters')
 ylim([0 1])
 

hold on
%title('LUT- radiation budget, mean rms=0.93')
h1=plot(w_leaf,R_est_sel,'r');
box on
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',12,'linew',1.5)
  xlim([350,2450])
  ylim([0 0.5])
  legend([h1(1)], 'total reflectance (simulated)');
xlabel('Wavelenght [nm]')
ylabel('Reflectance')
h2=plot(w_leaf,squeeze(C_sel(:,1,:)),'green');
legend([h1(1), h2(1)], 'total reflectance (simulated)', 'Canopy');
h3=plot(w_leaf,squeeze(S_sel(:,1,:)),'cyan');
legend([h1(1), h2(1),h3(1)], 'total reflectance (simulated)', 'Canopy',...
    'soil contributed in CRB');
h4=plot(w_leaf,total_soil,'k');
legend([h1(1), h2(1),h3(1),h4(1)], 'total reflectance (simulated)', 'Canopy',...
    'soil contributed in canopy','total soil');
 
 for i=1:size(ID_sel,1)
 hold on
 h1=plot(w_leaf,canopy_resampled_sel(:,i),'r',w_leaf,R_est_sel(:,i),'b');
h4=plot(w_leaf,soil_resampled(:,7),'k--')
 h3=plot(w_leaf,C_sel(:,1,i),'green',w_leaf,S_sel(:,1,i),'k');
 
 set(h1,'LineWidth',2)
 set(h3,'LineWidth',2)
 set(h4,'LineWidth',2)
 box on
 xlim([350,2450])
 ylim([0 0.7])
 xlabel('Wavelenght[nm]')
 ylabel('Reflectance')
 
 set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',12,'LineWidth',1.5)
 title(['LAI = ' num2str(LAI_Sel(i)) ';' 'rms = ' num2str(rms(i))...
     ';' 'veg = ' veg_type_sel{i}])
 text(500,0.6,{'PLL:','PLS:','PSL:','i0:'}, 'FontSize', 12)
 text(700,0.6,{num2str(round(P_sel(:,i),2))}, 'FontSize', 12)
  legend('Measured reflectance','Simulated reflectance',...
       'Soil reflectance','Canopy refelectance','soil contribution to CRB')
i
   pause
     close all
 end
 
 % plot the very CRB seperatley
 h1=plot(w_leaf,C_sel(:,1,1),'green',w_leaf,S_sel(:,1,1),'k');
 legend('Canopy reflectance','Soil contribution to CRB')
 xlabel('Wavelenght[nm]')
 ylabel('Reflectance')
 set(h1,'LineWidth',2)
 set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',12,'LineWidth',1.5)
 ShrubName_sel(96)
  
  
 %% Calculate DASF
 
 I_dasf=find(w_leaf>710 &w_leaf<790);
 LS_DASF_sage=LS_sage(I_dasf);
 LS_DASF_bt=LS_bt(I_dasf);
 LS_DASF_kyna=LS_knya(I_dasf,2)+LS_knya(I_dasf,3);
 S_ref_sel=squeeze(S_sel(:,1,:));
 for  i=1:size(ID_sel,1)
  
   BRF = squeeze(C_sel(I_dasf,1,i)); % if not correcting for S
   %BRF = squeeze(C_sel(I_dasf,1,i))-S_ref_sel(I_dasf,i);
   %BRF = canopy_resampled_sel(I_dasf,i); %if not correcting for soil
   if (strcmp(veg_type{i},'Sagebrush'))
   LS=LS_DASF_sage;
    else
     LS=LS_DASF_bt;
    end
LS_pnas=LS_DASF_kyna;

% ratio1= BRF./LS_DASF_sage;
% ratio2= BRF./LS_DASF_kyna;
% p1=fitlm(BRF,ratio1);
% p2=fitlm(BRF,ratio2);
% intercept1=p1.Coefficients.Estimate(1)
% intercept2=p2.Coefficients.Estimate(1)
% slope1 = p1.Coefficients.Estimate(2)
% slope2 = p2.Coefficients.Estimate(2)
% DASF1=intercept1/(1-slope1)
% DASF2=intercept2/(1-slope2)
[DASF(i),R2(i)]=DASF_fun(LS,BRF);
[DASF_pnas(i),R2_pnas(i)]=DASF_fun(LS_pnas,BRF);
%WC_all(:,i)=canopy_resampled_sel(:,i)/DASF(i);
%WC_all(:,i)= C_sel(:,1,i)/DASF(i);
WC_all(:,i)= (squeeze(C_sel(:,1,i))-S_ref_sel(:,i))/DASF(i);
 WC_all_pnas(:,i)= (squeeze(C_sel(:,1,i))-S_ref_sel(:,i))/DASF_pnas(i);

  if (sum(WC_all(:,i)<0)>0 | sum(WC_all(:,i)>1)>0|R2(i)<0.9|DASF(i)<=0)
      flag2(i)= 0;
  else
      flag2(i)=1;
  end
  
   if (sum(WC_all_pnas(:,i)<0)>0 | sum(WC_all_pnas(:,i)>1)>0|R2_pnas(i)<0.9)
      flag3_pnas(i)= 0;
  else
      flag3_pnas(i)=1;
  end
  
 end
 

 plot(DASF,DASF_pnas,'*')
 %ylim([0 0.5])
 %xlim([0 0.5])
%%

hist(R2(logical(flag2)))
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 xlabel('R2 values')
 ylabel('Frequency')

plot(w_leaf,WC_all(:,logical(flag2)))
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 xlabel('Wavelenght [nm]')
 ylabel('Reflectance')

DASF1=DASF(logical(flag2));
boxplot(DASF1)

shrubname1=ShrubName_sel(logical(flag2));
LAI1=LAI_Sel(logical(flag2));
P1=P_sel(:,logical(flag2));
N1=N_sel(logical(flag2));
Gap1=Gap_sel(logical(flag2));
WC1=WC_all(:,logical(flag2));
R1=canopy_resampled_sel(:,logical(flag2));
R_nosoil1=C_sel(:,1,logical(flag2));
Height1=Height_sel(logical(flag2));
I_Good=find(N1~=-999);
final_N=N1(I_Good);
Gap_final=Gap1(I_Good);
LAI_final=LAI1(I_Good);
Height_final=Height1(I_Good);
WC_final=WC1(:,I_Good);
R_final=R1(:,I_Good);
R_nosoil_final=R_nosoil1(:,I_Good);
P_final=P1(:,I_Good);
DASF_final=DASF1(I_Good);
shrubname_final=shrubname1(I_Good);

I_NIR=find(w_leaf>800 & w_leaf<850);
WC_NIR= mean(WC_final(I_NIR,:),1);
R_NIR=mean(R_final(I_NIR,:),1);
R_nosoil_NIR=mean(R_nosoil_final(I_NIR,:,1));


p1 = polyfit(final_N,WC_NIR',1); 
f = polyval(p1,final_N); 
p2=fitlm(final_N,WC_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(final_N,WC_NIR,'o',final_N,f,'-');
xlabel('N [g/100g]')
ylabel('NIR (800-850 [nm])')
text(1.8,0.7,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 
p1 = polyfit(LAI_final,R_NIR',1); 
f = polyval(p1,LAI_final); 
p2=fitlm(LAI_final,R_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(LAI_final,R_NIR,'o',LAI_final,f,'-');
xlabel('LAI [m/m]')
ylabel('BRF-NIR (800-850 [nm])')
text(3,0.35,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 
 p1 = polyfit(Height_final,R_NIR',1); 
f = polyval(p1,Height_final); 
p2=fitlm(Height_final,R_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(Height_final,R_NIR,'o',Height_final,f,'-');
xlabel('Height [m/m]')
ylabel('BRF-NIR (800-850 [nm])')
text(1.5,0.35,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)

  p1 = polyfit(DASF_final,R_NIR,1); 
f = polyval(p1,DASF_final); 
p2=fitlm(DASF_final,R_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(DASF_final,R_NIR,'o',DASF_final,f,'-');
xlabel('DASF')
ylabel('BRF-NIR (800-850 [nm])')
text(0.4,0.35,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)

 %% Export structure removed BRF to R
DASF_removed_canopy = [final_N LAI_final WC_final'];
save DASF_removed_canopy DASF_removed_canopy



%% check above using prospect leaf BRF

hist(R2_pnas(logical(flag3_pnas)))
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 xlabel('R2 values')
 ylabel('Frequency')

plot(w_leaf,WC_all_pnas(:,logical(flag3_pnas)))
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 xlabel('Wavelenght [nm]')
 ylabel('Reflectance')

DASF1=DASF_pnas(logical(flag3_pnas));
boxplot(DASF1)

shrubname1=ShrubName_sel(logical(flag3_pnas));
LAI1=LAI_Sel(logical(flag3_pnas));
P1=P_sel(:,logical(flag3_pnas));
N1=N_sel(logical(flag3_pnas));
Gap1=Gap_sel(logical(flag3_pnas));
WC1=WC_all(:,logical(flag3_pnas));
R1=canopy_resampled_sel(:,logical(flag3_pnas));
R_nosoil1=C_sel(:,1,logical(flag3_pnas));
Height1=Height_sel(logical(flag3_pnas));
I_Good=find(N1~=-999);
final_N=N1(I_Good);
Gap_final=Gap1(I_Good);
LAI_final=LAI1(I_Good);
Height_final=Height1(I_Good);
WC_final=WC1(:,I_Good);
R_final=R1(:,I_Good);
R_nosoil_final=R_nosoil1(:,I_Good);
P_final=P1(:,I_Good);
DASF_final=DASF1(I_Good);
shrubname_final=shrubname1(I_Good);

I_NIR=find(w_leaf>800 & w_leaf<850);
WC_NIR= mean(WC_final(I_NIR,:),1);
R_NIR=mean(R_final(I_NIR,:),1);
R_nosoil_NIR=mean(R_nosoil_final(I_NIR,:,1));


p1 = polyfit(final_N,WC_NIR',1); 
f = polyval(p1,final_N); 
p2=fitlm(final_N,WC_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(final_N,WC_NIR,'o',final_N,f,'-');
xlabel('N [g/100g]')
ylabel('NIR (800-850 [nm])')
text(1.8,0.7,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 
p1 = polyfit(LAI_final,R_NIR',1); 
f = polyval(p1,LAI_final); 
p2=fitlm(LAI_final,R_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(LAI_final,R_NIR,'o',LAI_final,f,'-');
xlabel('LAI [m/m]')
ylabel('BRF-NIR (800-850 [nm])')
text(3,0.35,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)
 
 p1 = polyfit(Height_final,R_NIR',1); 
f = polyval(p1,Height_final); 
p2=fitlm(Height_final,R_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(Height_final,R_NIR,'o',Height_final,f,'-');
xlabel('Height [m/m]')
ylabel('BRF-NIR (800-850 [nm])')
text(1.5,0.35,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)

  p1 = polyfit(DASF_final,R_NIR,1); 
f = polyval(p1,DASF_final); 
p2=fitlm(DASF_final,R_NIR);
R2_corr=p2.Rsquared.Ordinary;
plot(DASF_final,R_NIR,'o',DASF_final,f,'-');
xlabel('DASF')
ylabel('BRF-NIR (800-850 [nm])')
text(0.4,0.35,['R2 =' ' ' num2str(R2_corr)])
set(gcf,'color','w','Position', [250, 150,600,500]);
 set(gca,'FontSize',16)









