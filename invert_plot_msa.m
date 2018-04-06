% This is the final version of the optimization for PLOTS
% after trying many different methoods I decided to use CMAES algorithm for
% the optimization. In this version I first run the cmaes randomley 10
% times then the initial point with miimum RMSE would be used to run the
% cmaes. 
clear all
close all
clc
T=readtable('Plot_Spectrum.xlsx');
colName = T.Properties.VariableNames;
wave_tmp=colName(15:end);
for i=1:size(wave_tmp,2)
   a= wave_tmp(i);
   tmp1=erase(a{1},'x');
  plot_wave(i)= str2num(strrep(tmp1,'_','.'));
end
plot_ref= table2array(T(:,15:end))';
load RT_sim % This is sage leaf endmember send from Dar
load RT_sim_bt % This is bitterbrush leaf endmember send from Dar

I=find(plot_wave>400 & plot_wave<2400);
w_plot = plot_wave(I);
w_leaf = RT(:,1);
LS_sage=RT(:,2)+RT(:,3);
LS_bt=RT_BT(:,2)+RT_BT(:,3);

load('soilspectra.mat')
load('soil_wave.txt')
soil_ref=table2array(soilspectra(:,2:end))';
tmp2=mean(soil_ref(:,1:13),2);
soil_ref(:,16)=tmp2;

soil_resampled = interp1(soil_wave,soil_ref,w_plot,'spline');
leaf_sage_resampled = interp1(w_leaf,LS_sage,w_plot,'spline');
leaf_bt_resampled = interp1(w_leaf,LS_bt,w_plot,'spline');
soil_id=T.Soil_id;
veg_type=T.VegType;

%% Cmaes plot noise 
[m,n]=size(plot_ref);
for i=1:n
     SR=soil_resampled(:,soil_id(i));
    R=plot_ref(:,i)';
    if (strcmp(veg_type{i},'Bitterbrush'))
     LS=leaf_bt_resampled';
    else
     LS=leaf_sage_resampled';
    end
    
    opts.LBounds = 0; opts.UBounds = 1; 
    opts.Restarts=1; opts_int.MaxFunEvals=10000
    opts.Noise.on=3;
    [P(i,:),fmin(i)] = ...
        cmaes('chi2P5B_noise',0.5*ones(4,1),[],opts,LS,SR,R);
 
    P_est(1,1)=P(i,1);P_est(1,2)=P(i,2);P_est(2,1)=P(i,3);P_est(2,2)=0;
    w=[LS SR];
    a=[P(i,4) 1-P(i,4)]';
   
    [R_est(:,i),C(:,:,i),Ic(:,:,i)] = msa(w,P_est,a,1);
    %rms(i)=norm(R-R_est);
    [CB(:,:,i),CS(:,:,i),CBS(:,:,i),S(:,:,i),flag(i)] =CRB(P(i,:),LS,SR);
end

cmaes_noise_plot.P=P;
cmaes_noise_plot.R_est=R_est;
cmaes_noise_plot.C=C;
cmaes_noise_plot.IC=Ic;
cmaes_noise_plot.CB=CB;
cmaes_noise_plot.CS=CS;
cmaes_noise_plot.CBS=CBS;
cmaes_noise_plot.S=S;
cmaes_noise_plot.flag=flag;
cmeas_noise.fmin= fmin;
save cmaes_noise_plot cmaes_noise_plot

%% CMAES with noise with 10 iteration

% [m,n]=size(plot_ref);
% for i=1:n
% 
%     SR=soil_resampled(:,soil_id(i));
%     R=plot_ref(:,i)';
%     if (strcmp(veg_type{i},'Bitterbrush'))
%      LS=leaf_bt_resampled';
%     else
%      LS=leaf_sage_resampled';
%     end
%     
%     for j=1:10
%         i
%         j
%     IC(j,:)=rand(4,1)    
%     opts_int.LBounds = 0; opts_int.UBounds = 1; 
%     opts_int.Noise.on=1;  opts_int.MaxFunEvals=10000;
%     [P_int(:,j),fmin_int(j,1)] = ...
%         cmaes('chi2P5B_noise',IC(j,:),[],opts_int,LS,SR,R);
%  fmin_int(j,2)=j;
%     end
%     
%  fmin_sort=sortrows(fmin_int,1);
%  index = fmin_sort(1,2);
%  IC_opt = P_int(:,index);
%  opts.LBounds = 0; opts.UBounds = 1; 
%  opts.Restarts=2;  opts.MaxFunEvals=10000;
%  opts.Noise.on=1;
%  [P(:,i),fmin(i)] = cmaes('chi2P5B_noise',IC_opt,[],opts,LS,SR,R);
%  
%     P_out(1,1)=P(1,i);P_out(1,2)=P(2,i);P(2,1)=P(3,i);P_out(2,2)=0;
%     w=[LS SR];
%     a=[P(4,i) 1-P(4,i)]';
%    [R_est(:,i),C(:,:,i),Ic(:,:,i)] = msa(w,P_out,a,1);
%     %rms(i)=norm(R-R_est);
%     [CB(:,:,i),CS(:,:,i),CBS(:,:,i),S(:,:,i),flag(i)] =CRB(P(:,i),LS,SR);
% end
% 
% cmaes_plot_all.P=P;
% cmaes_plot_all.R_est=R_est;
% cmaes_plot_all.C=C;
% cmaes_plot_all.IC=Ic;
% cmaes_plot_all.CB=CB;
% cmaes_plot_all.CS=CS;
% cmaes_plot_all.CBS=CBS;
% cmaes_plot_all.S=S;
% cmaes_plot_all.flag=flag;
% cmaes_plot_all.fmin= fmin;
% save cmaes_plot_all cmaes_plot_all
% 
% %filepath =('/home/hdashti/scratch/msa/msa3/cmaes_plot_all.mat');
% %save(filepath,'cmaes_plot_all')
% %exit
% 
