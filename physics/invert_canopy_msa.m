% This is the codes for the inversion of msa at canopy scale...


%% Importing the data
clear all 
close all
clc
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
I=find(RT(:,1)>400 & RT(:,1)<2400);
w_leaf = RT(I,1);
LS_sage=RT(I,2)+RT(I,3);
LS_bt=RT_BT(I,2)+RT_BT(I,3);

load('soilspectra.mat')
load('soil_wave.txt')
soil_ref=table2array(soilspectra(:,2:end))';
soil_resampled = interp1(soil_wave,soil_ref,w_leaf,'spline');
canopy_resampled = interp1(canopy_wave,Canopy_ref,w_leaf,'spline');
soil_id=T.Soil_id;
veg_type=T.VegType;
Gap=T.Gap;
I_gap=find(T.Gap~=-999)
%% This is inversionusing the original code provided by the author of the paper

[m,n]=size(canopy_resampled);
for i=1:n
   i
    SR=soil_resampled(:,soil_id(i));
    R=canopy_resampled(:,i);
    if (strcmp(veg_type{i},'Sagebrush'))
     LS=LS_sage;
    else
     LS=LS_bt;
    end
    p_med = msa2si(LS,SR,R,[150 1000]);
    P(i,:)=median(p_med);
    P_IC(1,1)=P(i,1);P_IC(1,2)=P(i,2);P_IC(2,1)=P(i,3);P_IC(2,2)=0;
    w=[LS SR];
    a=[P(i,4) 1-P(i,4)]';
    [R_est(:,i),C(:,:,i),Ic(:,:,i)] = msa(w,P_IC,a,1);
    [CB(:,:,i),CS(:,:,i),CBS(:,:,i),S(:,:,i),flag(i)] =CRB(P(i,:),LS,SR);
end
msa2_results.P=P;
msa2_results.R_est=R_est;
msa2_results.C=C;
msa2_results.IC=Ic;
msa2_results.CB=CB;
msa2_results.CS=CS;
msa2_results.CBS=CBS;
msa2_results.S=S;
msa2_results.flag=flag;
save msa2_results msa2_results
%% This is for the inversion of msa using cmaes algorithm. 
[m,n]=size(canopy_resampled);
for i=1:n
     i 
    SR=soil_resampled(:,soil_id(i));
    R=canopy_resampled(:,i);
    
    if (strcmp(veg_type{i},'Sagebrush'))
     LS=LS_sage;
    else
     LS=LS_bt;
    end
    opts.LBounds = 0; opts.UBounds = 1; 
    opts.Restarts=3;
    opts.Noise.on=1;
    [P(:,i),fmin(i)] = ...
        cmaes('chi2P5B_noise',0.5*ones(4,1),[],opts,LS,SR,R);
 
    P(1,1)=P(1,i);P(1,2)=P(2,i);P(2,1)=P(3,i);P(2,2)=0;
    w=[LS SR];
    a=[P(4,i) 1-P(4,i)]';
    j=1;
    [R_est(:,i),C(:,:,i),Ic(:,:,i)] = msa(w,P,a,j);
    %rms(i)=norm(R-R_est);
    [CB(:,:,i),CS(:,:,i),CBS(:,:,i),S(:,:,i),flag(i)] =CRB(P(:,i),LS,SR);
end

cmaes_noise.P=P;
cmaes_noise.R_est=R_est;
cmaes_noise.C=C;
cmaes_noise.IC=Ic;
cmaes_noise.CB=CB;
cmaes_noise.CS=CS;
cmaes_noise.CBS=CBS;
cmaes_noise.S=S;
cmaes_noise.flag=flag;
cmeas_noise.fmin= fmin;
save cmaes_noise cmaes_noise

