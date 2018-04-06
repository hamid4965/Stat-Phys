clear all
clc
close all
cd('N:\Data02\bcal\Personal\hamid\StatPaper\final_working\PNAS\Bayesian')
%%
load plot_N_log.mat
load plotatr.mat

a1=plotatr.SiteName=="LP";
a2=plotatr.Year==2015;
Ind=logical(a1.*a2);
test_site=data(Ind,:);
cal_data=data(~Ind,:);

Ind2 = randperm(size(cal_data,1));
data_shuffle = cal_data(Ind2,:); 
Y=data_shuffle(:,1);
X=data_shuffle(:,2:size(data_shuffle,2));
Y_test=test_site(:,1);
X_test=test_site(:,2:size(test_site,2));
[N,p]=size(cal_data);
Indices = crossvalind('Kfold', N, 10);
%% Run the Bayesian model

   STANDARDISE=1; options.standardise=1;% Standarize the input.
   interaction=[0]; options.interaction = interaction; % our model is y ~ x1 + x2 +.... no interaction xi* xj is considered if set to 0
   order=[1]; options.order = order;% If set to 1, only the linear term is conisdered. No qudratic terms or higher-order terms like x^2 are included in the model
   alpha_1=0.01;alpha_2=0.01; options.alpha_1=alpha_1; options.alpha_2=alpha_2; %see below
   SAVE_SAMPLES = 1; options.save=SAVE_SAMPLES; 
% Please only adjsuts the following parameters   
   k_max = 100; options.k_max = k_max;% The maximum number of covaraites allowed in a model
   mcmc_samples = 60000; options.mcmc_samples = mcmc_samples; % the length of the MCMC chain
   burn_in = 10000; options.burn_in =burn_in; % the length of initial samples discarded as burn-in
for i=1:10
 
 y_val = Y(Indices==i);
 x_val = X(Indices==i,:);
 y_cal = Y(Indices~=i);
 x_cal = X(Indices~=i,:);
 [test_set_predictions,chain_stats,MODEL,mxsx,mnmx,SIG2] =...
     bayes_lm([y_cal, x_cal],[y_val, x_val],options)
 %[te_c, stat_c, MODEL_c, ms_c,mnmx_c, sig2_c] = bayes_lm([y_cal, x_In_cal],[y_val, x_val],options);
 [r2(i) cv(i)] = rsquare(y_val,test_set_predictions.pred_store)
end
   
[test_set_predictions_CS,chain_stats_CS,MODEL_CS,mxsx_CS,mnmx_CS,SIG2_CS] =...
     bayes_lm([Y, X],[Y_test, X_test],options);
[r2_CS cv_CS] = rsquare(Y_test,test_set_predictions_CS.pred_store);

LP15_plot_N_smooth.r2=r2;
LP15_plot_N_smooth.cv=cv;
LP15_plot_N_smooth.r2_CS=r2_CS;
LP15_plot_N_smooth.cv_CS=cv_CS;
save LP15_plot_N_smooth LP15_plot_N_smooth

mean(r2)
mean(cv)
r2_CS
cv_CS

   
   
   
   
   

