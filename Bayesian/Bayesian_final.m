clear all
clc
close all
cd('N:\Data02\bcal\Personal\hamid\StatPaper\final_working\PNAS\Bayesian')
%% Prepare data 

load plot_N_smooth.mat
%load plotatr.mat
%load leaf_wave_derv1.mat
%w=leaf_wave_derv1;
%rng('default')
%rng(1)
Ind = randperm(size(data,1));
data_shuffle = data(Ind,:); 
Y=data_shuffle(:,1);
X=data_shuffle(:,2:size(data_shuffle,2));
[N,p]=size(data);
Indices = crossvalind('Kfold', N, 10);
%% Use the x ~ y data above to test the bayesian linear model
%[test_set_predictions, chain_stats,MODEL,mxsx,mnmx,SIG2] = bayes_lm(data,test,options)
% INPUTS:
%		data - training data: first column is Y (the response) remaining columns are 
%									covariates, X
%		test - test data: same format as 'data'
%     options - this is an optional input object containing user defined settings for the 
%					program. You can run the program without this input if you wish to use 
%					our default settings. See below for further details.

% set up options
   STANDARDISE=1; options.standardise=1;% Standarize the input.
   interaction=[0]; options.interaction = interaction; % our model is y ~ x1 + x2 +.... no interaction xi* xj is considered if set to 0
   order=[1]; options.order = order;% If set to 1, only the linear term is conisdered. No qudratic terms or higher-order terms like x^2 are included in the model
   alpha_1=0.01;alpha_2=0.01; options.alpha_1=alpha_1; options.alpha_2=alpha_2; %see below
   SAVE_SAMPLES = 1; options.save=SAVE_SAMPLES; 
% Please only adjsuts the following parameters   
   k_max = 100; options.k_max = k_max;% The maximum number of covaraites allowed in a model
   mcmc_samples = 60000; options.mcmc_samples = mcmc_samples; % the length of the MCMC chain
   burn_in = 10000; options.burn_in =burn_in; % the length of initial samples discarded as burn-in

% In the following example, I use the same data for both traning and
% testing. Normally, they should be different. For example, you want to
% make predicitons for an indepedent dtest sample [yte, xte]
 %yte=y;
 %xte=x;

 
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
 %prediction_leaf_N_log_bayes.Model=MODEL;
 prediction_plot_LAI_derv1log_bayes.r2=r2;
 prediction_plot_LAI_derv1log_bayes.cv=cv;
 prediction_plot_LAI_derv1log_bayes.test_set_predictions=test_set_predictions;
 
save prediction_plot_LAI_derv1log_bayes.mat prediction_plot_LAI_derv1log_bayes
%mean(r2)
%mean(cv)

%% Let us look at the mean model cofficients
 % Often enough, MODEL is an large object array. So the following code is
 % to extract those subfields of MODEL we are only interested in. Here I
 % extracted the mean and sampled model cofficients of each iteration (
 % beta_mean and beta) as well as the indices of covariantes selected in
 % each sampled model (var).


%    STANDARDISE=1; options.standardise=1;% Standarize the input.
%    interaction=[0]; options.interaction = interaction; % our model is y ~ x1 + x2 +.... no interaction xi* xj is considered if set to 0
%    order=[1]; options.order = order;% If set to 1, only the linear term is conisdered. No qudratic terms or higher-order terms like x^2 are included in the model
%    alpha_1=0.01;alpha_2=0.01; options.alpha_1=alpha_1; options.alpha_2=alpha_2; %see below
%    SAVE_SAMPLES = 1; options.save=SAVE_SAMPLES; 
% Please only adjsuts the following parameters   
%    k_max = 100; options.k_max = k_max;% The maximum number of covaraites allowed in a model
%    mcmc_samples = 20000; options.mcmc_samples = mcmc_samples; % the length of the MCMC chain
%    burn_in = 10000; options.burn_in =burn_in; % the length of initial samples discarded as burn-in
   
   
 y_te=Y;
 x_te=X;
 [test_set_VS,chain_stats_VS,MODEL_VS,mxsx_VS,mnmx_VS,SIG2_VS]=...
     bayes_lm([y_cal, x_cal],[y_te, x_te],options)  
 
  ML={};
   for j=1:numel(MODEL_VS)
       M=MODEL_VS{j};
       ML{j}.var=[M.basis.var];
       ML{j}.beta=M.beta;
       ML{j}.beta_mean=M.beta_mean;  
   end
  % clear MODEL
   % average over the sampled beta_mean to get the mean model coefficients
   beta_bma=zeros(1,size(X,2)+1);
   for j=1:numel(ML)
       M=ML{j};
       cov_indx=[M.var];
       B=M.beta_mean;
       for jj=1:length(cov_indx)
           jjj=cov_indx(jj);
           if jjj == 0
               beta_bma(1)=beta_bma(1)+B(1);
           else
               beta_bma(jjj)=beta_bma(jjj)+B(jj);
           end
       end
       
       for kk=1:(p-1)
          sel(kk,j) = any(cov_indx == kk);
           
       end
       
   end
  s = sum(sel,2);
  p_sel = s/numel(ML);
  ave_p=mean(p_sel);
  std_p=std(p_sel);
  Ind_sel=p_sel>ave_p+std_p;
beta_bma=beta_bma/numel(ML);
clf
plot(beta_bma(2:end),'-r')

%selected_bayes.MODEL_VS=MODEL_VS;
selected_bayes.p_sel = p_sel;
selected_bayes.beta = beta_bma(2:end);
save SelectedBands_plot_LAI_derv1log_bayes selected_bayes








