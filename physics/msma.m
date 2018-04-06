function [F,X] = msma(Y,X0,F0,COVY,COVX,COVF,OPT)
% MSMA - Spectral mixture analysis w endmember optimization
%
% Synopsis:
%  [F,X] = MSMA(Y,X0)
%  [F,X] = MSMA(Y,X0)
%  [F,X] = MSMA(Y,X0,F0,COVY,COVX,COVF,OPT)
%
% Inputs:
%  Y - N-by-P matrix of P mixed pixels from a multispectral image of N
%    spectral bands
%  X0 - N-by-M matrix with M initial endmembers of N spectral bands.
%  F0 - M-by-1 or M-by-K matrix of K initial fraction vectors. Default is
%     least square estimate. You may specify an empty matrix to use
%     default.
%  COVY - N-by-N matrix with variance/covariance of data Y. Use a scalar to
%    specify zero covariances and constant variances. Default is 
%    COVY = 0.01. You may specify an empty matrix to use default.
%  COVX - M-by-M matrix with variance/covariance of endmembers. Use a
%    scalar to specify zero covariances and constant variances. Default is 
%    COVX = 0.01. You may specify an empty matrix to use default.
%  COVF - M-by-M matrix with variance/covariance of endmembers. Use a
%    scalar to specify zero covariances and constant variances. Default is 
%    COVF = 0.01. You may specify an empty matrix to use default.
%  OPT - 4-length vector with the following optimization parameters:  
%   OPT(1) = maximum number of iterations (Default 50)
%   OPT(2) = minimum change in X at each iteration (Default 1e-6)
%   OPT(3) = minimum change in F at each iteration (Default 1e-6)
%   OPT(4) = minimum RMSE of predicted pixels (Default 1e-6)
%
% Output:
%  F - M-by-K matrix of K fraction vectors for the M end members.
%  X - N-by-M matrix with M final (optimized) endmembers of N spectral
%      bands. 
%
% Note: This method uses damped least square nonlinear inversion to
% simulataneously estimate the fractions and endmembers of a linear
% mixture model formed for every pixels. Thus, the total parameters
% estimated are M*(N+1+P) whereas the total equations formed are P*(N+1)
% 
% Reference:
% Stefanie Tompkins, John F. Mustard, Carld M. Pieters, and Donald W. 1997
% Forsyth "Optimization of Endmembers Mixture Analysis". Remote Sensing of
% Environment 59:472-489.
%
% Example: This example shows use MSMA to estimate an unknon endmember. The
% mixed pixels are synthesized from three end members using random
% fractions. The endmember X(:,1) is replaced by the mixed spectra Y(:,1)
% for the initial solution.
% 
%  load FieldSpecRef.mat                   % load field spectral reflectance 
%  [n,m] = size(X);                      % size of endmembers
%  p = 100;                                % number of mixed pixels
%  F = rand(m,p);                          % true fractional cover
%  F(:) = F./repmat(sum(F,1),m,1);         % fractions must sum to one
%  Y = X*F+0.005*randn(n,p);                % simulated mixed pixels w gaussian noise
%  [Fsma,Xmsa] = msma(Y,[Y(:,1),X(:,2:3)],[],cov(Y'),[5 .1 .1],[5 .1 .1]);
%  plot(w,X,'-',w,Xmsa,':')
%  xlabel('Wavelength [nm]')
%  ylabel('Reflectance')

minF = -0.05;       % zero minus tolerance
maxF = 1.05;        % one plus tolerance
stdout = 1;
 
if nargin < 7,
  maxiter = 50;       % maximum number of iterations
  dXmin = 1e-6;       % minimum change in X at each iteration
  dFmin = 1e-6;       % minimum change in F ad each iteration
  minRMSE = 1e-6;     % minimum RMSE of predicted pixels
else
  maxiter = OPT(1);       % maximum number of iterations
  dXmin = OPT(2);       % minimum change in X at each iteration
  dFmin = OPT(3);       % minimum change in F ad each iteration
  minRMSE = OPT(4);     % minimum RMSE of predicted pixels 
end

[n,m] = size(X0);
if size(Y,1) ~= n,
  error(' X0 and Y matrix dimension mismatch!')
end
p = size(Y,2);

% add the sum constraint as equation
Y = [Y; ones(1,p)];
X0 = [X0; ones(1,m)];

% Set initial solution
if nargin < 3,
  F0 = [];
elseif size(F0,2) == 1,         % if m-by-1
  F0 = F0(:,ones(1,p));         % m-by-p
end

if isempty(F0),
  F0 = inv(X0'*X0)*X0'*Y;       % m-by-p    
end

% Set variances/covariance matrices
if nargin < 6,
  COVF = 0.01;
  if nargin < 5,
    COVX = 0.01;
    if nargin < 4, % set data covariance
      %COVY = cov(Y');           % (n+1)-by-(n+1)
      COVY = 0.01;
    end
  end
end

if length(COVY) == 1, % zero covariances/constant variances
  COVY = diag(COVY(ones(1,n+1)));   % (n+1)-by-(n+1)
elseif min(size(COVY)) == 1, % variances specified/zero covariances
  COVY = diag(COVY);   % (n+1)-by-(n+1)  
end
COVY(n+1,n+1) = 0.01;   % set variance of sum2one constraint

if length(COVX) == 1, % zero covariances/constant variances
  COVX = diag(COVX(ones(1,m)));   % m-by-m
elseif min(size(COVX)) == 1, % variances specified/zero covariances
  COVX = diag(COVX);
end
if length(COVF) == 1, % zero covariances/constant variances
  COVF = diag(COVF(ones(1,m)));     % m-by-m
elseif min(size(COVF)) == 1, % variances specified/zero covariances
  COVF = diag(COVF);                % m-by-m    
end

msg = sprintf('Estimating %d parameters using %d data...',m*(n+1+p),p*(n+1));
wb = waitbar(0,msg);

% initialize variables for iterative update of parameters
Z0 = [X0' F0]; Z = Z0;                      % m by n+1+p
ICOVY = kron(sparse(eye(p)),inv(COVY));     % (n+1)*p by (n+1)*p

ICOVZ = [kron(sparse(eye(n+1)),inv(COVX)) sparse(m*(n+1),m*p); sparse(m*p,m*(n+1)) kron(sparse(eye(p)),inv(COVF))]; % m(n+1+p) by m(n+1+p) 
%ICOVZ = [kron(sparse(eye(m)),inv(COVX)) sparse(m*(n+1),m*p); sparse(m*p,m*(n+1)) kron(sparse(eye(p)),inv(COVF))]; % m(n+1+p) by m(n+1+p) 
k1 = 1:n+1; k2 = n+2:n+1+p; % column inidices to endmembers and fractions
niter = 0; RMSE = Inf; dX = Inf; dF = Inf;
while RMSE > minRMSE & niter < maxiter & dX > dXmin & dF > dFmin,
  % apply non-negativity constraint to fractions  
  Z(:,k2) = max(min(Z(:,k2), maxF),minF);
  
  % calculate the jacobian 
  J = permute(reshape(kron(Z(:,k2),eye(n+1)),[(n+1),m*(n+1),p]),[1 3 2]);
  J = [reshape(J,(n+1)*p,m*(n+1)) kron(sparse(eye(p)),Z(:,k1)')]; % p*(n+1) by m*(n+1+p) = #eqs by #pars
  % calculate model residuals
  R = Y-Z(:,k1)'*Z(:,k2);
  
  % calculate parameter increment
  dZ = reshape(inv(J'*ICOVY*J+ICOVZ)*(J'*ICOVY*R(:)-ICOVZ*(Z(:)-Z0(:))),m,n+1+p);
  Z(:) = Z+dZ;  % update parameters
  
  % update stop conditions
  dX = max(max(abs(dZ(:,k1)./Z(:,k1))));                % max relative error of endmember
  dF = max(max(abs(dZ(:,k2)./Z(:,k2))));                % max relative error of fraction
  RMSE = sqrt(mean2(R.^2));                             % root mean square error
  niter = niter + 1;                                    % number of interation 
  waitbar(niter/maxiter,wb);
  fprintf(stdout,'iter/maxiter = %d/%d,',niter,maxiter)
  fprintf(stdout,'RMSE/min RMSE = %f/%f,',RMSE,minRMSE)
  fprintf(stdout,'dX/min dX = %f/%f,',dX,dXmin)    
  fprintf(stdout,'dF/min dF = %f/%f\n',dF,dFmin)
end
close(wb)

X = Z(:,1:n)';
F = Z(:,k2);