function [a,P,resnorm] = msafit(M,Y,P,a,method,UBP)
% MSAFIT - Multiple Scattering Approximation (MSA) Model Fitting
%
% Synopsis:
%   [X,RES] = MSAFIT(EM,Y,P)
%   [X,RES] = MSAFIT(EM,Y,P,X0)
%   [X,P,RES] = MSAFIT(EM,Y,P0,X0,'local')
%   [X,P,RES] = MSAFIT(EM,Y,P0,X0,'global')
%   [X,P,RES] = MSAFIT(EM,Y,P0,X0,'global',BKGEM)
%
% Inputs:
%   EM - N-by-M matrix with M end-members of N spectral bands
%   Y - N-by-K matrix of P mixed pixels from a multispectral image of N spectral bands
%   P - M-by-M matrix or M-by-M-by-K array of recollision probability matrices (RPMs).
%       The sum of columns of each RPM must be between zero and one.
%   P0 - M-by-M matrix or M-by-M-by-K array of initial recollision
%       probabilities. Set P0 = [] to use a random initialization of the RPM.
%   X0 - M-by-1 or M-by-K matrix of K initial fraction vectors. Set X0 = []
%       to use a random initialization of the optimization
%   'local' - perform a per-pixel esimation of P starting at P0. 
%       For each mixed spectra the recollision probabilities are estimated
%       independently of other mixed spectra. For estimation to be robust
%       the number of spectral bands must be N >= M^2+M-1. 
%   'global' - performs a global estimation of P stating at P0. The recollision
%       probability is estimated iteratively and considering all mixed pixels.
%       The method is based on the average of local probabilities, so that
%       MSAFIT is called several times with the 'local' option. This method
%       may be time consumming.
%   BKGEM - optionally specify indices background endmembers which has zero
%   recollision probability.
%
% Outputs:
%   X - M-by-K matrix of K fraction vectors for the M end members
%   P - M-by-M-by-K array of recollision probability matrices for each mixed
%       pixel, where P(:,:,k) corresponds to the k-th pixel. If the flag 'global' 
%       is used, then K = 1, and the same P is used for all the pixels.
%   RES - Residual norm of the fitting, i.e., RES =NORM(MSA(EM,R,X)-Y).^2
%
% Note:
%   - This function requires the optimization toolbox
%   - A multiple scattering model is used (See MSA for details on the mixing model)

[n,m] = size(M);    % dimension of end-members
if size(Y,1) ~= n,
  error(' EM and Y matrix dimension mismatch!')
end
k = size(Y,2);      % number of pixels to unmix

%OPT = optimset('Display','off','LargeScale','off');
%OPT = optimset('Display','off','LargeScale','on','Algorithm','interior-point');
OPT = optimset('Display','off','LargeScale','on','Algorithm','active-set');

if nargin < 5,  % optimize for fractions only, P must be provided
  if nargin < 4,
    a = repmat(1/m,m,1);                     % initial fractional
  end
  if size(a,2) == 1, a = a(:,ones(1,k));  end
  if size(P,3) == 1,  P = P(:,:,ones(1,k)); end  
  
  Aeq = ones(1,m);
  Beq = 1;
  LB = zeros(m,1);
  UB = ones(m,1);
  %t0 = clock;
  wb = waitbar(0,'Running optimization rutine. Please wait...'); 
  for i = 1:k, % for each pixel
    [a(:,i),resnorm(i)] = fmincon(@obj1,a(:,i),[],[],Aeq,1,LB,UB,[],OPT,Y(:,i),M,P(:,:,i));
    waitbar(i/k,wb);
  end
  close(wb)
  %etime(clock,t0)
  % pass the output;
  if nargout > 1,  P = resnorm.^2; end
  return;
end

if strcmp(method,'local'), % optimizes both fractions and TPMs locally
  
  % Set initial solutions
  if isempty(a), 
    a = rand(m,1);                             % initial fractional
    a = a./sum(a);                            % fractions must sum to one            
  end
  if isempty(P),
    q = rand(m,1);                             % initial escaping probability
    P = rand(m,m);                             % initial re-collision probabilities
    P(:) = P.*repmat((1-q)./sum(P,2),1,m);     % make sure SUM(P,2) = 1-q
  end
  if size(a,2) == 1, a = a(:,ones(1,k));  end
  if size(P,3) == 1,  P = P(:,:,ones(1,k)); end

  % Set constraints
  A = convmtx(ones(1,m), m^2+1); A = A(1:m:end,:);
  B = ones(1,m);
  Aeq = A(m+1,:);
  A = A(1:m,:);
  LB = zeros(m,m+1);
  if nargin < 6, 
    UB = ones(m,m+1);
  else
    UB = [UBP' ones(m,1)];
  end 
  %t0 = clock;
  wb = waitbar(0,'Running optimization rutine. Please wait...');     
  for i = 1:k, % for each pixel    
    x = [P(:,:,i)' a(:,i)];
    [x,resnorm(i)] = fmincon(@obj2,x,A,B,Aeq,1,LB,UB,[],OPT,Y(:,i),M);
    P(:,:,i) = x(1:m,1:m)';
    a(:,i) = x(1:m,m+1);
    waitbar(i/k,wb);
  end
  close(wb)  
  %etime(clock,t0)
  %if nargout > 2, resnorm = resnorm.^2; end % return squared norm
  
elseif strcmp(R,'global'),  % optimizes both fractions and TPMs globally
    
  maxerrP = 1e-6;
  maxiter = 10;
  [a,P] = msafit(M,Y,'local'); 
  P = mean(P,3);  % First estimation
  errP = inf;
  niter = 0;
  while errP > maxerrP & niter < maxiter,
    P0 = P;
    a0 = msafit(M,Y,P0,a); % refine abundance
    [a,P] = msafit(M,Y,P0,a0,'local');
    P = mean(P,3);          % Next aproximation
    errP = norm(P-P0);
    niter = niter+1;    
  end
  resnorm = obj2([P a],Y,M).^2;
end

function f = obj1(x,Y,M,P)
% OBJ - Objective function (P fix)

[n,m] = size(M);    % dimension of end-members
f = norm(Y-msa(M,P,x));

function f = obj2(x,Y,M)
% OBJ - Objective function (P variable)

[n,m] = size(M);    % dimension of end-members
P = x(1:m,1:m)';
a = x(1:m,m+1);
f = norm(Y-msa(M,P,a));