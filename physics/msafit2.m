function [a,P,res] = msafit(M,Y)
% MSAFIT2 - MSA model fitting using a 2nd order approx.
%
% Synopsis:
%   [X,P,RES] = MSAFIT2(EM,Y)
%
% Inputs:
%   EM - N-by-2 matrix with 2 end-members of N spectral bands
%   Y - N-by-K matrix of P mixed pixels from a multispectral image of N
%   spectral bands
%
% Outputs:
%    X - M-by-K matrix of K fraction vectors for the M end members
%    P - M-by-M-by-K array of recollision probability matrices for each mixed
%        pixel, where P(:,:,k) corresponds to the k-th pixel. If the flag 'global' 
%        is used, then K = 1, and the same P is used for all the pixels.
%    RESNORM - Residual norm of the fitting, i.e., RES =
%    NORM(MSA(EM,R,X)-Y).^2

[n,m] = size(M);
if m ~= 2,
  error('Number of endmembers must be 2!') 
end
M = [M M(:,1).^2 M(:,1).*M(:,2) M(:,2).^2 zeros(n,1)];

k = size(Y,2); % number of mixed pixels

[X,res] = lsu(M,Y,'fullcon');

P11 = X(3,:)./X(1,:);
P22 = X(5,:)./X(2,:);
C(1,:) = X(4,:)+X(2,:).*(1-P22)+X(1,:).*(1-P11);
C(2,:) = X(2,:).*(X(2,:)-(1-P22))-X(1,:).*(X(1,:)+3*(1-P11))-2*X(4,:);
C(3,:) = X(4,:)+3*X(1,:).*(X(1,:)+(1-P11));
C(4,:) = -X(1,:).*(3*X(1,:)+1-P11);
C(5,:) = X(1,:).^2;

P = zeros(m,m,k);
a = repmat(nan,m,k);
myeps = 1e-6;
for i = 1:k,
  r = roots(C(:,i));
  r = real(r(r >= -myeps & r <= 1+myeps & abs(imag(r)) < myeps));
  if length(r) > 1,
    q = X(1:2,repmat(k,1,length(r)))./[r'; 1-r'];
    r = r(all(q >= -myeps & q <= 1+myeps,1));
  end
  if ~isempty(r),
    a(:,i) = [r(1); 1-r(1)];
  end
end
P(1,1,:) = P11;
P(2,1,:) = 1-P22-X(2,:)./a(2,:);
P(1,2,:) = 1-P11-X(1,:)./a(1,:);
P(2,2,:) = P22;
