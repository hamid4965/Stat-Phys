function [R,C,Ic] = msa(w,P,a,j)
% MSA - Compute the multiple scatering albedo from several end members
%
% Synopsis:
%   R = MSA(EM,P,F)
%   [R,CRB] = MSA(EM,P,F,J)
%   [R,CRB,IC] = MSA(EM,P,F,J)
%
% Inputs:
%   EM - N-by-M matrix with M endmembers with N spectral bands. Endmembers
%       are spectral reflectance (or scattering) coefficients.
%   P - M-by-M matrix (or M-by-M-by-K array) of recollision probabilities (RPM). 
%       P(i,j) specifies the recollision probability from the i-th
%       end-member to the j-th end-member. Then if P is diagonal, materials
%       are assumed in segregated patters so that light does not interact
%       with two or more end-members before being reflected away. 
%       SUM(P,2) <= 1. K is the number of signatures (or pixels) to be modelled.
%   F - M-by-1 vector (or M-by-K matrix) with fraction of interecepted light 
%       by each end-member. Note that SIZE(F,2) == SIZE(P,3).
%   J - optionally pass indices of canopy endmembers EM(:,J) to use in
%       canopy absortance computation. If only one mixed signature is specified
%       (K=1), you may specify several canopy endmembers (LENGTH(J) >= 1).
%       However, for many signatures (K > 1) you have to specify one canopy
%       endmember per signature (LENGTH(J)==K). 
%
% Outputs:
%  R - N-by-K matrix of above-canopy reflectance. For each column in F,
%    there is a column in R. A = 1-R is  the matrix of above-canopy
%    absortance. 
%  CRB - N-by-3-by-K matrix of canopy reflectance-absortance-transmitance. For each column in F,
%    there is a column in T. Canopy absortance is the portion of light
%    that is absorbed by canopy endmembers.
%  IC - optionally returns the interaction coefficients, i.e., canopy
%    absortance contributions to endmember absortances. 
%
% Notes: 
%    - The model used here considers constant recollision probabilities from
%    one end-member to other end-members, which does not happen in reality.
%    - Lambertian surfaces are considered, so that the albedo is
%    is independent of direction of incoming light.
%    The albedo for a particular wavelength can be expressed:
%      Ri = q'*inv(I-wi*P)*wi*A
%    where 
%    - q(i) is the probability of scaping recollision from the i-th
%    endmember, and it is given by
%      q = 1-sum(P,2)
%    - wi is the i-end-member albedo (reflectacn + transmitance), i.e., 
%      wi = diag(EM(i,:))
%
% Example 1: This example shows that when the RPM is constant, the mixing
%   is esencially linear, but with transformed end-members.
% 
%   p = 10000;
%   n = 2;
%   m = 3;
%   E = rand(n,m);                       % endmembers
%   F = rand(m,p);                          % fractional conver
%   F(:) = F./repmat(sum(F,1),m,1);         % fractions must sum to one
%   P = rand(m,m);                          % RPM
%   c = rand(m,1)./sum(P,2);
%   P(:) = P.*c(:,ones(1,m));               % ensure SUM(P,2) <= 1
%   X = msa(E,P,F);                         % mixed pixels
%   Et = msa(E,P,eye(m));                   % transformed endmembers
%   H = convhull(Et(1,:),Et(2,:));
%   plot(X(1,:),X(2,:),'.',E(1,:),E(2,:),'o',Et(1,H(:,1)),Et(2,H(:,1)),'k-*','MarkerFaceColor','k')
%   xlabel('band 1')
%   ylabel('band 2')
%   legend({'Mixed pixel','Original endmember','Transformed endmember'})
%   title('MSA mixture with constant RPM')
%
% Example 2: This example shows non-linear mixing by alowing the RPM vary
% for various pixels
%   p = 10000;
%   n = 2;
%   m = 2;
%   E = rand(n,m);                          % endmembers
%   F = rand(m,p);                          % fractional cover
%   F(:) = F./repmat(sum(F,1),m,1);         % fractions must sum to one
%   P = rand(m,m,p);                        % RPM
%   c = rand(m,1,p)./sum(P,2);
%   P(:) = P.*c(:,ones(1,m),:);             % ensure SUM(P,2) <= 1
%   X = msa(E,P,F);                         % mixed pixels
%   E0 = [E monomial(E,2) monomial(E,3) monomial(E,4) zeros(n,1)]';    % and infinity order end member
%   H = convhull(E0(:,1),E0(:,2));
%   plot(X(1,:),X(2,:),'.',E(1,:),E(2,:),'o',E0(H(:,1),1),E0(H(:,1),2),'k-+','MarkerFaceColor','k')
%   xlabel('band 1')
%   ylabel('band 2')
%   legend({'Mixed pixel','Original endmember','Higher-order endmembers'})
%   title('MSA mixture with variable RPM')
%
% Example 3: This example shows iso-fraction trajectories with
% variaing scapping probabilities. Both endmembers have the same scapping
% probability. Proportions of intra- vs inter-recollisions varies between
% trajectories as well.
%
% p = 101;                             % points in each traectory
% m = 2;                               % number of endmembers
% E = [.2 .9; .9 .5];                       % endmembers
% P = zeros(m,m,p);                    % RPM
% q = [0:p-1]/(p-1);                   % scapping probabilities
% k = 1;
% x1 = zeros(p,25);
% x2 = zeros(p,25);
% for  i = 0:4,
%   F = [i 4-i]'/4;                        % fractional conver
%   for  j = 0:4,
%     P(1,1,:) = (1-q)*(4-j)/4;
%     P(1,2,:) = (1-q)*j/4;
%     P(2,2,:) = (1-q)*(4-j)/4;
%     P(2,1,:) = (1-q)*j/4;
%     X = msa(E,P,F(:,ones(1,p)));           % mixed pixels
%     x1(:,k) = X(1,:)';
%     x2(:,k) = X(2,:)';
%     k = k+1;
%   end
% end
%  plot(x1(:,3:5:end),x2(:,3:5:end),'-','LineWidth',2); hold on
%  text(x1(end,3:5:end)+.01,x2(end,3:5:end),{'\alpha=0','\alpha=0.25','\alpha=0.5','\alpha=0.75','\alpha=_1'},'FontSize',11,'FontName','TimesNewRoman') 
%  plot(x1,x2,'k:',E(1,:),E(2,:),'ko','MarkerFaceColor','k')
%  text(E(1,:)+.01,E(2,:)+.01,{'X_1','X_2'},'FontSize',11,'FontName','TimesNewRoman','FontWeight','Bold')
%  plot(x1([0 25 50 75 100]+1,3:5:end)',x2([0 25 50 75 100]+1,3:5:end)','k-');
%  text(x1([0 25 50 75 100]+1,3)'+.01,x2([0 25 50 75 100]+1,3)'-0.01,{'q=0.0','q=0.25','q=0.50','q=0.75','q=1.0'},'FontSize',9,'FontName','TimesNewRoman','Rotation',atan2(diff(E(2,:)),diff(E(1,:)))*180/pi)
%  plot(reshape(x1(51,:),5,5),reshape(x2(51,:),5,5),'k.-')
%  text(x1(51,20+[1 3 5])+.01,x2(51,20+[1 3 5]),{'p_{1,1}=q','p_{1,1}=q/2','p_{1,1}=0'},'FontSize',9,'FontName','TimesNewRoman','Rotation',atan2(diff(x2(51,[21 22])),diff(x1(51,[21 22])))*180/pi+90,'HorizontalAlignment','left');
%  set(gca,'xlim',[0 1],'ylim',[0,1])
%   xlabel('band 1')
%   ylabel('band 2')
%   title('MSA mixing space structure')

Cflag = nargout > 1 && nargin > 3;
len = size(a,2); % number of pixels
[n,m] = size(w);      % number of end members;

if size(P,3) == 1, % constant RPM
  P = P(:,:,ones(len,1)); 
end

q = squeeze(1-sum(P,2));
I = eye(m);
R = zeros(n,len);
Ic = zeros(m,n,len);
C = [];
if Cflag,   
  C = zeros(n,3,len); % canopy radiation budget
  qC = zeros(size(q));  qC(j,:) = q(j,:);
  aC = zeros(size(w));  aC(:,j) = 1-w(:,j);
end

for k = 1:len, % for each pixel
  for i = 1:n, % for each wavelength
    Ic(:,i,k) = (I-P(:,:,k)'*diag(w(i,:)))\a(:,k);
    R(i,k) = q(:,k)'*diag(w(i,:))*Ic(:,i,k);
    if Cflag, % optionally compute canopy radiation budget
      C(i,1,k) = qC(:,k)'*diag(w(i,:))*Ic(:,i,k);
      C(i,2,k)= aC(i,:)*Ic(:,i,k);
      C(i,3,k) = 1-C(i,1,k)-C(i,2,k);
    end
  end
end
%Ic = Ic';
