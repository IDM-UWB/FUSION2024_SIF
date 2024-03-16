function [xp,Pp] = ukfp(xf,Pf,alpha,beta,kappa,Q,k)
% UKF time update (for cyclus form)
% [xp,Pp]=ukfp(xf,Pf,kappa,Q,k)

% sigma points calculation
[chi,wm,wc] = msp(xf,Pf,alpha,beta,kappa);
[a,b]=size(chi);

% predictive sigma points
chip = zeros(a,b);
for j=1:1:b      
    chip(:,j) = ffunct4(chi(:,j),zeros(a,1),k);
end
   
% predictive mean
xp = zeros(a,1);
for j=1:1:b
  xp = xp + wm(j)*chip(:,j);
end

% predictive covariance matrix
Pp = zeros(a);
for j=1:1:b
  Pp = Pp + wc(j)*(chip(:,j)-xp)*(chip(:,j)-xp)';
end

Pp = Pp + Q;