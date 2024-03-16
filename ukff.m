function [xf,Pf] = ukff(xp,Pp,z,alpha,beta,kappa,R,k)
% UKF measurement update (for cyclus form)
% [xf,Pf]=ukff(xp,Pp,z,kappa,R,k)

% sigma points calculation
[chip, wm, wc] = msp(xp,Pp,alpha,beta,kappa);
[a,b] = size(chip);
c = size(z,1);


% measurement predictive s-points
dzeta = zeros(c,b);
for j=1:1:b      
    dzetap(:,j) = hfunct4(chip(:,j),zeros(c,1),k);
end
dzetap(1,:) = bearingWrapping(dzetap(1,:));

% measurement prediciton
zp = zeros(c,1);
for j=1:1:b
  zp = zp + wm(j)*dzetap(:,j);
end

% measurement prediction cov. matrix
Pzp = zeros(c);
for j=1:1:b
  Pzp = Pzp + wc(j)*(dzetap(:,j)-zp)*(dzetap(:,j)-zp)';
end
    
Pep = Pzp + R;    

% state and measurement prediction cov. matrix
Pxzp = zeros(a,c);
for j=1:1:b
  Pxzp = Pxzp + wc(j)*(chip(:,j)-xp)*(dzetap(:,j)-zp)';
end

%filtering update
K = Pxzp*inv(Pep);
e = z - zp;
xf = xp + K*e;
% Pf = Pp - K*Pep*K';
Happrox = Pxzp'/Pp;
Pf = (eye(a)-K*Happrox)*Pp*(eye(a)-K*Happrox)'+K*R*K';
   
% %weight/likelihood
% psi = 1/sqrt(2*pi*(Pep))*exp(-0.5*((z-zp)^2)/(Pep));