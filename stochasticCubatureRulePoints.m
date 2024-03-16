function [SCRSigmaPoints,weights] = stochasticCubatureRulePoints(nx,order)
% computation of cubature points and weights for the stochastic integration
% rule of order 1,3,5,

% ----------------------------------------------------------------------
% Identification and Decision Making Group, Faculty of Applied Sciences,
% University of West Bohemia, Pilsen, Czech Republic

% https://idm.kky.zcu.cz/

% References:
% J. Dunik, O. Straka, M. Simandl, and E. Blasch: „Random-Point-Based Filters: Analysis and Comparison in Target Tracking,“ IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 2, pp. 1403-1421, 2015.
% J. Dunik, O. Straka, and M. Simandl: „Stochastic Integration Filter,“ IEEE Transactions on Automatic Control, vol. 58, no. 6, pp. 1561-1566, 2013.

% (c) 2021

% ----------------------------------------------------------------------

switch order
    case 1
        X = randn(nx,1);
        SCRSigmaPoints = [X,-X];
        weights = [0.5 0.5];
    case 3
        CRSigmaPoints = [zeros(nx,1),eye(nx),-eye(nx)];
        rho = sqrt(chi2rnd(nx+2));
        Q = RandOrthMat(nx);
        SCRSigmaPoints = Q*rho*CRSigmaPoints;
        weights=[1-nx/rho^2,0.5*ones(1,2*nx)/rho^2];
    case 5
        % generating random values
        r=sqrt(chi2rnd(2*nx+7));
        q=betarnd(nx+2,3/2);
        rho=r*sin(asin(q)/2);
        delta=r*cos(asin(q)/2);
        
        % calculating weights
        c1up=nx+2-delta^2;
        c1do=rho^2*(rho^2-delta^2);
        c2up=nx+2-rho^2;
        c2do=delta^2*(delta^2-rho^2);
        cdo=2*(nx+1)^2*(nx+2);
        c3=(7-nx)*nx^2;
        c4=4*(nx-1)^2;
        coef1=c1up*c3/cdo/c1do;
        coef2=c2up*c3/cdo/c2do;
        coef3=c1up*c4/cdo/c1do;
        coef4=c2up*c4/cdo/c2do;
        
        weights=[1-nx*(rho^2+delta^2-nx-2)/(rho^2*delta^2),ones(1,2*nx+2)*coef1,ones(1,2*nx+2)*coef2,ones(1,nx*(nx+1))*coef3,ones(1,nx*(nx+1))*coef4];
        
        %calculating sigma points
        Q = RandOrthMat(nx);
        v = zeros(nx,nx+1);
        for i=1:nx
            v(i,i)=sqrt((nx+1)*(nx-i+1)/nx/(nx-i+2));
            for j=i+1:nx+1
                v(i,j)=-sqrt((nx+1)/((nx-i+1)*nx*(nx-i+2)));
            end
        end
        v = Q*v;
        
        y=zeros(nx,nx*(nx+1)/2);
        cnt=0;
        for j=1:nx+1
            for i=1:j-1
                cnt=cnt+1;
                y(:,cnt)=(v(:,j)+v(:,i))/norm(v(:,j)+v(:,i),2);
            end
        end
        
        SCRSigmaPoints=[zeros(nx,1),-rho*v,rho*v,-delta*v,+delta*v,-rho*y,rho*y,-delta*y,delta*y];
    otherwise
        disp('We are working hard, but ... this rule order is currently missing.')
end
end
