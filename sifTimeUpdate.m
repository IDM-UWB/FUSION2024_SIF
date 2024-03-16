function [xp,Pp]=sifTimeUpdate(filtMean,filtVar,wMean,Q,Time,Nmin,Nmax,Eps,SIorder)
% time update (prediction) of SIF

% ----------------------------------------------------------------------
% Identification and Decision Making Group, Faculty of Applied Sciences,
% University of West Bohemia, Pilsen, Czech Republic

% https://idm.kky.zcu.cz/

% References:
% J. Dunik, O. Straka, M. Simandl, and E. Blasch: „Random-Point-Based Filters: Analysis and Comparison in Target Tracking,“ IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 2, pp. 1403-1421, 2015.
% J. Dunik, O. Straka, and M. Simandl: „Stochastic Integration Filter,“ IEEE Transactions on Automatic Control, vol. 58, no. 6, pp. 1561-1566, 2013.

% (c) 2021

% ----------------------------------------------------------------------
% inputs:
% - filtering state mean
% - filtering state covariance matrix
% - state noise mean
% - state noise covariance matrix
% - actual time instant
% - minimal number of iterations of stochastic integration rule (SIR)
% - maximal number of iterations of SIR
% - allowed threshold for integration error
% - order of SIR (orders 1, 3, 5 are currently supported)

% outputs:
% - predictive state mean
% - predictive state covariance matrix

% - initialisation
nx = size(filtMean,1);

Sp = chol(filtVar)';
Ix  = zeros(nx,1);
Vx  = zeros(nx,1);
IPx  = zeros(nx);
VPx  = zeros(nx);
N = 0; % number of iterations
% - SIR recursion for state predictive moments computation
% -- until either required number of iterations is reached or threshold is reached
while N<Nmin || all([N<Nmax, any([(norm(Vx)> Eps) (norm(VPx)> Eps)])])
    N = N+1;
    % -- cubature points and weights computation (for standard normal PDF)
    [SCRSigmaPoints,w] = stochasticCubatureRulePoints(nx,SIorder);
    % -- points transformation for given filtering mean and covariance
    %    matrix
    xpoints = bsxfun(@plus,Sp*SCRSigmaPoints,filtMean);
    % -- points transformation via dynamics  (deterministic part)
    if nx==2
        fpoints = ffunct(xpoints,wMean,Time);
    elseif nx==4
        fpoints = ffunct4(xpoints,wMean,Time);
    end
    % -- stochastic integration rule for predictive state mean and covariance
    %    matrix
    SumRx = fpoints*w';
    SumRPx=(w(ones(nx,1),:).*fpoints)*fpoints';
    % --- update mean Ix
    Dx = (SumRx-Ix)/N;
    Ix = Ix+Dx;
    Vx = (N-2)*Vx/N+Dx*Dx';
    % --- update covariance matrix IPx
    DPx = (SumRPx-IPx)/N;
    IPx = IPx+DPx;
    VPx = (N-2)*VPx/N+DPx*DPx';
end
xp = Ix;
Pp = IPx - Ix*Ix';
Pp = Pp + Q;
Pp = (Pp+Pp')/2;

end % function timeUpdate
