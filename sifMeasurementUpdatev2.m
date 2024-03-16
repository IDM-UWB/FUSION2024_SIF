function [xf,Pf]=sifMeasurementUpdatev2(predMean,predVar,vMean,R,Time,Measurement,Nmin,Nmax,Eps,SIorder)
% Stochastic Interation Filter (SIF) - part of implementation illustration
% measurement update (filtering) of SIF

% ----------------------------------------------------------------------
% Identification and Decision Making Group, Faculty of Applied Sciences,
% University of West Bohemia, Pilsen, Czech Republic

% https://idm.kky.zcu.cz/

% References:
% J. Dunik, O. Straka, M. Simandl, and E. Blasch: „Random-Point-Based Filters: Analysis and Comparison in Target Tracking,“ IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 2, pp. 1403-1421, 2015.
% J. Dunik, O. Straka, and M. Simandl: „Stochastic Integration Filter,“ IEEE Transactions on Automatic Control, vol. 58, no. 6, pp. 1561-1566, 2013.

% (c) 2024

% ----------------------------------------------------------------------

% inputs:
% - predictive state mean
% - predictive state covariance matrix
% - measurement noise mean
% - measurement noise covariance matrix
% - actual time instant
% - minimal number of iterations of stochastic integration rule (SIR)
% - maximal number of iterations of SIR
% - allowed threshold for integration error
% - order of SIR (orders 1, 3, 5 are currently supported)

% outputs:
% - filtering state mean
% - filtering state covariance matrix

% - initialisation
nx = size(predMean,1);
nz = size(Measurement,1);

epMean = predMean;
epVar = predVar;
Sp = chol(epVar)';
Iz  = zeros(nz,1);
Vz  = zeros(nz,1);
IPz  = zeros(nz);
VPz  = zeros(nz);
IPxz = zeros(nx,nz);
VPxz  = zeros(nx,nz);
N = 0; % number of iterations
% - SIR recursion for measurement predictive moments computation (mean in
% one recurstion, covariance matrices in other recursion)
% -- until either required number of iterations is reached or threshold is reached
while N<Nmin || all([N<Nmax, norm(Vz)> Eps])
    N = N+1;
    % -- cubature points and weights computation (for standard normal PDF)
    [SCRSigmaPoints,w] = stochasticCubatureRulePoints(nx,SIorder);
    % -- points transformation for given filtering mean and covariance
    %    matrix
    xpoints = bsxfun(@plus,Sp*SCRSigmaPoints,epMean);
    % xpoints = bsxfun(@plus,SCRSigmaPoints,0*epMean);
    % -- points transformation via measurement equation (deterministic part)
    hpoints = hfunct4(xpoints,vMean,Time);
    hpoints(1,:) = bearingWrapping(hpoints(1,:));
    % -- stochastic integration rule for predictive measurement mean and covariance
    %    matrix and predictive state and measurement covariance matrix
    SumRz = hpoints*w';
    % --- update mean Iz
    Dz = (SumRz-Iz)/N;
    Iz = Iz+Dz;
    Vz = (N-2)*Vz/N+Dz*Dz';
end
zp = Iz;
N = 0; % number of iterations
while N<Nmin || all([N<Nmax, any([(norm(VPz)> Eps) (norm(VPxz)> Eps)])])
    N = N+1;
    % -- cubature points and weights computation (for standard normal PDF)
    [SCRSigmaPoints,w] = stochasticCubatureRulePoints(nx,SIorder);
    % -- points transformation for given filtering mean and covariance
    %    matrix
    xpoints = bsxfun(@plus,Sp*SCRSigmaPoints,epMean);
    % -- points transformation via measurement equation (deterministic part)
    hpoints = hfunct4(xpoints,vMean,Time);
    hpoints(1,:) = bearingWrapping(hpoints(1,:));
    % -- stochastic integration rule for predictive measurement mean and covariance
    %    matrix and predictive state and measurement covariance matrix
    SumRPz = (hpoints-zp)*((hpoints-zp).*w(ones(nz,1),:))';
    SumRPxz = (xpoints(1:nx,:)-predMean)*((hpoints-zp).*w(ones(nz,1),:))';
    % --- update covariance matrix IPz
    DPz = (SumRPz-IPz)/N;
    IPz = IPz+DPz;
    VPz = (N-2)*VPz/N+DPz*DPz;
    % --- update cross-covariance matrix IPxz
    DPxz = (SumRPxz-IPxz)/N;
    IPxz = IPxz+DPxz;
    VPxz = (N-2)*VPxz/N+DPxz.^2;
end
% - measurement predictive moments
Pzp = IPz + R + Vz;
Pxzp = IPxz;
% - filter gain
K = Pxzp/Pzp;
% - filtering moments computation
filtMean = predMean+ K*(Measurement-zp);
Happrox = Pxzp'/predVar;
filtVar = (eye(nx)-K*Happrox)*predVar*(eye(nx)-K*Happrox)'+K*R*K';

xf = filtMean;
Pf = filtVar;

end %function measurement update