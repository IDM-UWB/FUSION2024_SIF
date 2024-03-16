function [xs, Ps] = sifRTSsmoothingUpdate(filtMean, filtVar, predMean, predVar, smoothMean_prev, smoothVar_prev, wMean, Time, Nmin, Nmax, Eps, SIorder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% - initialisation
nx = size(predMean,1);

efMean = filtMean;
efVar = filtVar;
Sf = chol(efVar)';
% Ix  = zeros(nx,1);
% Vx  = zeros(nx,1);
% IPx  = zeros(nx);
% VPx  = zeros(nx);
IPxx = zeros(nx,nx);
VPxx  = zeros(nx,nx);
N = 0; % number of iterations
% - SIR recursion for measurement predictive moments computation
% -- until either required number of iterations is reached or threshold is reached
while N<Nmin || all([N<Nmax, (norm(VPxx)> Eps)])
    N = N+1;
    % -- cubature points and weights computation (for standard normal PDF)
    [SCRSigmaPoints,w] = stochasticCubatureRulePoints(nx,SIorder);
    % -- points transformation for given filtering mean and covariance
    %    matrix
    xpoints = bsxfun(@plus,Sf*SCRSigmaPoints,efMean);
    % -- points transformation via measurement equation (deterministic part)
    fpoints = ffunct(xpoints,wMean,Time);
    % -- stochastic integration rule for predictive measurement mean and covariance
    %    matrix and predictive state and measurement covariance matrix
    % SumRx = fpoints*w';
    % SumRPx = fpoints*(fpoints.*w(ones(nx,1),:))';
    SumRPxx = xpoints(1:nx,:)*(fpoints.*w(ones(nx,1),:))';
    % % --- update mean Ix
    % Dx = (SumRx-Ix)/N;
    % Ix = Ix+Dx;
    % Vx = (N-2)*Vx/N+Dx.^2;
    % % --- update covariance matrix IPx
    % DPx = (SumRPx-IPx)/N;
    % IPx = IPz+DPx;
    % VPx = (N-2)*VPx/N+DPx.^2;
    % --- update cross-covariance matrix IPxx
    DPxx = (SumRPxx-IPxx)/N;
    IPxx = IPxx+DPxx;
    VPxx = (N-2)*VPxx/N+DPxx.^2;
end
% - measurement predictive moments
% xps = Ix;
Pxxps = IPxx - filtMean*predMean'; 
% - smoother gain
Ks = Pxxps/predVar;
% Ks = filtVar*diag([0.9, 1])'*inv(predVar);
% - smoothing moments computation
smoothMean = filtMean + Ks*(smoothMean_prev - predMean);
smoothVar = filtVar - Ks*(predVar - smoothVar_prev)*Ks';
smoothVar = (smoothVar+smoothVar')/2;
% - check for negative eigenvalues in the covariance matrix (from numerical reasons)
mieig = min(eig(smoothVar));
if mieig <= 0
    smoothVar =  smoothVar - mieig*eye(nx) - realmin;
end


%-------------
xs = smoothMean;
Ps = smoothVar;
end

