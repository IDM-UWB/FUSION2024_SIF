function [xpred] = ffunct4(x,w,k)
% state dynamics

xpred = [1 1 0 0; 0 1 0 0; 0 0 1 1; 0 0 0 1]*x + w;
