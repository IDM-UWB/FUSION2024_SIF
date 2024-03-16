function zwrapped = bearingWrapping(z)
% Stochastic Interation Filter (SIF) - part of implementation illustration
% wrapping of bearings measurement

zwrapped = mod(z+pi,2*pi)-pi;