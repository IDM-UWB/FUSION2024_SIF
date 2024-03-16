function [z] = hfunct4(x,v,k)
% measurement equation

sx = 50;
sy = 0;

z = [atan2((x(3,:)-sy),(x(1,:)-sx)); sqrt( (x(1,:)-sx).^2+(x(3,:)-sy).^2)] + v;

