function [rIn, rOut] = bumpySphereFunction(d1, d2, X, Y, Z)
%BUMPYSPHEREFUNCTION Summary of this function goes here
%   Detailed explanation goes here

% Convert the grid of XYZ points into spherical polars
[theta, phi, ~] = cart2sph(X,Y,Z); 

% Properties of the spherical surface
u = 5; 
lim = pi/2 * ( (u-1)/(u) ); 
tomodify = (abs(phi) < lim); 
s1 = sin(u*(phi - theta)).^20; 
s2 = sin(u*(phi + theta)).^20; 

% Each theta,phi pair maps to a unique r1,r2 pair that define the radii of the
% spherical surfaces
rIn = d1 * (1 - 0.1 * max(s1, s2)); 
rIn(~tomodify) = d1; 

rOut = d2 * rIn; 

end

