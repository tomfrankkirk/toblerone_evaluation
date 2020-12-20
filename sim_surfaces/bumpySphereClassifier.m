function [ labels ] = bumpySphereClassifier(d1, d2, X, Y, Z) 


% Calculate the radii of the provided points. 
r = sqrt(X.^2 + Y.^2 + Z.^2); 

% Expected radii for the surfaces of the spheres
[rIn, rOut] = bumpySphereFunction(d1, d2, X, Y, Z); 

WMfltr = (r < rIn); 
GMfltr = ((r >= rIn) & (r < rOut));
CSFfltr = (r >= rOut);

labels = zeros(size(WMfltr));
labels(GMfltr) = 1;
labels(WMfltr) = 2;
labels(CSFfltr) = 3;

assert(sum(WMfltr(:)) + sum(GMfltr(:)) + sum(CSFfltr(:)) == numel(WMfltr), ... 
  'Error: not all voxels assigned a label.'); 

end 