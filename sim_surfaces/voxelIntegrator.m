function [ truePVs ] = voxelIntegrator ( imageSize, voxSize, vox2world, ...
  supersample, func, voxList )

truePVs = zeros(prod(imageSize),3);
[I,J,K] = ind2sub(imageSize, 1:prod(imageSize)); 
voxIJK = [I', J', K']; 

for v = voxList 
  
  voxCent = voxIJK(v,:); 
  sX = ((1 / (2*supersample)) : 1 / supersample : 1 - (1 / (2*supersample))) - 0.5;
  [sX, sY, sZ] = ndgrid(sX, sX, sX);
  subVoxCents = [sX(:), sY(:), sZ(:)] + voxCent;
  
  % Map into world space
  subVoxCents(:,4) = 2; 
  svC = vox2world * (subVoxCents' - 1); 
  svC = svC(1:3,:)';
  
  svX = reshape(svC(:,1), supersample, supersample, supersample);
  svY = reshape(svC(:,2), supersample, supersample, supersample) ;
  svZ = reshape(svC(:,3), supersample, supersample, supersample) ;
  
  labels = func(svX, svY, svZ); 
  thisVoxVols = [0 0 0]; 
  thisVoxVols(1) = sum(labels(:) == 1);
  thisVoxVols(2) = sum(labels(:) == 2);
  thisVoxVols(3) = sum(labels(:) == 3);
  truePVs(v,:) = thisVoxVols / numel(labels);
  
end

end
