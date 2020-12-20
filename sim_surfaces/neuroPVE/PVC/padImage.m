function [ img ] = padImage( img )
%PADIMAGE Summary of this function goes here
%   Detailed explanation goes here

imgSize = size(img); 
imgSize = imgSize(1:3); 
[dSize, maxD] = max(imgSize);
dims = [1,2,3];
dims = dims(dims~=maxD);
for d = dims
  padAmount = dSize - imgSize(d);
  if d == 1
    img(end:end+padAmount,:,:) = 0;
  elseif d ==2
    img(:,end:end+padAmount,:) = 0;
  elseif d == 3
    img(:,:,end:end+padAmount) = 0;
  end
end
    
end



