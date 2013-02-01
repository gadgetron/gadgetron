function [k] = ismrm_transform_image_to_kspace(img, dim)
%
%  [k] = ismrm_transform_image_to_kspace(img, dim)
%
%  Fourier transform from image space to k-space space along a given or all 
%  dimensions
%
%  INPUT:
%    - img     [x,y,..]      : image space data
%    - dim     vector        : Vector with dimensions to transform
%
%  OUPUT:
%    - k       [kx,ky,...]   : Data in k-space (along transformed dimensions)
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%   Philip Beatty (philip.beatty@sri.utoronto.ca)
%

if nargin < 2,
    dim = [];
end    
   
if isempty(dim),
    k = fftshift(fftn(ifftshift(img))) ./ sqrt(numel(img));
else
   k = img;
   for d=1:length(dim),
      k = fftshift(fft(ifftshift(k,dim(d)),[],dim(d)),dim(d)) ./ sqrt(size(k,d)); 
   end
end

return