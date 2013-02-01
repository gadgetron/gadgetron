function [img] = ismrm_transform_kspace_to_image(k, dim)
%
%  [img] = ismrm_transform_kspace_to_image(k, dim)
%
%  Fourier transform from k-space to image space along a given or all 
%  dimensions
%
%  INPUT:
%    - k       [kx,ky,..]    : k-space data
%    - dim     vector        : Vector with dimensions to transform
%
%  OUPUT:
%    - img    [x,y,...]      : Data in image space (along transformed
%                                                   dimensions)
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
    img = fftshift(ifftn(ifftshift(k))) .* sqrt(numel(k));
else
   img = k;
   for d=1:length(dim),
      img = fftshift(ifft(ifftshift(img,dim(d)),[],dim(d)),dim(d)) .* sqrt(size(img,d)); 
   end
end

return