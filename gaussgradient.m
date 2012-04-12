function [gx,gy,gt]=gaussgradient(IM,size,sigma)
% GAUSSGRADIENT Gradient using first order derivative of Gaussian.
%  [gx,gy,gt]=gaussgradient(IM,sigma) outputs the gradient image gx,gy and
%  gt image IM using a 3D Gaussian kernel. Sigma is the standard deviation of
%  this kernel along all the three directions.
% IM
%   3D image cube
% size
%   The size of the Gaussian kernel
% sigma
%   The variance of Gaussian function

    % Generate 3-D Gaussian kernels along x/y/t directions
    halfsize = double(round(size/2));
    sigma = double(sigma);

%     X = double(repmat(repmat((1:size)',[1,size]),[1,1,size]))-ones(size,size,size)*halfsize;
%     Y = permute(X,[2 1 3]);
%     T = permute(X,[3 2 1]);
%     hx = -exp(-1.0/(2.0*sigma^2)*(X.^2+Y.^2+T.^2)).*X;
%     hy = -exp(-1.0/(2.0*sigma^2)*(X.^2+Y.^2+T.^2)).*Y;
%     ht = -exp(-1.0/(2.0*sigma^2)*(X.^2+Y.^2+T.^2)).*T;

    X = double(repmat(repmat((1:size)',[1,size]),[1,1,size]))-ones(size,size,size)*halfsize;
    hx = -exp(-1.0/(2.0*sigma^2)*(X.^2)).*X;
    hy = permute(hx,[2 1 3]);
    ht = permute(hx,[3 2 1]);
    
%     hSobel = fspecial('sobel');
%     hx = repmat(hSobel, [1 1 3]);
%     hy = repmat(hSobel', [1 1 3]);
%     ht = permute(hx, [3 2 1]);

    % 3D filtering
    gx = imfilter(IM,hx,'replicate','conv');
    gy = imfilter(IM,hy,'replicate','conv');
    gt = imfilter(IM,ht,'replicate','conv');

end