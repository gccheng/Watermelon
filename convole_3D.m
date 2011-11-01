function [CData] = convole_3D(data,szConv,varargin)
% CONVOLE_3D 3D convolution.
% data
%   The 3 dimensional data to be convolved
    
if (3==nargin)
    GH = varargin{1};
else
    %load DGKernel;
    halfsize = double(round(szConv/2.0));
    sigma = double(szConv/3.0);
    X = double(repmat(repmat((1:szConv)',[1,szConv]),[1,1,szConv]))-ones(szConv,szConv,szConv)*halfsize;
    Y = permute(X,[2 1 3]);
    T = permute(X,[3 2 1]);
    GH = 1.0/((2*pi)^(3.0/2.0)*sigma^3)*exp(-1.0/(2.0*sigma^2)*(X.^2+Y.^2+T.^2));
end

CData = imfilter(data,GH,'replicate','conv');
    
end
