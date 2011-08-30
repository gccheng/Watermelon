function [CData] = convole_3D(data,func,varargin)
% CONVOLE_3D 3D convolution.
% data
%   The 3 dimensional data to be convolved
% func
%   The kernel function of convolution. There's only one implementation
%   now, i.e., 'Gauss'

    if (strcmp(func,'Gauss')) && (numel(varargin)~=2)
        error('Gaussian kernal need two parameters: size and sigma!');
    end

    % Generate 3-D Gaussian kernels along x/y/t directions
    sigma = varargin{2};
    size = varargin{1};
    halfsize = round(size/2);
    G = zeros(size, size, size);
    for i=1:size
        for j=1:size
            for k=1:size
                u=[i-halfsize-1 j-halfsize-1 k-halfsize-1];
                G(i,j,k)=gauss(u(2),sigma)*gauss(u(1),sigma)*gauss(u(3),sigma);
            end
        end
    end

    % Added 2010-3-4
    G = G/sqrt(sum(sum(sum(abs(G).*abs(G)))));
    
    % 3D filtering
    CData = imfilter(data,G,'replicate','conv');
end
