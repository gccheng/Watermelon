function [CData] = convole_3D(data,szConv)
% CONVOLE_3D 3D convolution.
% data
%   The 3 dimensional data to be convolved
% func
%   The kernel function of convolution. There's only one implementation
%   now, i.e., 'Gauss'

    sigma = szConv/3.0;
    
%     % Generate 3-D Gaussian kernels along x/y/t directions;
%     halfsize = round(szConv/2);
%     G = zeros(szConv, szConv, szConv);
%     for i=1:szConv
%         for j=1:szConv
%             for k=1:szConv
%                 u=[i-halfsize-1 j-halfsize-1 k-halfsize-1];
%                 G(i,j,k)=gauss(u(2),sigma)*gauss(u(1),sigma)*gauss(u(3),sigma);
%             end
%         end
%     end
% 
%     % Added 2010-3-4
%     GH = G/sqrt(sum(sum(sum(abs(G).*abs(G)))));
% 
% %    load DGKernel;
% 
    GH = fspecial('gauss', [szConv,szConv], sigma);
    CData = imfilter(data,GH,'replicate','conv');
end
