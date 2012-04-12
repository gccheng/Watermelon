function [IX,IY,IT] = partial_derivative_3D(data, gsize)
% PARTIAL_DERIVATIVE_3D 

    if ndims(data)~=3
       error('data is not in NxMxP form'); 
    end
    
    gsigma = gsize/2.0; %gsize = 3; 
    [IX, IY, IT] = gaussgradient(data, gsize, gsigma);
    
end