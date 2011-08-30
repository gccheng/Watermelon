function [IX,IY,IT] = partial_derivative_3D(data)
% PARTIAL_DERIVATIVE_3D 

    if ndims(data)~=3
       error('data is not in NxMxP form'); 
    end
    if nargout~=3
       error('incorrect number of output parameters'); 
    end
    
    gsize = 5; gsigma = 1.0;
    [IX, IY, IT] = gaussgradient(data, gsize, gsigma);
    
end