function y = gauss(x,sigma)
% GAUSS 1D Gaussian function
    y = exp(-x^2/(2*sigma^2)) / (sigma*sqrt(2*pi));
end