function y = jonswap (x, gam)

% JONSWAP  unscaled JONSWAP spectrum
%
% usage: y = jonswap(x, gam)
%
% x   : nondimensional frequency, divided by the peak frequency
% gam : peak enhancement factor, optional parameter (DEFAULT 3.3)
% y   : nondimensional relative spectral density, equal to one at the peak
%

if nargin < 2,
  gam = 3.3;
end;

xa  = abs(x);
zxa = find (xa == 0);
if length(zxa) > 0,
  xa(zxa) = 1e-20 * ones(size(xa(zxa)));
end;

sigma = (xa < 1) * 0.07 + (xa >= 1) * 0.09;
fac1  = xa .^ (-5);
fac2  = exp (-1.25*(xa.^(-4)));
fac3  = (gam .* ones(size(xa))) .^ (exp (-((xa-1).^2) ./ (2*sigma.^2)));

y = fac1 .* fac2 .* fac3;

y = y / max(y);

%
