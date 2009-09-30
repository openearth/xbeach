function tmm10 = tref(tp,fnyq,gam)

% compute Tm-1,0 wave period from peak period

df = 0.005;
f = df:df:fnyq;
fp = 1/tp;

% JONSWAP  unscaled JONSWAP spectrum
%
% usage: y = jonswap(x, gam)
%
% x   : nondimensional frequency, divided by the peak frequency
% gam : peak enhancement factor, optional parameter (DEFAULT 3.3)
% y   : nondimensional relative spectral density, equal to one at the peak
%
xa  = f./fp;
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

tmm10 = sum(y./f)/sum(y)
