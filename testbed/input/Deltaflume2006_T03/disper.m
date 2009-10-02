function k = disper(w, h, g);

% DISPER   Linear dispersion relation
%
% usage:   k  = disper(w,h,g)     or
%
%          k  = wave number             (2 * pi / wave length)
%          w  = wave angular frequency  (2 * pi / wave period)
%          h  = water depth
%          g  = gravitational acceleration constant, optional (DEFAULT 9.81)
%
%          absolute error in k*h < 5.0e-16 for all k*h
%

%          programmer: G. Klopman, Delft Hydraulics, 6 Dec 1994

if nargin < 3,
  g = 9.81;
end;

w2 = (w.^2) .* h ./ g;
q  = w2 ./ (1 - exp (-(w2.^(5/4)))) .^ (2/5);

for j=1:2,
  thq     = tanh(q);
  thq2    = 1 - thq.^2;
  a       = (1 - q .* thq) .* thq2;
  b       = thq + q .* thq2;
  c       = q .* thq - w2;
  arg     = (b.^2) - 4 .* a .* c;
  arg     = (-b + sqrt(arg)) ./ (2 * a);
  iq      = find (abs(a.*c) < 1.0e-8 * (b.^2));
  arg(iq) = - c(iq) ./ b(iq);
  q       = q + arg;
end;

k = sign(w) .* q ./ h;

ik    = isnan (k);
k(ik) = zeros(size(k(ik)));