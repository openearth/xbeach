function [s]=round2(r,n)
s=num2str(round(r*10^n)/10^n);