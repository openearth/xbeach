function [vsed,vero,R]=profana(x,z0,z,level);
dz=z-z0;
[xi,zi]=intersections(x,z0,x,z);
R=min(xi(zi>level));
Q=max(xi(zi<level));
P=max(xi(xi<Q));
iR=find(x<R,1,'last');
iQ=find(x<Q,1,'last');
iP=find(x<P,1,'last');
% accretion area
vsed=.5*(x(iP+1)-P)*dz(iP+1);
vsed=vsed+.5*sum(( x(iP+2:iQ)- x(iP+1:iQ-1)).* ...
                 (dz(iP+2:iQ)+dz(iP+1:iQ-1)));
vsed=vsed+.5*(Q-x(iQ))*dz(iQ);
% erosion area
vero=.5*(x(iQ+1)-Q)*dz(iQ+1);
vero=vero+.5*sum(( x(iQ+2:iR)- x(iQ+1:iR-1)).* ...
                 (dz(iQ+2:iR)+dz(iQ+1:iR-1)));
vero=vero+.5*(R-x(iR))*dz(iR);