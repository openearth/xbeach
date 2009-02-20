function [r2,slope,eps]=compskill(comp,meas);
c=interp1(comp(:,1),comp(:,2),meas(:,1));
m=meas(:,2);
r2=corrcoef(c,m);
r2=r2(1,2);
slope=mean(c)/mean(m);
eps=std(c-m)/mean(m);
