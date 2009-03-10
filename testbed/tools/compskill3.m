function [r2,slope,eps,bss]=compskill(comp,meas,runid,testid,varname);
c=interp1(comp(:,1),comp(:,2),meas(:,1));
c(isnan(c))=0;
m=meas(:,2);
r2=corrcoef(c,m);
r2=r2(1,2);
slope=mean(c)/mean(m);
eps=std(c-m)/abs(mean(m));
bss=1-(std(c-m))^2/(std(m))^2
fi=fopen([varname '.err'],'wt');
fprintf(fi,'%s %s %6.3f %6.3f %6.3f %6.3f\n',runid,testid,r2,slope,eps,bss)
fclose(fi)
