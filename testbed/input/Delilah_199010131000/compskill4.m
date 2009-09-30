function [r2,sci,relbias,bss]=compskill4(comp,meas,runid,testid,varname);
c=interp1(comp(:,1),comp(:,2),meas(:,1));
m=meas(:,2)
cred=c(~isnan(c)&abs(m)>.05*max(abs(m)))
mred=m(~isnan(c)&abs(m)>.05*max(abs(m)))
r2=mean((cred-mean(cred)).*(mred-mean(mred)))/(std(mred)*std(cred));
rmserr=sqrt(mean((cred-mred).^2));
rmsm=sqrt(mean(mred.^2));
sci=rmserr/max(rmsm,abs(mean(mred)))
relbias=mean(cred-mred)/max(rmsm,abs(mean(mred)))
bss=1-(std(cred-mred))^2/var(mred);
fi=fopen([varname '.err'],'wt');
fprintf(fi,'%s %s %s %s %s %s\n','runid','testid','r2','sci','relbias','bss')
fprintf(fi,'%s %s %6.3f %6.3f %6.3f %6.3f\n',runid,testid,r2,sci,relbias,bss)
fclose(fi)
