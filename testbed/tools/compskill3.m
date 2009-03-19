function [r2,slope,eps,bss]=compskill3(comp,meas,runid,testid,varname);
c=interp1(comp(:,1),comp(:,2),meas(:,1));
m=meas(:,2)
cred=c(~isnan(c)&abs(m)>.05*max(abs(m)))
mred=m(~isnan(c)&abs(m)>.05*max(abs(m)))
%r2=corrcoef(cred,mred);
%r2=r2(1,2);
r2=mean((cred-mean(cred)).*(mred-mean(mred)))/(std(mred)*std(cred));
if abs(mean(mred))<std(mred)
    slope=mean((cred-mean(cred)).*(mred-mean(mred)))/mean((mred-mean(mred)).^2);
    eps=std(cred-mred)/std(mred);
else
    slope=mean(cred)/mean(mred);
    eps=std(cred-mred)/abs(mean(mred));
end
bss=1-(std(cred-mred))^2/var(mred);
fi=fopen([varname '.err'],'wt');
fprintf(fi,'%s %s %6.3f %6.3f %6.3f %6.3f\n',runid,testid,r2,slope,eps,bss)
fclose(fi)
