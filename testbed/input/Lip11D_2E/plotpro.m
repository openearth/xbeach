function plotpro()
[runid,testid,datadir]=testinfo
fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fid,[nx+1,ny+1],'double');
y=fread(fid,[nx+1,ny+1],'double');
fclose(fixy)
fid=fopen('zs.dat','r');
fi2=fopen('zb.dat','r');
fi3=fopen('H.dat','r');
clear z
for i=1:nt;
    zb=fread(fi2,[nx+1,ny+1],'double');
    z{i}=zb(:,2);
    if i==1
        z0=z;
    end
end
load ([datadir 'bot2e'])
x=x+1;
figure;
subplot(221)
plot(xp,zp(:,1),'k-',x(:,2),z{661},'b--',xp,zp(:,2),'b-','linewidth',2)
legend('0 hr','1 hr comp.','1 hr meas.','location','northwest')
xlabel('x (m)');ylabel('z_b (m)');axis([100 200 2 6])
subplot(222)
plot(xp,zp(:,1),'k-',x(:,2),z{1021},'b--',xp,zp(:,3),'b-','linewidth',2)
legend('0 hr','2 hr comp.','2 hr meas.','location','northwest')
xlabel('x (m)');ylabel('z_b (m)');axis([100 200 2 6])
subplot(223)
plot(xp,zp(:,1),'k-',x(:,2),z{1741},'b--',xp,zp(:,5),'b-','linewidth',2)
legend('0 hr','4 hr comp.','4 hr meas.','location','northwest')
xlabel('x (m)');ylabel('z_b (m)');axis([100 200 2 6])
subplot(224)
plot(xp,zp(:,1),'k-',x(:,2),z{nt},'b--',xp,zp(:,9),'b-','linewidth',2)
legend('0 hr','8 hr comp.','8 hr meas.','location','northwest')
xlabel('x (m)');ylabel('z_b (m)');axis([100 200 2 6])
pname = ['..\..\report\' testid '_' runid '_fig' '2' '.jpg'];
eval(['print -djpeg ' pname]);
fclose(fid)
fclose(fi2)
