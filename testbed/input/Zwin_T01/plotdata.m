[runid,testid,datadir]=testinfo
xmin=-20;xmax=50;ymin=-20;ymax=20;
fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fid,[nx+1,ny+1],'double');
y=fread(fid,[nx+1,ny+1],'double');
fclose(fixy)
fid=fopen('zb.dat','r');
fiz=fopen('zs.dat','r');
fiu=fopen('u.dat','r');
fiv=fopen('v.dat','r');
fic=fopen('cc.dat','r');
fifx=fopen('Fx.dat','r');
fify=fopen('Fy.dat','r');
t=[1:nt]*0.5-0.5;
first=1
for i=1:nt;
    i
    zb=fread(fid,[nx+1,ny+1],'double');
    z=fread(fiz,[nx+1,ny+1],'double');
    u=fread(fiu,[nx+1,ny+1],'double');
    v=fread(fiv,[nx+1,ny+1],'double');
    cc=fread(fic,[nx+1,ny+1],'double');
    fy=fread(fify,[nx+1,ny+1],'double');
    zst(:,i)=z(:,ny/2+1);
    zbt(:,i)=zb(:,ny/2+1);
    ut(:,i)=u(:,ny/2+1);
    vt(:,i)=v(:,ny/2+1);
    cct(:,i)=cc(:,ny/2+1);
    Bt(i)=min(y(zb==3.3&abs(x)<20&y>0))-max(y(zb==3.3&abs(x)<20&y<0));
end;

fclose(fid)
fclose(fiz)
fclose(fiu)
fclose(fiv)
a=load([datadir 'fig1\ms1.xyz']);
t1=a(:,1)/1000;w1=a(:,2)/10000;
for i=1:12;w1(2:end-1)=.25*w1(1:end-2)+.5*w1(2:end-1)+.25*w1(3:end);end
a=load([datadir 'fig1\ms2.xyz']);
t2=a(:,1)/1000;w2=a(:,2)/10000;
a=load([datadir 'fig1\ms3.xyz']);
t3=a(:,1)/1000;w3=a(:,2)/10000;
a=load([datadir 'fig1\ms4.xyz']);
t4=a(:,1)/1000;w4=a(:,2)/10000;
a=load([datadir 'fig1\ms5.xyz']);
t5=a(:,1)/1000;w5=a(:,2)/10000;
figure(4);
subplot(311)
plot(t2,w2,'k.',t,zst(35,:),'k-',t3,w3,'b.',t,zst(39,:),'b-',t4,w4,'r.',t,zst(103,:),'r-','linewidth',2)
legend('MS2 obs','MS2 comp','MS3 obs','MS3 comp','MS4 obs','MS4 comp')
a=load([datadir 'fig2\ms2.xyz']);
t2=a(:,1)/1000;v2=a(:,2)/10000;
a=load([datadir 'fig2\ms3.xyz']);
t3=a(:,1)/1000;v3=a(:,2)/10000;
a=load([datadir 'fig2\ms4.xyz']);
t4=a(:,1)/1000;v4=a(:,2)/10000;
a=load([datadir 'fig2\ms5.xyz']);
t5=a(:,1)/1000;v5=a(:,2)/10000;
subplot(312)
plot(t2,v2,'k.',t,ut(35,:),'k-',t3,v3,'b.',t,ut(39,:),'b-',t4,v4,'r.',t,ut(103,:),'r-','linewidth',2)
legend('MS2 obs','MS2 comp','MS3 obs','MS3 comp','MS4 obs','MS4 comp')
Bmeas=load([datadir 'B.txt'])
subplot(313)
plot(Bmeas(:,1),Bmeas(:,2),'.k',t,Bt,'k-','linewidth',2)
legend('B obs','B comp')