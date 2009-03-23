function plotdata()
clear all
ix=65;iy=40
fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fid,[nx+1,ny+1],'double');
y=fread(fid,[nx+1,ny+1],'double');
fclose(fixy)
figure(1);
fib=fopen('zb.dat','r');
fid=fopen('H.dat','r');
fiz=fopen('zs.dat','r');
fiu=fopen('u.dat','r');
fiv=fopen('v.dat','r');
fify=fopen('Fy.dat','r');
ii=0;
for i=1:nt;
    zb=fread(fib,[nx+1,ny+1],'double');
    f=fread(fid,[nx+1,ny+1],'double');
    z=fread(fiz,[nx+1,ny+1],'double');
    u=fread(fiu,[nx+1,ny+1],'double');
    v=fread(fiv,[nx+1,ny+1],'double');
    fy=fread(fify,[nx+1,ny+1],'double');
    if i>387
    ii=ii+1;
    Hrmst(ii,1:ny+1)=f(ix,1:ny+1);
    zst(ii,1:ny+1)=z(ix,1:ny+1);
    ut(ii,1:ny+1)=u(ix,1:ny+1);
    vt(ii,1:ny)=v(ix,1:ny);
        title(num2str(i));
    subplot(141);pcolor(x,y,f);shading interp;caxis([0 2]);colorbar;axis equal;title('H')
    set(gca,'xticklabel',[0 500]);set(gca,'xtick',[0 500]);
    subplot(142);pcolor(x,y,z);shading interp;caxis([-.5 .5]);colorbar;axis equal;title('zs')
    set(gca,'xticklabel',[0 500]);set(gca,'xtick',[0 500]);set(gca,'yticklabel',[]);    
    subplot(143);pcolor(x,y,u);shading interp;caxis([-1 1]);colorbar;axis equal;title('u')
    set(gca,'xticklabel',[0 500]);set(gca,'xtick',[0 500]);set(gca,'yticklabel',[]);    
    subplot(144);pcolor(x,y,v);shading interp;caxis([-2 2]);colorbar;axis equal;title('v');
    set(gca,'xticklabel',[0 500]);set(gca,'xtick',[0 500]);set(gca,'yticklabel',[]);drawnow;    
    fname=strcat(num2str(1000+i),'.jpg')
    print('-djpeg',fname)
    end
end;
fclose(fib)
fclose(fid)
fclose(fiz)
fclose(fiu)
fclose(fiv)
fclose(fify)

figure(2);
subplot(421);plot(Hrmst');
subplot(423);plot(zst');
subplot(425);plot(ut');
subplot(427);plot(vt');
subplot(422);plot(Hrmst(:,iy));
subplot(424);plot(zst(:,iy));
subplot(426);plot(ut(:,iy));
subplot(428);plot(vt(:,iy));
