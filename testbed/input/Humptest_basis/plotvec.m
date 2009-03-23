fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fid,[nx+1,ny+1],'double');
y=fread(fid,[nx+1,ny+1],'double');
fclose(fixy)
ntint=2000
figure(1);
fid=fopen('Hrms.dat','r');
fiz=fopen('zb.dat','r');
fiu=fopen('u.dat','r');
fiv=fopen('v.dat','r');
for i=1:nt
    f=fread(fid,[nx+1,ny+1],'double');
    z=fread(fiz,[nx+1,ny+1],'double');
    if i==1
        z0=z;
    end
    u=fread(fiu,[nx+1,ny+1],'double');
    v=fread(fiv,[nx+1,ny+1],'double');
    if i>0;%mod(i,1)==0
        title(num2str(i));
    subplot(121);pcolor(x,y,f);shading interp;caxis([0 4]);colorbar;axis([0 700 0 800]);axis equal;
    subplot(122);pcolor(x,y,z-z0);shading interp;caxis([-1 1]);colorbar;axis([0 700 0 800]);axis equal;
    hold on 
    quiver(x,y,u,v,1); hold off;title(num2str(i));drawnow
    end
end;
fclose(fid)

