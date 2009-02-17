fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fixy,[nx+1,ny+1],'double');
y=fread(fixy,[nx+1,ny+1],'double');
fclose(fixy)
figure(2);
fid=fopen('Hrms.dat','r');
fiz=fopen('zs.dat','r');
fiu=fopen('u.dat','r');
fiv=fopen('v.dat','r');
for i=1:nt;
    i
    f=fread(fid,[nx+1,ny+1],'double');
    z=fread(fiz,[nx+1,ny+1],'double');
    u=fread(fiu,[nx+1,ny+1],'double');
    v=fread(fiv,[nx+1,ny+1],'double');
    if i>0&mod(i,1)==0
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

fclose(fid)
fclose(fiz)
fclose(fiu)
fclose(fiv)

