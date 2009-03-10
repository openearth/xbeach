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
fih=fopen('hh.dat','r');
fiu=fopen('u.dat','r');
fiv=fopen('v.dat','r');
fic=fopen('ccg.dat','r');
fisx=fopen('Sug.dat','r');
fisy=fopen('Svg.dat','r');
t=[1:nt]*0.5-0.5;
first=1
for i=1:nt;
    i
    zb=fread(fid,[nx+1,ny+1],'double');
    zs=fread(fiz,[nx+1,ny+1],'double');
    hh=fread(fih,[nx+1,ny+1],'double');
    u=fread(fiu,[nx+1,ny+1],'double');
    v=fread(fiv,[nx+1,ny+1],'double');
    cc=fread(fic,[nx+1,ny+1],'double');
    sx=fread(fisx,[nx+1,ny+1],'double');
    sy=fread(fisy,[nx+1,ny+1],'double');
    zst(:,i)=z(:,ny/2+1);
    zbt(:,i)=zb(:,ny/2+1);
    ut(:,i)=u(:,ny/2+1);
    vt(:,i)=v(:,ny/2+1);
    cct(:,i)=cc(:,ny/2+1);
    subplot(231);pcolor(x,y,zb);shading flat;colorbar;axis([xmin xmax ymin ymax]);title('zb')
    subplot(232);pcolor(x,y,zs);shading flat;colorbar;axis([xmin xmax ymin ymax]);title('zs')
    subplot(233);pcolor(x,y,hh);shading flat;colorbar;axis([xmin xmax ymin ymax]);title('hh')
    subplot(234);pcolor(x,y,u);shading flat;colorbar;axis([xmin xmax ymin ymax]);title('u')
    subplot(235);pcolor(x,y,v);shading flat;colorbar;axis([xmin xmax ymin ymax]);title('v')
    subplot(236);pcolor(x,y,cc);shading flat;colorbar;axis([xmin xmax ymin ymax]);title('cc')
    drawnow;pause
end