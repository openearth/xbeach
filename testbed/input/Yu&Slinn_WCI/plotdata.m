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
fid=fopen('H.dat','r');
fiz=fopen('zb.dat','r');
fiu=fopen('u.dat','r');
fiv=fopen('v.dat','r');
for i=1:nt;
    i
    f=fread(fid,[nx+1,ny+1],'double');
    z=fread(fiz,[nx+1,ny+1],'double');
    u=fread(fiu,[nx+1,ny+1],'double');
    v=fread(fiv,[nx+1,ny+1],'double');
    
    if i==nt
        pcolor(x,y,f);shading interp;caxis([0 1.6]);colorbar;axis equal;
        axis([0 700 0 500]);set(gca,'fontsize',14)
        hold on;quiver(x,y,50*u,50*v,0,'k');
        [cs,h]=contour(x,y,z,[-20:2:-6 -4:1:0],'w','linewidth',2);
        clabel(cs,'fontsize',12)
        drawnow;hold off
        title('with wave current interaction');
        xlabel('X (m)','fontsize',14)
        ylabel('Y (m)','fontsize',14)
	PN = fliplr(pwd);
	[runid,R] = strtok(PN,'\'); runid = fliplr(runid);
	[testid,R] = strtok(R,'\'); testid = fliplr(testid);
	[dataid,R] = strtok(R,'\'); dataid = fliplr(dataid);
	pname = [dataid '_' testid '_' runid '_fig' '1' '.jpg'];
	eval(['print -djpeg ' pname]);
    end

end;

fclose(fid)
fclose(fiz)
fclose(fiu)
fclose(fiv)

