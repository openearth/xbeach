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
    if i>0&mod(i,1)==0
        %          figure(2);
        %             subplot(221);pcolor(x,y,zb);shading interp;axis ([xmin xmax ymin ymax]);title('H')
        %             caxis([-4 4]);colorbar;
        %             subplot(222);pcolor(x,y,z);shading interp;axis ([xmin xmax ymin ymax]);title('zs')
        %             caxis([0.5 3]);colorbar;
        %             subplot(223);pcolor(x,y,u);shading interp;axis ([xmin xmax ymin ymax]);title('u')
        %             caxis([-5 5]);colorbar;
        %             subplot(224);pcolor(x,y,v);shading interp;axis ([xmin xmax ymin ymax]);title('v');
        %             caxis([-5 5]);colorbar;
        % %     figure(1);plot(u(:,end-4:end));legend('end-4','end-3','end-2','end-1','end')
        % %        title(num2str(i));
        % %     drawnow;pause
        % %     fname=strcat(num2str(1000+i),'.jpg')
        % %     print('-djpeg',fname)
        z(z-zb<0.001)=zb(z-zb<0.001)-.001;
        if first
            figure(1);set(gcf,'color','w')
            subplot('position',[0 .2 1 .8])
            b=surf(x(40:110,17:85),y(40:110,17:85),zb(40:110,17:85));
            material dull
            hold on;s=surf(x(40:110,17:85),y(40:110,17:85),z(40:110,17:85));shading interp;
            axis([-50 100 -50 50 -3 4]);daspect ([10 10 1])
            set(gca,'linewidth',2')
            set(b,'facecolor',[1 .8 0]);
            set(s,'facecolor','b','facealpha',0.3);
            material shiny
            daspect([10 10 1]);caxis([-5 5])
            camorbit(90,0);camlight;lighting phong;
            first=0;pause
        else
            set(b,'zdata',zb(40:110,17:85))
            set(s,'zdata',z(40:110,17:85))
            set(s,'cdata',u(40:110,17:85))
            camorbit(90/130,0)
            %camorbit(2,0)
        end
        hold off;
    end
    Bt(i)=min(y(zb==3.3&abs(x)<20&y>0))-max(y(zb==3.3&abs(x)<20&y<0));
    subplot('position',[0.15 0.1 .7 .1]);
    plot(t(1:i),Bt(1:i),'-',t(i),Bt(i),'o','linewidth',2);axis([0 70 0 50]);xlabel('t (min)');ylabel('Width (m)');
    set(gca,'linewidth',2')
    drawnow
    % fname=strcat(num2str(1000+i),'.jpg')
    % print('-djpeg',fname)
end;
% for i=1:360;
%     camorbit(1,0);
%     drawnow
% end
     fname=strcat('Zwin_t01.jpg')
     print('-djpeg',fname)

fclose(fid)
fclose(fiz)
fclose(fiu)
fclose(fiv)

