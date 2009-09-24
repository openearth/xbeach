figure(1);
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
fi4=fopen('u.dat','r');
fi6=fopen('Qb.dat','r');

for i=1:nt;
    f=fread(fid,[nx+1,ny+1],'double');
    z=fread(fi2,[nx+1,ny+1],'double');
    hrms=fread(fi3,[nx+1,ny+1],'double');
    u=fread(fi4,[nx+1,ny+1],'double');
    Qb=fread(fi6,[nx+1,ny+1],'double');
    if i==1
        z0=z;
    end

    if mod(i,1)==0&i>300
        subplot(211);
        plot(f,'b-','linewidth',2);title(num2str(i));hold on;
        plot(z,'r-','linewidth',2);
        plot(z0,'k-','linewidth',2);
        plot(f+hrms/2,'g-','linewidth',2);
        plot(f-hrms/2,'g-','linewidth',2);
        hold off
%        axis([150 200 0 6]);
%        axis([0 250 3 6]);
        drawnow;
        subplot(212)
        plot(u,'b-','linewidth',2);title(num2str(i));hold on;
        plot(Qb,'k-','linewidth',2);
        hold off
 %       axis([150 200 -2 2]);
 %        axis([0 250 -2 2]);
        drawnow;
    end
end;
fclose(fid)