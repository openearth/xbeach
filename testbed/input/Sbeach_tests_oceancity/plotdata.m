function plotdata()
figure(3);
fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fid,[nx+1,ny+1],'double')-300;
y=fread(fid,[nx+1,ny+1],'double');
fclose(fixy)
fid=fopen('zs.dat','r');
fi2=fopen('zb.dat','r');
fi3=fopen('H.dat','r');
%fi4=fopen('u.dat','r');
fi5=fopen('ue.dat','r');

for i=1:nt;
    f=fread(fid,[nx+1,ny+1],'double');
    z=fread(fi2,[nx+1,ny+1],'double');
    hrms=fread(fi3,[nx+1,ny+1],'double');
    %u=fread(fi4,[nx+1,ny+1],'double');
    ue=fread(fi5,[nx+1,ny+1],'double');
    if i==1
        z0=z;
    end
    zs0t(i)=f(1,2);
    if mod(i,10)==0&i>500
        subplot(211);
        plot(x,f,'b-','linewidth',2);title(num2str(i));hold on;
        plot(x,z,'r-','linewidth',2);
        plot(x,z0,'k-','linewidth',2);
        plot(x,f+hrms/2,'g-','linewidth',2);
        plot(x,f-hrms/2,'g-','linewidth',2);
        hold off
        axis([-300 1000 -15 10]);
%         axis([0 250 3 6]);
        drawnow;%pause
        subplot(212)
        %plot(x,u,'b-','linewidth',2);title(num2str(i));hold on;
        plot(x,ue,'r-','linewidth',2);
        hold off
         axis([-300 1000 -1 1]);
% %         axis([0 250 -2 2]);
         drawnow;
    end
end;
fclose(fid)
t0=[1:length(zs0t)]/120;