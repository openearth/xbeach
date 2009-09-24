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
fi5=fopen('ue.dat','r');
i0=0;
count=0
for i=1:nt;
    f=fread(fid,[nx+1,ny+1],'double');
    ft(i)=f(2,1);
    z=fread(fi2,[nx+1,ny+1],'double');
    hrms=fread(fi3,[nx+1,ny+1],'double');
    u=fread(fi4,[nx+1,ny+1],'double');
    ue=fread(fi5,[nx+1,ny+1],'double');
    if i==1
        z0=z;
    end
    if i>2&ft(i-1)<ft(i-2)&ft(i)>ft(i-1)&i0==0
        i0=i
    end

    if mod(i-i0,20)==0&i0~=0&count<20
        count=count+1
        subplot(211);
        f(u==0)=nan;
        u(f-z<.01)=nan;
        plot(x,f,x,z,'b-','linewidth',2);
        set(gca,'ylim',[min(min(f)),max(max(f))]);
        title(num2str(i));hold on;
        %plot(x,z,'r-','linewidth',2);
        %plot(x,z0,'k-','linewidth',2);
        %plot(x,f+hrms/2,'g-','linewidth',2);
        %plot(x,f-hrms/2,'g-','linewidth',2);
        hold on
%        axis([150 250 3 6]);
        % axis([0 250 3 6]);
        drawnow;
        subplot(212)
        plot(x,u,'b-','linewidth',2);title(num2str(i));hold on;
        %plot(x,ue,'r-','linewidth',2);
        hold on
%        axis([150 250 -2 2]);
         %axis([0 250 -2 2]);
        drawnow;
    end
end;
fclose(fid)