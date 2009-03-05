clear all;close all
load bot2e.mat;
x=[0:198];
figure(1);
ii=0
for i=1:19;
    if sum(zp(:,i)~=0)
        ii=ii+1;
        z=interp1(xp,zp(:,i),x);
        z(isnan(z))=0;
        a=num2str(100+i-1);
        a=a(2:3);
        fname=['2E' a '.dep']
        fi=fopen(fname,'wt');
        fprintf(fi,'%6.2f ',z);fprintf(fi,'\n');
        fprintf(fi,'%6.2f ',z);fprintf(fi,'\n');
        fprintf(fi,'%6.2f ',z);fprintf(fi,'\n');
        fclose(fi)
        if i==1
            z0=min(z,5.9);
        end
        if i==9
            tek(:,1)=x;
            tek(:,2)=z-z0;
            tekal('write','sedero.tek',tek)
        end
        level=4.6;
        [xi,zi]=intersections(x,z0,x,z);
        [vsed(ii),vero(ii),R(ii)]=profana(x,z0,z,level);
        time(ii)=(i-1);
        plot(x,z0,x,z,xi,zi,'o');
        title(fname);
    end
end