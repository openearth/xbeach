function plotdata()
clear variables
close all
[runid,testid,datadir]=testinfo


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
fiz=fopen('zs.dat','r');
fiu=fopen('u.dat','r');
fib=fopen('zb.dat','r');
fih=fopen('H.dat','r');

zs=zeros(nt,nx+1);
Hrms=zeros(nt,nx+1);
uu=zeros(nt,nx+1);
dep=zeros(nt,nx+1);
htot=zeros(nt,nx+1);
g=9.81;

for i=1:nt;
    i/nt;
    z=fread(fiz,[nx+1,ny+1],'double');
    b=fread(fib,[nx+1,ny+1],'double');
    h=fread(fih,[nx+1,ny+1],'double');
    u=fread(fiu,[nx+1,ny+1],'double');    
    

    zs(i,:) = z(:,2)';
    Hrms(i,:) = h(:,2)';
    uu(i,:) = u(:,2)';
    dep(i,:) = b(:,2)';
    htot(i,:) = zs(i,:)-dep(i,:); 
    
%     plot(zs(i,:))
%     hold on
%     plot(uu(i,:),'r')
%     grid on
%     pause

%     if i>0; %mod(i,1)==0
%     figure(1)
%     clf
%     hmin=0.01;
%         ll = surf(x,y,b); colormap(copper); shading interp;
%     caxis([-20 2]);
%     hold on; 
%     z(z<=b+0.01)=NaN;
%     hh = surf(x,y,z); 
%     set(hh,'FaceColor',[0 0 1],'LineStyle','None'); 
%     alpha(hh,0.5);
%     ii = surf(x,y,h+z); 
%     set(ii,'FaceColor',[0 1 1],'LineStyle','None'); 
%     alpha(ii,0.5);
%     vv=axis;
%     axis([vv(1:4) -5 3]);
%     
%     %view(2); 
%     title(i)
%     cmapbrighten(-0.5);
%    
%     LIGHTING GOURAUD
    
%    end
end;
% 


meaneta=load([datadir 'meaneta.asc']);
plot(meaneta(:,1),meaneta(:,2),'g*');
hold on
meanzs = mean(zs(500:end,:));
plot(x,meanzs,'g')
title('blue = hm0hi, red=hm0lowtot, green = steadysetdown');

df = 1/1800;
[dum, zs] = filterhjb(zs,0.15,df,1); % low pass
[zs, dum] = filterhjb(zs,0.03,df,1); % hipass
[dum, uu] = filterhjb(uu,0.15,df,1); % low pass
[uu, dum] = filterhjb(uu,0.03,df,1); %hi pass
[dum, htot] = filterhjb(htot,0.15,df,1); % low pass



Hm0mean = sqrt(2)*sqrt(mean(Hrms(500:end,:).^2,1));
plot(x,Hm0mean,'b')
hold on
Hm0hi=load([datadir 'Hm0_high.asc']);
plot(Hm0hi(:,1),Hm0hi(:,2),'b*');

Hm0lo=load([datadir 'Hm0_low.asc']);
plot(Hm0lo(:,1),Hm0lo(:,2),'r*');
Hm0low = 4*(std(zs(500:end,:)));
plot(x,Hm0low,'r')


fp=0.3;
k=disper(2*pi*fp,htot);
cg= 0.5*(1+2*k.*htot./sinh(2*k.*htot))*2*pi*fp./k;

g=9.81;
zetai=(sqrt(g*htot).*(detrend(zs))+detrend(uu).*htot)./...
   (sqrt(g*htot)+cg);
zetar=(cg.*(detrend(zs))-detrend(uu).*htot)./ ...
   (sqrt(g*htot)+cg);
Hm0lowin = 4*(std(detrend(zetai(500:end,:))));
Hm0lowout = 4*(std(detrend(zetar(500:end,:))));
plot(x,Hm0lowin,'c--')
plot(x,Hm0lowout,'k--')
Hrmsloin=load([datadir 'Hrms_inc_LF.asc']);
plot(Hrmsloin(:,1),Hrmsloin(:,2)*sqrt(2),'c*');
Hrmsloout=load([datadir 'Hrms_out_LF.asc']);
plot(Hrmsloout(:,1),Hrmsloout(:,2)*sqrt(2),'k*');
axis([0 30 -.05 .15])
grid on
fclose(fiz)
fclose(fib)
pname = ['..\..\report\',testid '_' runid '_fig1' '.jpg'];
eval(['print -djpeg ' pname]);
