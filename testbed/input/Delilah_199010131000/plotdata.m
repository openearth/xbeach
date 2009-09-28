function plotdata()
clear all;close all
[runid,testid,datadir]=testinfo
%sdir=pwd;
%cd w:\postbox\mccall\duck1990\19901013\
load ([datadir 'measuredata.mat']);
% load spectrumdata5.mat; 
load ([datadir 'Pdatanew.mat']);

%cd i:\santarosa\the'sis tests'\test11_duck_run18_13_10_newexe\
sfigs=1;
[nt,x,y]=getdimensions;
Hmean=0*x;
fid=fopen('H.dat','r');

start=600;
stop=3600;
% start=850;
% stop=1750;
le=stop-start;

% nsx= [166 162 157 153 150 146 142 132 117];  % regular
nsx= [112 108 103  99  96  92  88  78  63];  % irregular
% irreg grid 44
% small grid reg 34
% large grid reg 54
measurerow=44;
dx=5;
progressbar(0);
for i=1:stop
    hrms=fread(fid,size(x),'double');
    progressbar(i/stop);
    if i>start
        Hmean=Hmean+hrms.^2;
    end
end

Hmean=sqrt(Hmean/(stop-start));

frewind(fid);
VarH=0*x;
progressbar(0);
for i=1:stop
    hrms=fread(fid,size(x),'double');
    progressbar(i/stop);
    if i>start
        VarH=VarH+(Hmean-hrms).^2;
    end
end

VarH=VarH/(stop-start);
VarH=VarH.^0.5;

fclose(fid);
figure(1);
hold on;
% for i=1:71
% plot(x(:,1),Hmean(:,i),'b.');
% end
for i=1:9
eval(['plot(repmat(x(nsx(i),1),1,length(measuredata.H.times)),measuredata.H.Hs_measurement_inst',num2str(i),'/sqrt(2),''ro'')']);
end
plot(x(:,1),Hmean(:,measurerow),'color','k','linewidth',1.5);
plot(x(:,1),Hmean(:,measurerow)+VarH(:,measurerow),'color','k','linewidth',1);
plot(x(:,1),Hmean(:,measurerow)-VarH(:,measurerow),'color','k','linewidth',1);
title('RMS-wave height');
xlabel('cross shore direction (m)');
ylabel('H_{rms} (m)');
hold off;
if sfigs
   PN = fliplr(pwd);
   [runid,R] = strtok(PN,'\'); runid = fliplr(runid);
   [testid,R] = strtok(R,'\'); testid = fliplr(testid);
   [dataid,R] = strtok(R,'\'); dataid = fliplr(dataid);
   pname = [dataid '_' testid '_' runid '_hrmspoints.jpg'];
   eval(['print -djpeg ' pname]);
end

Umean=0*x;
Vmean=0*x;
fidu=fopen('ue.dat','r');
fidv=fopen('ve.dat','r');

progressbar(0);
for i=1:stop
    u=fread(fidu,size(x),'double');
    v=fread(fidv,size(x),'double');
    progressbar(i/stop);
    if i>start
        Umean=Umean+u;
        Vmean=Vmean+v;        
    end
end

Umean=Umean/(stop-start);
Vmean=Vmean/(stop-start);
Smean=pyth(Umean,Vmean);

frewind(fidu);frewind(fidv);
VarS=0*x;
VarU=0*x;
VarV=0*x;
progressbar(0);
for i=1:stop
    u=fread(fidu,size(x),'double');
    v=fread(fidv,size(x),'double');
    s=pyth(u,v);
    progressbar(i/stop);
    if i>start
        VarU=VarU+(Umean-u).^2;
        VarV=VarV+(Vmean-v).^2;
        VarS=VarS+(Smean-s).^2;
    end
end

VarU=VarU/(stop-start);
VarU=VarU.^0.5;
VarV=VarV/(stop-start);
VarV=VarV.^0.5;
VarS=VarS/(stop-start);
VarS=VarS.^0.5;

fclose(fidu);fclose(fidv);

figure(2);
hold on;
% for i=1:71
% plot(x(:,1),Umean(:,i),'b.');
% end
for i=1:9
eval(['plot(repmat(x(nsx(i),1),1,length(measuredata.U.times)),measuredata.U.u_measurement_inst',num2str(i),',''ro'')']);
end
plot(x(:,1),Umean(:,measurerow),'color','k','linewidth',1.5);
plot(x(:,1),Umean(:,measurerow)+VarU(:,measurerow),'color','k','linewidth',1);
plot(x(:,1),Umean(:,measurerow)-VarU(:,measurerow),'color','k','linewidth',1);
title('Cross shore velocity');
xlabel('cross shore direction (m)');
ylabel('velocity (m/s)');
ylim([-0.4 1.4]);
hold off;
if sfigs
   pname = [dataid '_' testid '_' runid '_upoints.jpg'];
   eval(['print -djpeg ' pname]);
end

figure(3);
hold on;
% for i=1:71
% plot(x(:,1),Vmean(:,i),'b.');
% end
for i=1:9
eval(['plot(repmat(x(nsx(i),1),1,length(measuredata.V.times)),-measuredata.V.v_measurement_inst',num2str(i),',''ro'')']);
end
plot(x(:,1),Vmean(:,measurerow),'color','k','linewidth',1.5);
plot(x(:,1),Vmean(:,measurerow)+VarV(:,measurerow),'color','k','linewidth',1);
plot(x(:,1),Vmean(:,measurerow)-VarV(:,measurerow),'color','k','linewidth',1);
title('Longshore velocity');
xlabel('cross shore direction (m)');
ylabel('velocity (m/s)');
ylim([-0.5 2.5]);
hold off;
if sfigs
   pname = [dataid '_' testid '_' runid '_vpoints.jpg'];
   eval(['print -djpeg ' pname]);
end

figure(4);
hold on;
% for i=1:71
% plot(x(:,1),Smean(:,i),'b.');
% end
for i=1:9
eval(['plot(repmat(x(nsx(i),1),1,length(measuredata.V.times)),pyth(measuredata.U.u_measurement_inst',num2str(i),',measuredata.V.v_measurement_inst',num2str(i),'),''ro'')']);
end
plot(x(:,1),Smean(:,measurerow),'color','k','linewidth',1.5);
plot(x(:,1),Smean(:,measurerow)+VarS(:,measurerow),'color','k','linewidth',1);
plot(x(:,1),Smean(:,measurerow)-VarS(:,measurerow),'color','k','linewidth',1);
title('Total velocity');
xlabel('cross shore direction (m)');
ylabel('velocity (m/s)');
ylim([0 2.5]);
hold off;
if sfigs
   pname = [dataid '_' testid '_' runid '_tvpoints.jpg'];
   eval(['print -djpeg ' pname]);
end

figure(5);
pcolor(x,y,pyth(Umean,Vmean));shading flat;
% pcolor(x,y,Umean);shading flat;
%plotcontour(1);
hold on
quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),Umean(1:2:end,1:2:end),Vmean(1:2:end,1:2:end),'w');
plot(x(nsx,measurerow),y(nsx,measurerow),'mo');
hold off
title('Mean current');
%colorbarrob('[m/s]');
colorbar
xlabel('cross shore position (m)');
ylabel('longshore position (m)');
if sfigs
   pname = [dataid '_' testid '_' runid '_velspatial.jpg'];
   eval(['print -djpeg ' pname]);
end



%%

fidz=fopen('zs.dat','r');
progressbar(0);
ZS=zeros(le,9);
mm=0;
for i=1:stop
    zs=fread(fidz,size(x),'double');
    progressbar(i/stop);
    if i>start
        mm=mm+1;
        ZS(mm,:)=zs(nsx,measurerow);
    end
end
ZS=detrend(ZS);
Pnew=detrend(Pnew);
T = stop-start;
df = 1/T;

for i=1:9
    temp=fft(ZS(:,i),(le),1)/(le);
    temp2=fft(Pnew(1:le,i),(le),1)/(le);
    Y(:,i)=temp;
    Y2(:,i)=temp2;
end
n=length(Y);
for i=1:9
    powerxb(:,i)=2*(Y(1:floor(n/2),i)).*conj(Y(1:floor(n/2),i))/df;
    power(:,i)=  2*(Y2(1:floor(n/2),i)).*conj(Y2(1:floor(n/2),i))/df;
end
nyquist = 0.5;
freqxb = (1:n/2)/(n/2)*nyquist;
dfxb=freqxb(2)-freqxb(1);
freduced=[0.005:0.005:0.5];

for ii=1:9
    for i=1:length(freduced)-1
        temp=find(freqxb<=freduced(i+1) & freqxb>freduced(i));
        Sxbreduced(i,ii)=sum(powerxb(temp,ii))*dfxb/0.005;
        Sreduced(i,ii)=sum(power(temp,ii))*dfxb/0.005;
    end
end



% dt = 1;
% T = stop-start;
% df = 1/T;
% fnyq = 1/2/dt;
% f = 0:df:fnyq;
% Fp = fft(ZS,m,1)/m;
% for icounter=1:9
%     Sxb(:,icounter) = 2*Fp(:,icounter).*conj(Fp(:,icounter))/df;
% end

% for ii=1:9
%     for i=1:length(freduced)-1
%         temp=find(f<=freduced(i+1) & f>freduced(i));
%         Sxbreduced(i,ii)=sum(Sxb(temp,ii));
%     end
% end

% Sxbreduced=Sxbreduced*0.005/df;

for i=1:9
    figure
%     bar(freduced(2:end-1),[Sreduced(2:end,i) Sxbreduced(2:end,i)]);
    plot(freqxb,power(:,i)/(freqxb(2)-freqxb(1)),'b');
% %     hold on;plot(freqxb(2:end),powerxb(:,i)/(freqxb(2)-freqxb(1)),'r');hold off;
    hold on;plot(freqxb,powerxb(:,i)/(freqxb(2)-freqxb(1)),'r');
    xlim([0 0.1]);
%     ylim([0 0.4]);
    legend('measure','XBeach','location','Northeast');
    title(['Low frequency spectrum at location ',num2str(i*10)]);
    xlabel('frequency (Hz)');
    ylabel('variance density');
    if sfigs
      pname = [dataid '_' testid '_' runid,'spectrum',num2str(i),'.jpg'];
      eval(['print -djpeg ' pname]);
    end
    
end

figure
title('H_{rms} infragravity');
hold on;
ifrange=find(freqxb>=0.005 & freqxb<=0.05);
for i=1:9
    hrmslowm(i)=sqrt(8*sum(power(ifrange,i))*df);
    hrmslowxb(i)=sqrt(8*sum(powerxb(ifrange,i))*df);
    plot(x(nsx(i),1),hrmslowm(i),'b*',x(nsx(i),1),hrmslowxb(i),'ro');
end
hold off
title('H_{rms} infragravity');
xlabel('cross shore direction (m)');
ylabel('H_{rms,ifg} (m)');
legend('measured','modelled','location','SouthWest');
if sfigs
   pname = [dataid '_' testid '_' runid '_hrmsinf.jpg'];
   eval(['print -djpeg ' pname]);
end
pname = ['..\..\report\',testid '_' runid '_Hrmslo' '.jpg'];
eval(['print -djpeg ' pname]);

