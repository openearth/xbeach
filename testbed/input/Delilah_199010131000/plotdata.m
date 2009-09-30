z = [ ];
day = 13;
jd  = 15;   % hour interval (6 corresponds to 7th hour, see line 10) 
fch = 1.0;
nc =  1;
g = 9.81;
rho = 1023;

cl={'ks' 'k' 'k' 'k--' 'ks' 'k' 'ks' 'k' 'k' 'kx'};
% cl={'bs' 'b' 'b' 'r--' 'rs' 'r' 'gs' 'g' 'k' 'kx'}; 


for jjd = 1:1:nc
   jd = jd + 1
if jd == 24
day = day + 1
jd = 0
end

XB=getdimensions;
zsxb=readvar('zs.dat',XB);
Hxb=readvar('H_mean.dat',XB);
vxb=readvar('ve_mean.dat',XB);
vxbm=mean(vxb,3);
Hxbm=sqrt(mean(Hxb.^2,3));
zsxbm=zsxb(:,41,:);zsxbm=squeeze(zsxbm);
Snn=0;
Snnf=0;
hrmslo=0;
ni=3;
anastart=0;
% anastart=0;
for i=1:ni
    st=(i-1)*(3600/ni)+1+anastart;
    en=i*(3600/ni)+anastart;
    [Snni Snnfi fxb hrmsloi]=makespectrum(zsxbm(:,st:en),1,1,1,10,0);
    Snn=Snn+1/ni*Snni;
    Snnf=Snnf+1/ni*Snnfi;
    hrmslo=hrmslo+1/ni*hrmsloi;
end
    
xpoints=fliplr([580 655 705 725 745 760 780 805 825]);
nsx=(round(xpoints/5)+1);


load del_hrms;
Hrms_m = interp2(xmeter,thm,Hrms,xmeter,day+(jd+0.5)/24);  % evaluated at the middle of the hour
f1=figure;
s1=subplot(411);
plot(950-xmeter,Hrms_m,cl{1},'markersize',5);
hold on
p1=plot(XB.x(:,41),Hxbm(:,41),cl{2});set(p1,'linewidth',1.5);
hold off
xlabel('x (m)')
ylabel('H_{rms,high} (m)')

% find corresponding low freq wave height
if jd+1<10
fname = ['sp13_0' num2str(jd) '_0' num2str(jd+1) ]
else
fname = ['sp13_' num2str(jd) '_' num2str(jd+1) ]
end
eval(['load ' fname])
f2=figure;
minspec=find(fxb>=0.005,1,'first')+2;
S=[];
for jp =1:9
    ss=subplot(3,3,jp);
    S=[S ss];
    plot(fm,Sn(:,jp),cl{3})
    hold on
    p1=plot(fxb(minspec:end),Snnf(minspec:end,nsx(jp)),cl{4});set(p1,'LineWidth',1.5);
    hold off
    title(['Gauge ',num2str(10*jp)]);
    ylabel('S_{\eta\eta} m^{2}/Hz')
    xlabel('F (Hz)')
    set(gca,'Xlim',[ 0 0.5])
    if jp<4
        set(gca,'ylim',[0 1]);
    elseif jp<7
        set(gca,'ylim',[0 2]);
    else
        set(gca,'ylim',[0 3.5]);
    end
    display(['Point ',num2str(jp)]);
    display(['hrms from whole spectrum = ',num2str(sqrt(8*trapz(fxb,Snn(:,nsx(jp)))))]);
    display(['hrms from water surface timeseries = ',num2str(sqrt(8)*std(detrend(squeeze(zsxbm(nsx(jp),:)))))]);
end
%  linkaxes(S);
% calculate low freq wave height 
fi = [0.005:0.005:0.05];
Sl = interp1(fm,Sn,fi);
Hrms_lo= sqrt(8*trapz(fi,Sl));   % assuming the transfer function equals one for infragravity waves

figure(f1)
s2=subplot(412);
plot(950-xmeter,Hrms_lo,cl{5},'markersize',5);
hold on
p1=plot(XB.x(:,41),hrmslo,cl{6});set(p1,'linewidth',1.5);
% p1=plot(XB.x(:,41),sqrt(8)*std(zsxbm,[],2),'r');set(p1,'linewidth',1.5);
hold off
set(gca,'ytick',[0 0.25 0.5]);
xlabel('x (m)')
ylabel('H_{rms,low} (m)')

% mean alongshore velocities

load vmean;


% add velocity
jm = find(x2>day+jd/24 & x2<day+(jd+1)/24)
jm = jm(1)

% linear interpolation for missing velocities

vmean3i(738:1016,3) = 0.5*(vmean3i(738:1016,2) + vmean3i(738:1016, 4));
vmean3i(360:500,7) = 0.5*(vmean3i(360:500,6) + vmean3i(360:500,8));
vmean3i(600:670,7) = 0.5*(vmean3i(600:670,5) + vmean3i(600:670,8));
vmean3i(500:1016,6) = 0.5*(vmean3i(500:1016,5) + vmean3i(500:1016, 7));

vm = 0.5*(vmean3i(jm,:)+vmean3i(jm+1,:))/100;
figure(f1)
s3=subplot(413);
plot(950-xmeter,vm,cl{7},'markersize',5);
hold on
p1=plot(XB.x(:,41),-vxbm(:,41),cl{8});set(p1,'linewidth',1.5);
hold off
xlabel('x (m)')
ylabel('V (m/s)')

s4=subplot(414);
zbmean0=readvar('zb_mean.dat',XB);
hold on
p1=plot(XB.x(:,41),zbmean0(:,41,1),cl{9});set(p1,'linewidth',2);
% p2=plot([0 838],[0 0],'b');
xlabel('x (m)')
ylabel('Elevation (m)')
plot(XB.x(nsx,41),0*ones(9,1),cl{10})
hold off
clear zbmean0

set([s1 s2 s3 s4],'xlim',[550 850]);
set(s1,'ylim',[0 2]);
set(s2,'ylim',[0 0.5]);
set(s3,'ylim',[-1.5 0.5]);


end