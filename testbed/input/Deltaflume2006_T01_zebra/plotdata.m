clear s m

%% coefficients and settings
fid = fopen(['datanam.txt']);
Tnam = fgetl(fid);
nvar = str2num(fgetl(fid));
for i = 1:nvar
    line = fgetl(fid);
    [nam,remain] = strtok(line,' ');
    [dummy,remain] = strtok(remain,'=');
    par.(nam) = str2num(remain(2:end));
end
fclose(fid);
%% measurements
load([Tnam]);
m=s;
clear s;

%% simulation results

% dimensions
fid = fopen('dims.dat','r'); 
temp = fread(fid,[7,1],'double'); 
nt = temp(2);
nx = temp(3)+1;
ny = temp(4)+1;
kmax = temp(5);
ngd = temp(6);
nd = temp(7);
fclose(fid);

% read grid coordinates 
fid = fopen('xy.dat','r');
xw = fread(fid,[nx,ny],'double'); 
yw = fread(fid,[nx,ny],'double');
fclose(fid);

x = [70 100 130 150 170:5:205]; % locations for q3d-flow output
ind{1} = 1:1:nx;
for i = 1:length(x)
    ind{2}(i) = find(xw(:,2)==x(i));
end
ind{3} = ind{2};

nam = [{'zb'};{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'cc'};{'ceq'};{'Su'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'Tbore'};{'Beta'};{'R'};{'DR'};{'kb'};{'c'};];
%nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'dzav'};{'z0bed'}]; %
nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'Fx'};{'dzav'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'BR'};{'Tbore'};{'R'};{'DR'};{'D'};{'kb'};{'kblong'};{'throl'};{'c'};{'gwhead'};{'gwlevel'};{'gwu'};{'gww'};]; %
%nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'ccg'};{'ceqg'};{'Sug'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'BR'};{'Tbore'};{'R'};{'DR'};{'D'};{'kb'};{'c'};]; %
nam2 = [{'ccg'};{'ceqsg'};{'ceqbg'};{'Susg'};{'Subg'};{'depo_ex'};{'ero'};];
nam3 = {'dzbed','sig'};                 dim3 = [nd kmax; length(ind{1}) length(ind{2});];
nam4 = {'pbbed','veloc','ccq3d'};       dim4 = [nd kmax kmax; ngd 3 ngd; length(ind{1}) length(ind{2}) length(ind{2});];

% read XBeach output
ts = 0:1:nt-1;
for j = 1:length(nam)
    temp = zeros(nx,ny,nt);
    fid = fopen([nam{j},'.dat'],'r');
    for i = 1:nt
        temp(:,:,i) = fread(fid,[nx,ny],'double');  % all data
    end
    fclose(fid);
    s.(nam{j}) = zeros(nx,nt);
    s.(nam{j}) = squeeze(temp(:,2,:));
end

for j = 1:length(nam2)
    temp = zeros(nx,ngd,nt);
    fid = fopen([nam2{j},'.dat'],'r');
    for i=1:nt
        for ii=1:ngd
            ttemp = fread(fid,[nx,ny],'double');
            temp(:,ii,i) = ttemp(:,2);               
        end
    end
    fclose(fid);
    s.(nam2{j}) = temp;
end

for j = 1:length(nam3)
    temp = zeros(dim3(2,j),dim3(1,j),nt);
    fid = fopen([nam3{j},'.dat'],'r');
    for i=1:nt
        for ii=1:dim3(1,j)
            ttemp = fread(fid,[nx,ny],'double');
            temp(:,ii,i) = ttemp(ind{j},2);               
        end
    end
    fclose(fid);
    s.(nam3{j}) = temp;
end

for j = 1:length(nam4)
    temp = zeros(dim4(3,j),dim4(1,j),dim4(2,j),nt);
    fid = fopen([nam4{j},'.dat'],'r');
    for i=1:nt
        for ii=1:dim4(2,j)
            for jj=1:dim4(1,j)
                ttemp = fread(fid,[nx,ny],'double');
                temp(:,jj,ii,i) = ttemp(ind{j},2);    
            end
        end
    end
    fclose(fid);
    s.(nam4{j}) = temp;
end

% ind = find(xw(:,2)==200);

% q3d-plot settings
fac1 = 2;
fac2 = 0.2;

temp2 = zeros(length(ind{2}),kmax+1,3,nt);
temp2(:,2:kmax+1,:,:) = s.veloc;
temp3 = zeros(length(ind{2}),kmax,ngd,nt);
temp3(:,2:kmax+1,:,:) = s.ccq3d*2650;
temp3(:,1,:,:) = temp3(:,2,:,:);
%

% movieplot

figure;
ca=[-0.2 1.2];
cmap=colormap;
c1r=interp1(ca,[1 length(cmap)],1);
c1=cmap(round(c1r),:);
c0r=interp1(ca,[1 length(cmap)],0);
c0=cmap(round(c0r),:);
i=1;
temp = zeros(length(xw),nd+1);
temp(:,2:end) = cumsum(s.dzbed(:,:,i),2);
zlayer = repmat(s.zb(:,i),1,nd+1)-temp;
xt = repmat(xw(:,2),1,nd+1);
temp(:,2:end) = s.pbbed(:,:,1,i);
temp(:,1) = temp(:,2);
% p1=    scatter(xt(:),zlayer(:),[],temp(:)); view(2); colorbar; caxis([0 1]);hold on;
p1=    pcolor(xt,zlayer,temp); view(2); colorbar; caxis(ca); hold on; 
p2=    plot(xw,s.zb(:,1),'k--','LineWidth',1.5); 
p3=    plot(xw,s.zb(:,i),'k-','LineWidth',1.5);
p4=    plot(xw,s.zs(:,i),'b','LineWidth',1.5);
p5=    plot(xw,2650*s.ccg(:,1,i),'color',c1,'LineWidth',1.5);
p6=    plot(xw,2650*s.ccg(:,2,i),'color',c0,'LineWidth',1.5);
p7=    plot(xw,s.H(:,i),'c','LineWidth',1.5);
p8=    plot(xw,s.ue(:,i),'g','LineWidth',1.5);
for j = 1:length(ind{2})    
   htemp = [0 (s.sig(j,:,i))*s.hh(ind{2}(j),i)];
   p9(j) = plot(xw(ind{2}(j),2)+fac1*(squeeze(temp2(j,:,1,i))+repmat(s.u(ind{2}(j),i),1,kmax+1)),s.zb(ind{2}(j),i)+htemp,'g-o');
   p10(j) = plot(xw(ind{2}(j),2)+fac1*squeeze(temp3(j,:,1,i)),s.zb(ind{2}(j),i)+htemp,'color',c1,'Marker','o');  
   p11(j) = plot(xw(ind{2}(j),2)+fac1*squeeze(temp3(j,:,2,i)),s.zb(ind{2}(j),i)+htemp,'color',c0,'Marker','o'); 
end
p12=  plot(xw,s.gwlevel(:,i),'y','LineWidth',1.5);  
p13 = plot(xw,s.zs(:,i)-s.throl(:,i),'b','LineWidth',1.5);  
%p9 =   plot(x,m.z(end,:),'k','LineWidth',1.5);
% axis([41 220 -4.5 2]);
axis([100 210 -2 2]);
for i = nt-100:nt; %nt-500:nt
    temp = zeros(length(xw),nd+1);
    temp(:,2:end) = cumsum(s.dzbed(:,:,i),2);
    zlayer = repmat(s.zb(:,i),1,nd+1)-temp;
    temp(:,2:end) = s.pbbed(:,:,1,i);
    temp(:,1) = temp(:,2);
    set(p1,'ydata',zlayer,'cdata',temp);
    % set(p1,'ydata',zlayer(:),'cdata',temp(:));
    set(p3,'ydata',s.zb(:,i));
    set(p4,'ydata',s.zs(:,i));
    set(p5,'ydata',2650*s.ccg(:,1,i));
    set(p6,'ydata',2650*s.ccg(:,2,i));
    set(p7,'ydata',s.H(:,i));
    set(p8,'ydata',s.ue(:,i));
    for j = 1:length(ind{2})    
       htemp = [0 (s.sig(j,:,i))*s.hh(ind{2}(j),i)];
       set(p9(j),'xdata',xw(ind{2}(j),2)+fac1*(squeeze(temp2(j,:,1,i))+repmat(s.u(ind{2}(j),i),1,kmax+1)));
       set(p9(j),'ydata',s.zb(ind{2}(j),i)+htemp); hold on;
       set(p10(j),'xdata',xw(ind{2}(j),2)+fac2*squeeze(temp3(j,:,1,i)));
       set(p10(j),'ydata',s.zb(ind{2}(j),i)+htemp); hold on;
       set(p11(j),'xdata',xw(ind{2}(j),2)+fac2*squeeze(temp3(j,:,2,i)));
       set(p11(j),'ydata',s.zb(ind{2}(j),i)+htemp); hold on;
    end
    set(p12,'ydata',s.gwlevel(:,i));
    set(p13,'ydata',s.zs(:,i)-s.throl(:,i));
    title(num2str(i));
    pause(0.1); hold off;
    
end
print final_profile.png -dpng

[X,T] = meshgrid(xw(:,2),ts);
range = [1 ts(end) 100 210];
figure;
subplot(311);
pcolor(T,X,s.kblong'); shading flat; caxis([0 1]); colorbar; 
axis(range);
subplot(312);
pcolor(T,X,s.kb'+s.kblong'); shading flat; caxis([0 1]); colorbar; 
axis(range);
subplot(313);
% pcolor(T,X,2650*squeeze(sum(s.ccg,2))'); shading flat; colorbar;
pcolor(T,X,log(2650*squeeze(sum(s.ccg,2))')); shading flat; colorbar;
axis(range);

% setup
s.setup = mean(s.zs-max(0,s.zb),2);

% long wave height
h = s.zs-s.zb;
windowSize = 60; % ten seconds window 
fs = fspecial('average',windowSize);
hm = [];
for i=1:length(xw)
    hm(i,:) = imfilter(h(i,:),fs,'symmetric');
end
zs = h-hm;

s.Hrms_lfm = sqrt(8)*std(zs')';
s.Hrms_lfm2 = sqrt(8)*std(h')';
s.Hrms_hfm = sqrt(mean(s.Hrms.^2,2));
% s.Hrms_lfm = sqrt(8)*std(detrend((s.zs-s.zb)'))';
s.Hrms_m = sqrt(s.Hrms_hfm.^2 + s.Hrms_lfm.^2);
s.Urms_hfm = sqrt(mean(s.urms.^2,2));
s.Urms_lfm = std(s.u')';
s.Urms_m = sqrt(s.Urms_hfm.^2 + s.Urms_lfm.^2);
s.Uon_m = 1/sqrt(2)*mean(s.uon,2);
s.Uoff_m = 1/sqrt(2)*mean(s.uoff,2);
s.Um = mean(s.ue,2);
s.Ua = mean(s.ua,2);
s.Skm = mean(s.Sk,2);
s.Asm = mean(s.As,2);
s.Cm = mean(s.cc,2);
s.Ceqm = mean(s.ceq,2);
s.usR = mean(2*s.R./(par.rho.*max(h,0.1).*s.c),2);
s.usW = mean(0.125*par.g*s.Hrms.^2./(max(h,0.1).*s.c),2);

% transform Su to Sz
Stemp = zeros(size(s.Su));
Stemp(1,:) = s.Su(1,:);
Stemp(2:end,:) = s.Su(1:end-1,:); %0.5*(s.Su(1:end-1,:)+s.Su(2:end,:));

% compare measured and simulated profiles
int = [];
for i =1:length(m.pint) int(i) = str2num(m.pint{i}); end;
ind = find(int<nt*par.morfac/3600);
Tout = round(int(ind)*3600/par.morfac/par.tint)+1;

dz = s.zb(:,Tout(ind(end)))-s.zb(:,1);
Tsim = Tout(end)*par.morfac;
s.Sdzbtot = (1-par.np)*flipud(cumsum(flipud(dz)))*(xw(2,2)-xw(1,2))/Tsim;
s.Sdzavtot = (1-par.np)*flipud(cumsum(flipud(s.dzav(:,end))))*(xw(2,2)-xw(1,2))/Tsim;
s.Scutot = mean(Stemp,2);

% Analysis long waves
w = 2*pi/par.Tp;
hh = max(mean(h,2),0.01); % max(h,0.01); 
k = disper(w,hh,par.g);
cg = w./k*0.5.*(1+2*k.*hh./sinh(2*k.*hh));
cg = sqrt(par.g*hh);
% cg = sqrt(par.g*hh);
% Guza
s.etain=zeros(nx,nt);
s.etaout=zeros(nx,nt);
s.uin=zeros(nx,nt);
for i = 1:nx
    s.etain(i,:) = (detrend(zs(i,:)).*sqrt(9.81.*hh(i,:))+detrend(s.u(i,:)).*hh(i,:))./(cg(i,:)+sqrt(9.81.*hh(i,:)));
    s.etaout(i,:) = (detrend(zs(i,:)).*sqrt(9.81.*hh(i,:))-detrend(s.u(i,:)).*hh(i,:))./(cg(i,:)+sqrt(9.81.*hh(i,:)));
    s.uin(i,:) = par.g./cg(i,:).*s.etain(i,:);
end
s.Hrms_lfm_in = sqrt(8)*std(s.etain')';
s.Hrms_lfm_out = sqrt(8)*std(s.etaout')';
% phase averaged work on long waves (Battjes, 2004)
s.RR = -mean(s.uin.*s.Fx,2);
% correlation short and longh waves
s.rho = []; s.pval = [];
for i = 1:nx
    [s.rho(i), s.pval(i)] = corr(detrend(s.zs(i,:))',s.Hrms(i,:).^2');
end
%% figures;

% hydrodynamics
figure(1);
subplot(221);
plot(xw,s.Hrms_m,'k-','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.Hrms.^2)),'ks','LineWidth',1.5);
plot(xw,s.Hrms_hfm,'k--','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.Hrms_hf.^2)),'k^','LineWidth',1.5);
plot(xw,s.Hrms_lfm,'k-.','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.Hrms_lf.^2)),'kv','LineWidth',1.5);
xlabel('x [m]'); ylabel('H_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.25]);
subplot(222);
plot(xw,s.setup,'k-','LineWidth',1.5); hold on; plot(m.xps,mean(m.setup),'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('\eta_{mean} [m]'); axis([min(min(xw)) max(max(xw)) -0.1 0.3]);
subplot(223);
plot(xw,s.Urms_m,'k-','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.urms1.^2)),'ks','LineWidth',1.5); % plot(m.xm1,m.urms3(:,1),'ys','LineWidth',1.5); 
plot(xw,s.Urms_hfm,'k--','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.urms_hf1.^2)),'k^','LineWidth',1.5); % plot(m.xm1,m.urms_hf3(:,1),'y^','LineWidth',1.5); 
plot(xw,s.Urms_lfm,'k-.','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.urms_lf1.^2)),'kv','LineWidth',1.5); % plot(m.xm1,m.urms_lf3(:,1),'yv','LineWidth',1.5);
xlabel('x [m]'); ylabel('U_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
subplot(224);
plot(xw,s.Um,'k-','LineWidth',1.5); hold on; plot(m.xm2,m.um,'ks','LineWidth',1.5);
plot(xw,s.Ua,'r-','LineWidth',1.5); plot(m.xps,0.25*mean(m.ua),'rs','LineWidth',1.5);
plot(xw,s.Um+s.Ua,'b','LineWidth',1.5); 
% plot(xw(:,2),s.usW,'g--','LineWidth',1.5); plot(xw(:,2),s.usR,'g-','LineWidth',1.5); plot(xw(:,2),-(s.usR+s.usW),'g-.','LineWidth',1.5);
xlabel('x [m]'); ylabel('U_{mean} & U_{wave assymetry} [m/s]'); axis([min(min(xw)) max(max(xw)) -0.75 0.75]);
print('hydrodynamics.png','-dpng');

% morphodynamics
figure(2);
subplot(221);
plot(xw,s.Cm*par.rhos,'k-','LineWidth',1.5); hold on; plot(m.xm2,m.cm,'ks','LineWidth',1.5); 
plot(xw,s.Ceqm*par.rhos,'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('C_{mean} [g/l]'); axis([min(min(xw)) max(max(xw)) 0 30]);
subplot(222);
plot(xw(:,2),s.Sdzavtot,'r--','LineWidth',1.5); hold on; plot(xw(:,2),s.Scutot,'r-.','LineWidth',1.5); 
plot(xw(:,2),s.Sdzbtot,'k--','LineWidth',1.5); plot(m.x{1},m.Sdzb{ind(end)},'k-','LineWidth',1.5);
xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]'); axis([0 max(max(xw)) -3E-4 1E-4]);
legend('S_{Avalanch}','S_{CU}','S_{total,simulated}','S_{total,measured}',3);
subplot(2,2,3:4);
for i = 1:length(Tout)
    plot(xw(:,2),s.zb(:,Tout(i)),'r--','LineWidth',1.5); hold on; plot(m.x{i},m.z{i},'k','LineWidth',1.5);
end
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([100 max(max(xw)) -2.5 2]);
print('morphodynamics.png','-dpng');

% Detailed wave analysis
figure(3);
subplot(221);
plot(xw(:,2),s.Hrms_lfm,'k','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.Hrms_lf.^2)),'kv','LineWidth',1.5);
plot(xw(:,2),s.Hrms_lfm_in,'b','LineWidth',1.5); 
plot(xw(:,2),s.Hrms_lfm_out,'b--','LineWidth',1.5);
plot(xw(:,2),sqrt(s.Hrms_lfm_in.^2+s.Hrms_lfm_out.^2),'r--','LineWidth',1.5);
% plot Greens' curve: h^-0.25 for outgoing waves
plot(xw(:,2),s.Hrms_lfm_out(1)*hh.^-0.25/hh(1).^-0.25,'b--');
plot(xw(:,2),s.Hrms_lfm_in(1)*hh.^-2.50/hh(1).^-2.50,'b');
legend('Long wave height','measurements','Incoming long wave height','Reflected long wave height','Green: H ~ h^{-1/4}','Longuet-Higgins & Stewart: H ~ h^{-5/2}',2);
xlabel('x [m]'); ylabel('H_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
subplot(222);
plot(xw,s.rho,'k','LineWidth',1.5); hold on; plot(m.xps,mean(m.rho),'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('\rho [-]'); axis([min(min(xw)) max(max(xw)) -1 1]);
subplot(223);
plot(xw,s.RR,'k','LineWidth',1.5);
xlabel('x [m]'); ylabel('R [W/m^2]'); axis([min(min(xw)) max(max(xw)) 1.25*min(s.RR) 1.25*max(s.RR)]);
subplot(224);
Beta = atan(s.Asm./s.Skm); B = sqrt(s.Asm.^2+s.Skm.^2);
BetaFW = mean(atan(m.As1./m.Sk1)); BFW = mean(sqrt(m.As1.^2+m.Sk1.^2));
% BetaMF = atan(1/3*(AsuMF(:,1)+AsuMF(:,2)+AsuMF(:,3))./(1/3*(SkuMF(:,1)+SkuMF(:,2)+SkuMF(:,3))));
% BMF = sqrt(1/3*(AsuMF(:,1)+AsuMF(:,2)+AsuMF(:,3)).^2+(1/3*(SkuMF(:,1)+SkuMF(:,2)+SkuMF(:,3))).^2);
plot(xw(:,2),s.Skm,'b','LineWidth',1.5); hold on; plot(m.xps,m.Sk1,'bs','LineWidth',1.5); 
plot(xw(:,2),s.Asm,'r','LineWidth',1.5); plot(m.xps,m.As1,'rs','LineWidth',1.5); 
plot(xw(:,2),Beta,'g','LineWidth',1.5); plot(m.xps,BetaFW,'gs','LineWidth',1.5);
plot(xw(:,2),B,'k','LineWidth',1.5); plot(m.xps,BFW,'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('S_K & A_S [-]'); axis([min(min(xw)) max(max(xw)) -2 2]);
legend('Skewness','Assymetry',3);
print('detailed_waves.png','-dpng');

% detailed analysis sediment suspensions
figure(4)
subplot(221); plot(xw(:,2),mean(s.Tbore,2)/par.Tp,'k','LineWidth',1.5); 
axis([min(min(xw)) max(max(xw)) 0 1]); xlabel('x [m]'); ylabel('T_{BORE}/T_P [-]');
subplot(222); plot(xw(:,2),mean(s.Beta,2),'k','LineWidth',1.5); hold on; plot(xw(:,2),mean(tan(asin(s.Beta)),2),'r','LineWidth',1.5);
axis([min(min(xw)) max(max(xw)) 0 2.0]); xlabel('x [m]'); ylabel('sin(\beta) [-] & \partial\eta/\partialx_{max} [-]');
subplot(223); plot(xw(:,2),mean(2*s.R.*par.Tp./max(0.01,(par.rho.*s.c.*s.Hrms)),2),'k','LineWidth',1.5); hold on;
plot(xw(:,2),mean(max(0.01,s.Hrms),2),'r','LineWidth',1.5); hold on;
axis([min(min(xw)) max(max(xw)) 0 3]); xlabel('x [m]'); ylabel('mixing length [m]');
subplot(224); plot(xw(:,2),mean((s.DR/par.rho).^(2/3).*par.Tp./s.Tbore,2),'r','LineWidth',1.5); hold on;
              plot(xw(:,2),mean(s.kb,2),'b','LineWidth',1.5);
axis([min(min(xw)) max(max(xw)) 0 4]); xlabel('x [m]'); ylabel('k & k_b [m^2/s^2]');
print('detailed_sediment_suspensions.png','-dpng');

fclose all; close all; 
