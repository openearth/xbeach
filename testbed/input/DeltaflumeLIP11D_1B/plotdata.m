function plotdata();

clear s m

%% coefficients and settings
[runid,testid,datadir]=testinfo
fid = fopen([datadir 'datanam.txt']);
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
load([datadir,Tnam]);
%% measurements
load([datadir,Tnam]);
m=s;
clear s;

%% simulation results
nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'As'};{'Sk'};{'Tbore'};{'R'};{'DR'};{'kb'};{'c'};{'BR'};]; % 
nam2 = [{'ccg'};{'ceqbg'};{'ceqsg'};{'Subg'};{'Susg'};];

% dimensions
fid = fopen('dims.dat','r'); 
temp = fread(fid,[7,1],'double'); 
nt = temp(1);
nx = temp(2)+1;
ny = temp(3)+1;
ntheta = temp(4);
kmax = temp(5);
ngd = temp(6);
nd = temp(7);
fclose(fid);

% read grid coordinates 
fid = fopen('xy.dat','r');
xw = fread(fid,[nx,ny],'double'); 
yw = fread(fid,[nx,ny],'double');
fclose(fid);

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
s.Hrms_hfm = sqrt(mean(s.H.^2,2));
s.Hrms_m = sqrt(s.Hrms_hfm.^2 + s.Hrms_lfm.^2);
s.Urms_hfm = sqrt(mean(s.urms.^2,2));
s.Urms_lfm = std(s.u')';
s.Urms_m = sqrt(s.Urms_hfm.^2 + s.Urms_lfm.^2);
s.Um = mean(s.ue,2);
s.Ua = mean(s.ua,2);
s.Skm = mean(s.Sk,2);
s.Asm = mean(s.As,2);
s.Cm = mean(squeeze(s.ccg),2);
s.ceq = squeeze(s.ceqsg);
s.Ceqsm = mean(s.ceq,2);
s.Ceqbm = mean(squeeze(s.ceqbg),2);
s.Su = squeeze(s.Susg) + squeeze(s.Subg);
s.usR = mean(2*s.R./(par.rho.*max(h,0.1).*s.c),2);
s.usW = mean(0.125*par.g*s.H.^2./(max(h,0.1).*s.c),2);

% transform Su to Sz
Ssus1 = zeros(size(s.Susg));
Ssus1(1,:,:) = s.Susg(1,:,:);
Ssus1(2:end,:,:) = s.Susg(1:end-1,:,:); %0.5*(s.Su(1:end-1,:)+s.Su(2:end,:));
Sbed1 = zeros(size(s.Subg));
Sbed1(1,:,:) = s.Subg(1,:,:);
Sbed1(2:end,:,:) = s.Subg(1:end-1,:,:); %0.5*(s.Su(1:end-1,:)+s.Su(2:end,:));

% compare measured and simulated profiles
int = [];
for i =1:length(m.pint) int(i) = str2num(m.pint{i}); end;
ind = find(int<nt*par.morfac/3600);
Tout = round(int(ind)*3600/par.morfac/par.tint)+1;

dz = s.zb(:,Tout(ind(end)))-s.zb(:,1);
Tsim = Tout(end)*par.morfac;
s.Sdzbtot = (1-par.np)*flipud(cumsum(flipud(dz)))*(xw(2,2)-xw(1,2))/Tsim;
s.Sdzavtot = (1-par.np)*flipud(cumsum(flipud(s.dzav(:,end))))*(xw(2,2)-xw(1,2))/Tsim;
s.Sustot = mean(Ssus1,3);
s.Subtot = mean(Sbed1,3);
s.Sutot = s.Subtot+s.Sustot

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
% s.rho = []; s.pval = [];
% for i = 1:nx
%     [s.rho(i), s.pval(i)] = corr(detrend(s.zs(i,:))',s.Hrms(i,:).^2');
% end
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
pname = ['..\..\report\',testid '_' runid '_hydrodynamics' '.jpg'];
eval(['print -djpeg ' pname]);

% morphodynamics
figure(2);
subplot(221);
plot(xw,s.Cm*par.rhos,'k-','LineWidth',1.5); hold on; plot(m.xm2,m.cm,'ks','LineWidth',1.5); 
plot(xw,s.Ceqsm*par.rhos,'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('C_{mean} [g/l]'); axis([0 max(max(xw)) 0 30]);
subplot(222);
plot(xw(:,2),s.Sdzbtot,'k--','LineWidth',1.5); hold on; plot(m.x{1},m.Sdzb{ind(end)},'k-','LineWidth',1.5);
plot(xw(:,2),s.Sdzavtot,'r--','LineWidth',1.5); plot(xw(:,2),s.Subtot,'r:','LineWidth',1.5); plot(xw(:,2),s.Sustot,'r-.','LineWidth',1.5); plot(xw(:,2),s.Sdzavtot+s.Sutot,'r');
xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]'); axis([100 max(max(xw)) -1E-4 1E-4]);
legend('S_{\Deltaz,simulated}','S_{\Deltaz,measured}','S_{Avalanching}','S_{bed}','S_{suspended}','S_{total,simulated}',2);
subplot(2,2,3:4);
for i = 1:length(Tout)
    plot(xw(:,2),s.zb(:,Tout(i)),'r--','LineWidth',1.5); hold on; plot(m.x{i},m.z{i},'k','LineWidth',1.5);
end
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([100 max(max(xw)) -2.5 2]);
pname = ['..\..\report\',testid '_' runid '_morphodynamics' '.jpg'];
eval(['print -djpeg ' pname]);

% Detailed wave analysis
figure(3);
plot(xw(:,2),s.Hrms_lfm,'k','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.Hrms_lf.^2)),'kv','LineWidth',1.5);
plot(xw(:,2),s.Hrms_lfm_in,'b','LineWidth',1.5); 
plot(xw(:,2),s.Hrms_lfm_out,'b--','LineWidth',1.5);
plot(xw(:,2),sqrt(s.Hrms_lfm_in.^2+s.Hrms_lfm_out.^2),'r--','LineWidth',1.5);
% plot Greens' curve: h^-0.25 for outgoing waves
plot(xw(:,2),s.Hrms_lfm_out(1)*hh.^-0.25/hh(1).^-0.25,'b--');
plot(xw(:,2),s.Hrms_lfm_in(1)*hh.^-2.50/hh(1).^-2.50,'b');
legend('Long wave height','measurements','Incoming long wave height','Reflected long wave height','Green: H ~ h^{-1/4}','Longuet-Higgins & Stewart: H ~ h^{-5/2}',2);
xlabel('x [m]'); ylabel('H_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
pname = ['..\..\report\',testid '_' runid '_detailed_waves' '.jpg'];
eval(['print -djpeg ' pname]);

end