%% coefficients and settings
fid = fopen('..\data\datanam.txt');
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
load(['..\data\',Tnam]);

%% simulation results
nam = [{'zb'};{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'cc'};{'ceq'};{'Su'}; {'Fx'};];

% dimensions
fid = fopen('dims.dat','r'); 
temp = fread(fid,[3,1],'double'); 
nt = temp(1);
nx = temp(2)+1;
ny = temp(3)+1;
fclose(fid);

% read grid coordinates 
fid = fopen('xy.dat','r');
xw = fread(fid,[nx,ny],'double'); 
yw = fread(fid,[nx,ny],'double');
% x = fread(fid,[nx,ny],'double'); 
% y = fread(fid,[nx,ny],'double');
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
s.Um = mean(s.ue,2);
s.Cm = mean(s.cc,2);
dz = s.zb(:,end)-s.zb(:,1);
Tsim = nt*par.morfac*par.tint;
s.Sdzb = (1-par.np)*flipud(cumsum(flipud(dz)))*(xw(2,2)-xw(1,2))/Tsim;
Tout = round([0 0.1 0.3 1.0 2.04 6.0]*3600/par.morfac/par.tint)+1;
s.zb_interval = s.zb(:,Tout);

% Analysis long waves
w = 2*pi/par.Tp;
hh = max(mean(h,2),0.01); % max(h,0.01); 
k = disper(w,hh,par.g);
% cg = w./k*0.5.*(1+2*k.*hh./sinh(2*k.*hh));
% h0 = setup(1)-s.zb(:,2);
cg = sqrt(par.g*hh);
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
s.R = -mean(s.uin.*s.Fx,2);
% correlation short and longh waves
s.rho = []; s.pval = [];
for i = 1:nx
    [s.rho(i), s.pval(i)] = corr(detrend(s.zs(i,:))',s.Hrms(i,:).^2');
end
%% figures;

% hydrodynamics
figure(1);
subplot(221);
plot(xw,s.Hrms_m,'k-','LineWidth',1.5); hold on; plot(xhrms,Hrms,'ks','LineWidth',1.5);
plot(xw,s.Hrms_hfm,'k--','LineWidth',1.5); hold on; plot(xhrms,Hrms_hf,'k^','LineWidth',1.5);
plot(xw,s.Hrms_lfm,'k-.','LineWidth',1.5); hold on; plot(xhrms,Hrms_lf,'kv','LineWidth',1.5);
xlabel('x [m]'); ylabel('H_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.25]);
subplot(222);
plot(xw,s.setup,'k-','LineWidth',1.5); hold on; plot(xsetup,setup,'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('\eta_{mean} [m]'); axis([min(min(xw)) max(max(xw)) -0.1 0.3]);
subplot(223);
plot(xw,s.Urms_m,'k-','LineWidth',1.5); hold on; plot(xurms,Urms,'ks','LineWidth',1.5); 
plot(xw,s.Urms_hfm,'k--','LineWidth',1.5); hold on; plot(xurms,Urms_hf,'k^','LineWidth',1.5);
plot(xw,s.Urms_lfm,'k-.','LineWidth',1.5); hold on; plot(xurms, Urms_lf,'kv','LineWidth',1.5);
xlabel('x [m]'); ylabel('U_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
subplot(224);
plot(xw,s.Um,'k-','LineWidth',1.5); hold on; plot(xum,Um,'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('U_{mean} [m]'); axis([min(min(xw)) max(max(xw)) -0.75 0.1]);
print('hydrodynamics.png','-dpng');

% morphodynamics
figure(2);
subplot(221);
plot(xw,s.Cm*par.rhos,'k-','LineWidth',1.5); hold on; plot(xc,Cm,'ks','LineWidth',1.5); 
xlabel('x [m]'); ylabel('C_{mean} [g/l]'); axis([min(min(xw)) max(max(xw)) 0 30]);
subplot(222);
plot(xw,s.Sdzb,'k--','LineWidth',1.5); hold on; plot(x,Sdzb,'k-','LineWidth',1.5);
xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]'); axis([min(min(xw)) max(max(xw)) -8E-4 1E-4]);
subplot(2,2,3:4);
for i = 1:length(Tout)
    plot(xw,s.zb_interval(:,i),'r--','LineWidth',1.5); hold on; plot(x,z(i,:),'k-','LineWidth',1.5);
end
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([150 max(max(xw)) -2.5 2]);
print('morphodynamics.png','-dpng');

% incoming versus outgoing long waves
figure(3);
subplot(221);
plot(xw,s.Hrms_lfm,'k','LineWidth',1.5); hold on; plot(xhrms,Hrms_lf,'kv','LineWidth',1.5);
plot(xw,s.Hrms_lfm_in,'b','LineWidth',1.5); 
plot(xw,s.Hrms_lfm_out,'b--','LineWidth',1.5);
plot(xw,sqrt(s.Hrms_lfm_in.^2+s.Hrms_lfm_out.^2),'r--','LineWidth',1.5);
xlabel('x [m]'); ylabel('H_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 0.75]);
subplot(222);
plot(xw,s.rho,'k','LineWidth',1.5); 
xlabel('x [m]'); ylabel('\rho [-]'); axis([min(min(xw)) max(max(xw)) -1 1]);
subplot(223);
plot(xw,s.R,'k','LineWidth',1.5);
xlabel('x [m]'); ylabel('R [W/m^2]'); axis([min(min(xw)) max(max(xw)) 1.25*min(s.R) 1.25*max(s.R)]);
print('long_waves.png','-dpng');

pause(1.0)

fclose all; close all;