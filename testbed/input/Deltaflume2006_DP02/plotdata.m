clear s

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
nam = [{'zb'};{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'cc'};{'ceq'};{'Su'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'Tbore'};{'Beta'};{'R'};{'DR'};{'kb'};{'c'};{'Ts'};];

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
s.Uon_m = 1/sqrt(2)*mean(s.uon,2);
s.Uoff_m = 1/sqrt(2)*mean(s.uoff,2);
s.Um = mean(s.ue,2);
s.Ua = mean(s.ua,2);
s.Skm = mean(s.Sk,2);
s.Asm = mean(s.As,2);
s.Cm = mean(s.cc,2);
s.Ceqm = mean(s.ceq,2);
dz = s.zb(:,end)-s.zb(:,1);
s.usR = mean(2*s.R./(par.rho.*max(h,0.1).*s.c),2);
s.usW = mean(0.125*par.g*s.Hrms.^2./(max(h,0.1).*s.c),2);

% transform Su to Sz
Stemp = zeros(size(s.Su));
Stemp(1,:) = s.Su(1,:);
Stemp(2:end,:) = s.Su(1:end-1,:); %0.5*(s.Su(1:end-1,:)+s.Su(2:end,:));

Tsim = nt*par.morfac*par.tint;
s.Sdzbtot = (1-par.np)*flipud(cumsum(flipud(dz)))*(xw(2,2)-xw(1,2))/Tsim;
s.Sdzavtot = (1-par.np)*flipud(cumsum(flipud(s.dzav(:,end))))*(xw(2,2)-xw(1,2))/Tsim;
s.Scutot = mean(Stemp,2);

% Get interval averaged values
Tout = round([0 0.1 0.3 1.0 2.04 6.0]*3600/par.morfac/par.tint)+1;
s.zb_interval = s.zb(:,Tout);
for i = 1:length(Tout)-1
    
    % ratio long wave height and short wave height
    s.Hrmslf(:,i) = sqrt(8)*std(zs(:,Tout(i):Tout(i+1))')';
    s.Hrmshf(:,i) = mean(s.Hrms(:,Tout(i):Tout(i+1)).^2,2);
    % ratio long wave orbital motion and short wave orbital motion
    s.Urmslf(:,i) = std(s.u(:,Tout(i):Tout(i+1))')';
    s.Urmshf(:,i) = mean(s.urms(:,Tout(i):Tout(i+1)).^2,2);
    % mean flow
    s.Um2(:,i) = mean(s.ue(:,Tout(i):Tout(i+1)),2);
    % mean sediment suspensions
    s.Cm2(:,i) = mean(s.cc(:,Tout(i):Tout(i+1)),2);
    % ration acual concentration and equilibrium concentration
    s.ccceq(:,i) = mean(s.cc(:,Tout(i):Tout(i+1)),2)./mean(s.ceq(:,Tout(i):Tout(i+1)),2);
    % total sediment transport from sediment suspensions and flow including bedslope terms and diffusion
    s.Scu(:,i) = mean(Stemp(:,Tout(i):Tout(i+1)),2);
    % total sediment transport from sediment suspensions and flow
    s.Scu2(:,i) = mean(s.ue(:,Tout(i):Tout(i+1)).*s.cc(:,Tout(i):Tout(i+1)).*s.hh(:,Tout(i):Tout(i+1)),2);
    % sediment transport due to bedslope effects
    
    % sediment transport due to diffusion
    
    % mean sediment transport from sediment suspensions and flow
    s.Scum(:,i) = s.Um2(:,i).*s.Cm2(:,i).*mean(s.hh(:,Tout(i):Tout(i+1)),2);
    % long wave related sediment transport from sediment suspensions and flow
    s.Scuf(:,i) = s.Scu2(:,i)-s.Scum(:,i);
    dz = s.zb(:,Tout(i+1))-s.zb(:,Tout(i));
    dzav = s.dzav(:,Tout(i+1))-s.dzav(:,Tout(i));
    % sediment transport from bed level changes
    s.Sdzb(:,i) = (1-par.np)*flipud(cumsum(flipud(dz)))*(xw(2,2)-xw(1,2))/((Tout(i+1)-Tout(i))*par.morfac*par.tint);
    % sediment transport due to avalanching
    s.Sdzav(:,i) = (1-par.np)*flipud(cumsum(flipud(dzav)))*(xw(2,2)-xw(1,2))/((Tout(i+1)-Tout(i))*par.morfac*par.tint);
    
end

% compute sediment supply from dune face as function of time;
Vflow(1)=0;
for i = 1:nt;
    ind = find(s.dzav(:,i)<0);
    Vav(i) = (1-par.np)*sum(s.dzav(ind,i))*(xw(2,2)-xw(1,2));
    ind = find(abs(s.Su(:,i))==max(abs(s.Su(:,i))));
    Vflow(i) = Vflow(max(1,i-1))+s.Su(ind,i);
end

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
% plot(xw,s.Uon_m,'b--','LineWidth',1.5); hold on; plot(xw,s.Uoff_m,'r--','LineWidth',1.5);    
plot(xw,s.Urms_lfm,'k-.','LineWidth',1.5); hold on; plot(xurms, Urms_lf,'kv','LineWidth',1.5);
xlabel('x [m]'); ylabel('U_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
subplot(224);
plot(xw,s.Um,'k-','LineWidth',1.5); hold on; plot(xum,Um,'ks','LineWidth',1.5);
plot(xw,s.Ua,'r-','LineWidth',1.5); plot(xhrms,0.25*ua,'rs','LineWidth',1.5); plot(xurms,0.25*1/3*(uaMF(:,1)+uaMF(:,2)+uaMF(:,3)),'rv','LineWidth',1.5);
plot(xw,s.Um+s.Ua,'b','LineWidth',1.5); 
plot(xw(:,2),s.usW,'g--','LineWidth',1.5); plot(xw(:,2),s.usR,'g-','LineWidth',1.5); plot(xw(:,2),-(s.usR+s.usW),'g-.','LineWidth',1.5);
xlabel('x [m]'); ylabel('U_{mean} & U_{wave assymetry} [m/s]'); axis([min(min(xw)) max(max(xw)) -0.75 0.75]);
print('hydrodynamics.png','-dpng');
% morphodynamics
figure(2);
subplot(221);
plot(xw,s.Cm*par.rhos,'k-','LineWidth',1.5); hold on; plot(xc,Cm,'ks','LineWidth',1.5); 
plot(xw,s.Ceqm*par.rhos,'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('C_{mean} [g/l]'); axis([min(min(xw)) max(max(xw)) 0 30]);
subplot(222);
plot(xw,s.Sdzbtot,'k--','LineWidth',1.5); hold on; plot(x,Sdzb,'k-','LineWidth',1.5);
plot(xw,s.Sdzavtot,'r--','LineWidth',1.5); plot(xw,s.Scutot,'r-.','LineWidth',1.5); plot(xw,s.Sdzav(:,end)+s.Scutot,'r','LineWidth',1.5);
xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]'); axis([0 max(max(xw)) -8E-4 1E-4]);
subplot(2,2,3:4);
for i = 1:length(Tout)
    plot(xw,s.zb_interval(:,i),'r--','LineWidth',1.5); hold on; plot(x,z(i,:),'k-','LineWidth',1.5);
end
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([100 max(max(xw)) -2.5 2]);
print('morphodynamics.png','-dpng');

% Detailed wave analysis
figure(3);
subplot(221);
plot(xw(:,2),s.Hrms_lfm,'k','LineWidth',1.5); hold on; plot(xhrms,Hrms_lf,'kv','LineWidth',1.5);
plot(xw(:,2),s.Hrms_lfm_in,'b','LineWidth',1.5); 
plot(xw(:,2),s.Hrms_lfm_out,'b--','LineWidth',1.5);
plot(xw(:,2),sqrt(s.Hrms_lfm_in.^2+s.Hrms_lfm_out.^2),'r--','LineWidth',1.5);
% plot Greens' curve: h^-0.25 for outgoing waves
plot(xw(:,2),s.Hrms_lfm_out(1)*hh.^-0.25/hh(1).^-0.25,'b--');
plot(xw(:,2),s.Hrms_lfm_in(1)*hh.^-2.50/hh(1).^-2.50,'b');
legend('Long wave height','measurements','Incoming long wave height','Reflected long wave height','Green: H ~ h^{-1/4}','Longuet-Higgins & Stewart: H ~ h^{-5/2}',2);
xlabel('x [m]'); ylabel('H_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
subplot(222);
plot(xw,s.rho,'k','LineWidth',1.5); hold on; plot(xhrms,rho,'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('\rho [-]'); axis([min(min(xw)) max(max(xw)) -1 1]);
subplot(223);
plot(xw,s.RR,'k','LineWidth',1.5);
xlabel('x [m]'); ylabel('R [W/m^2]'); axis([min(min(xw)) max(max(xw)) 1.25*min(s.RR) 1.25*max(s.RR)]);
subplot(224);
Beta = atan(s.Asm./s.Skm); B = sqrt(s.Asm.^2+s.Skm.^2);
BetaFW = atan(Asu./Sku); BFW = sqrt(Asu.^2+Sku.^2);
BetaMF = atan(1/3*(AsuMF(:,1)+AsuMF(:,2)+AsuMF(:,3))./(1/3*(SkuMF(:,1)+SkuMF(:,2)+SkuMF(:,3))));
BMF = sqrt(1/3*(AsuMF(:,1)+AsuMF(:,2)+AsuMF(:,3)).^2+(1/3*(SkuMF(:,1)+SkuMF(:,2)+SkuMF(:,3))).^2);
plot(xw(:,2),s.Skm,'b','LineWidth',1.5); hold on; plot(xhrms,Sku,'bs','LineWidth',1.5); plot(xurms,1/3*(SkuMF(:,1)+SkuMF(:,2)+SkuMF(:,3)),'rv','LineWidth',1.5);
plot(xw(:,2),s.Asm,'r','LineWidth',1.5); plot(xhrms,Asu,'rs','LineWidth',1.5); plot(xurms,1/3*(AsuMF(:,1)+AsuMF(:,2)+AsuMF(:,3)),'rv','LineWidth',1.5);
plot(xw(:,2),Beta,'g','LineWidth',1.5); plot(xhrms,BetaFW,'gs','LineWidth',1.5); plot(xurms,BetaMF,'gv','LineWidth',1.5);
plot(xw(:,2),B,'k','LineWidth',1.5); plot(xhrms,BFW,'ks','LineWidth',1.5); plot(xurms,BMF,'kv','LineWidth',1.5);
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

% sediment transports in XBeach
figure(5);
for i = 1:length(Tout)-1
    subplot(3,2,i)
    % transports from flow and sediment suspension
    plot(xw(:,2),s.Scu2(:,i),'g','LineWidth',2); hold on;
    plot(xw(:,2),s.Scum(:,i),'g--','LineWidth',2); plot(xw(:,2),s.Scuf(:,i),'g-.','LineWidth',2);
    % Sbedslopes/diffusion = Scu-Scu2
    plot(xw(:,2),s.Scu(:,i)-s.Scu2(:,i),'r--');
    % Sdzb = Sdzav+Scu
    plot(xw(:,2),s.Scu(:,i),'b','LineWidth',2);
     % plot(xw(:,2),s.Sdzav(:,i)+s.Scu(:,i),'r','LineWidth',2); 
     % transports from bed level changes
    plot(xw(:,2),s.Sdzb(:,i),'k','LineWidth',1.5); plot(xw,s.Sdzav(:,i),'k--','LineWidth',1.5); 
    xlabel('x [m]'); ylabel('S [m^3/m/s]');
    axis([150 max(max(xw)) 1.25*min(min(s.Scu2(:,i)),min(s.Sdzb(:,i))) 1.25*max(max(s.Scuf(:,i)),max(s.Scu(:,i)-s.Scu2(:,i)))]);
    if i==1
        legend('S = U*C','S = U_{mean}*C_{mean}','S = U_{long wave}*C_{long wave}','S_{XBeach} - U*C','S_{XBeach}','S_{\Deltazb}','S_{avalanch}',3);
    end
end
print('sediment_transports_intervals.png','-dpng');

% plot test averaged results
figure(6);
c = {'r','m','b','c','g'};
for i = 1:length(Tout)-1
    % wave ratio
    subplot(3,2,1)
    plot(xw(:,2),s.Hrmslf(:,i)./max(0.01,sqrt(s.Hrmshf(:,i).^2+s.Hrmslf(:,i).^2)),[c{i},'-'],'LineWidth',2); hold on;
    xlabel('x [m]'); ylabel('Hrms_{lf}/Hrms [-]'); axis([150 max(max(xw)) 0 1.25]);
    legend('Interval A (0.1 hour)','Interval B (0.2 hour)','Interval C (0.7 hour)','Interval D (1.04 hour)','Interval E (3.96 hour)',2);
    % urms ratio
    subplot(3,2,2)
    plot(xw(:,2),s.Urmslf(:,i)./max(0.01,sqrt(s.Urmshf(:,i).^2+s.Urmslf(:,i).^2)),[c{i},'-'],'LineWidth',2); hold on;
    xlabel('x [m]'); ylabel('Urms_{lf}/Urms [-]'); axis([150 max(max(xw)) 0 1.25]);
    % umean 
    subplot(3,2,3)
    plot(xw(:,2),s.Um2(:,i) ,c{i},'LineWidth',2); hold on;
    xlabel('x [m]'); ylabel('U_{mean} [m/s]'); axis([150 max(max(xw)) -1.5 0.1]);
    % cmean
    subplot(3,2,4)
    plot(xw(:,2),s.Cm2(:,i)*par.rhos ,c{i},'LineWidth',2); hold on;
    xlabel('x [m]'); ylabel('C_{mean} [g/l]'); axis([150 max(max(xw)) 0 1.25*par.rhos*(max(max(s.Cm2)))]);
    % cmean
    subplot(3,2,5)
    semilogy(xw(:,2),s.ccceq(:,i) ,c{i},'LineWidth',2); hold on;
    xlabel('x [m]'); ylabel('C/C_{equilibrium} [-]'); axis([150 max(max(xw)) 0.01 10*max(max(s.ccceq(:,i)))]);
    % Sediment transport
    subplot(3,2,6)
    plot(xw(:,2),s.Sdzb(:,i) ,c{i},'LineWidth',2); hold on;
    xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]'); axis([150 max(max(xw)) 1.25*min(min(s.Sdzb)) 1.25*max(max(s.Sdzb))]);
end
for i = 1:length(Tout)-1
    subplot(3,2,1); plot(xhrms,Hrms_lft(i,:)./max(0.01,Hrmst(i,:)),[c{i},'s'],'LineWidth',2); hold on;
end
for i = 1:length(interval);
    subplot(3,2,2); plot(xswf(i),Urms_lft(i)/Urmst(i),[c{interval(i)},'s'],'LineWidth',2); hold on;
    subplot(3,2,3); plot(xswf(i),Umt(i),[c{interval(i)},'s'],'LineWidth',2); hold on;
    subplot(3,2,4); plot(xswf(i),Cmt(i),[c{interval(i)},'s'],'LineWidth',2); hold on;
end
print('intervals.png','-dpng');

figure(7)
t = [0:1:nt-1]*par.tint*par.morfac/3600;
subplot(2,1,1);
plot(t,Vav,'b'); hold on; plot(t,Vflow,'r');
xlabel('t [hours]'); ylabel('V_{avalanch} & V_{flow} [m^3]'); 

save('T01_Sednew.mat','s');

% figure;
% for i = 1:nt
%     plot(xw(:,2),s.zb(:,1),'k--'); hold on;
%     plot(xw(:,2),s.zs(:,i),'b');
%     plot(xw(:,2),s.zb(:,i),'k');
%     plot(xw(:,2),s.ue(:,i),'r');
%     plot(xw(:,2),s.cc(:,i)*par.rhos/10,'g');
%     pause;
%     hold off;
% end

fclose all; close all;