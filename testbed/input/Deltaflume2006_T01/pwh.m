function pwh(s,m,par,nt,nx,ny,xw,yw);

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
s.Hrms_hfm = sqrt(mean(s.Hrms.^2,2));
s.Hrms_m = sqrt(s.Hrms_hfm.^2 + s.Hrms_lfm.^2);

% Analysis long waves
w = 2*pi/par.Tp;
hh = max(mean(h,2),0.01); % max(h,0.01); 
k = disper(w,hh,par.g);
cg = w./k*0.5.*(1+2*k.*hh./sinh(2*k.*hh));
cg = sqrt(par.g*hh); % cg = sqrt(par.g*hh);
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
    Rt = corrcoef(detrend(s.zs(i,:))',s.Hrms(i,:).^2');
    s.rho(i) = Rt(1,2);
end

% analyse short wave shape 
s.Skm = mean(s.Sk,2);
s.Asm = mean(s.As,2);
s.usR = mean(2*s.R./(par.rho.*max(h,0.1).*s.c),2);
s.usW = mean(0.125*par.g*s.Hrms.^2./(max(h,0.1).*s.c),2);
Beta = atan(s.Asm./s.Skm); B = sqrt(s.Asm.^2+s.Skm.^2);
BetaFW =atan(m.Asu./m.Sku); BFW = sqrt(m.Asu.^2+m.Sku.^2);
BetaMF =atan(m.AsuMF./m.SkuMF); BMF = sqrt(m.AsuMF.^2+m.SkuMF.^2);

error = [];

% wave height transformation
figure(1); set(gca,'FontSize',16);
plot(xw,s.Hrms_m,'k-','LineWidth',1.5); hold on; plot(m.xhrms,m.Hrms,'ks','LineWidth',1.5);
plot(xw,s.Hrms_hfm,'k--','LineWidth',1.5); hold on; plot(m.xhrms,m.Hrms_hf,'k^','LineWidth',1.5);
plot(xw,s.Hrms_lfm,'k-.','LineWidth',1.5); hold on; plot(m.xhrms,m.Hrms_lf,'kv','LineWidth',1.5);
plot(xw,s.setup,'k:','LineWidth',1.5); hold on; plot(m.xhrms,m.setup,'ko','LineWidth',1.5);
xlabel('x [m]'); ylabel('H_{rms} & \eta_{m} [m]'); axis([min(min(xw)) max(max(xw)) -0.1 1.20]);
print wave_height_transformation.png -dpng

% error(1) = 

% correlation wave groups and long waves
figure(2); set(gca,'FontSize',16);
plot(xw,s.rho,'k','LineWidth',1.5); hold on; plot(m.xhrms,m.rho,'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('\rho [-]'); axis([min(min(xw)) max(max(xw)) -1 1]);
print correlation_groups_long.png -dpng

% wave shape evolution: skewness and assymetry 
figure(3); set(gca,'FontSize',16);
plot(xw(:,2),s.Skm,'k-','LineWidth',1.5); hold on; plot(m.xhrms,m.Sku,'k^','LineWidth',1.5); % plot(m.xurms,m.SkuMF,'k.','LineWidth',1.5); 
plot(xw(:,2),s.Asm,'k--','LineWidth',1.5); plot(m.xhrms,m.Asu,'kv','LineWidth',1.5); % plot(m.xurms,m.AsuMF,'k.','LineWidth',1.5); 
xlabel('x [m]'); ylabel('S_K & A_S [-]'); axis([min(min(xw)) max(max(xw)) -2 2]);
print skewness&assyemtry.png -dpng

% wave shape evolution: non-linearity and phase
figure(4); set(gca,'FontSize',16);
plot(xw(:,2),Beta,'k--','LineWidth',1.5); hold on; plot(m.xhrms,BetaFW,'kv','LineWidth',1.5); % plot(m.xurms,BetaMF,'k.','LineWidth',1.5);
plot(xw(:,2),B,'k-','LineWidth',1.5); plot(m.xhrms,BFW,'k^','LineWidth',1.5); % plot(m.xurms,BMF,'k.','LineWidth',1.5);
xlabel('x [m]'); ylabel('\beta [-] & B [-]'); axis([min(min(xw)) max(max(xw)) -2 2]);
print nonlinearity&phase.png -dpng

% incoming and outgoing waves
figure(5); set(gca,'FontSize',16);
plot(xw(:,2),s.Hrms_lfm,'k','LineWidth',1.5); hold on; plot(m.xhrms,m.Hrms_lf,'kv','LineWidth',1.5);
plot(xw(:,2),s.Hrms_lfm_in,'k--','LineWidth',1.5); 
plot(xw(:,2),s.Hrms_lfm_out,'k-.','LineWidth',1.5);
% plot(xw(:,2),sqrt(s.Hrms_lfm_in.^2+s.Hrms_lfm_out.^2),'r--','LineWidth',1.5);
plot(xw(:,2),s.Hrms_lfm_out(1)*hh.^-0.25/hh(1).^-0.25,'k--'); % Greens law
plot(xw(:,2),s.Hrms_lfm_in(1)*hh.^-2.50/hh(1).^-2.50,'k-.');    % Longuet-Higgins and Stewart
legend('Long wave height','measurements','Incoming long wave height','Reflected long wave height','Green: H ~ h^{-1/4}','Longuet-Higgins & Stewart: H ~ h^{-5/2}',2);
xlabel('x [m]'); ylabel('H_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
print incoming&outgoingwaves.png -dpng

% roller energy, wave dissiipation and roller dissipation
figure(6); set(gca,'FontSize',16);
% plot(xw(:,2),mean(0.125*par.rho*par.g*s.Hrms.^2,2),'k','LineWidth',1.5);
plot(xw(:,2),mean(s.R,2),'k--','LineWidth',1.5); hold on;
plot(xw(:,2),mean(s.DR,2),'k-.','LineWidth',1.5);
plot(xw(:,2),mean(s.D,2),'k-','LineWidth',1.5);
xlabel('x [m]'); ylabel('E_R [J/m^2] & D, D_R [J/m^2/s]'); axis([160 max(max(xw)) 0 125]); %1.2*max(mean(s.R,2))
print roller_energy.png -dpng

figure(7); set(gca,'FontSize',16);
plot(xw(:,2),mean(s.RR,2),'k-','LineWidth',1.5);
xlabel('x [m]'); ylabel('R [W/m^2]');
print energy_transfer.png -dpng
