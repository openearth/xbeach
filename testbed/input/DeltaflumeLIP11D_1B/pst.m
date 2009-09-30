function pst(s,m,par,nt,nx,ny,xw,yw);

h = s.zs-s.zb;
s.Cm = mean(s.cc,2);
s.Ceqm = mean(s.ceq,2);
t = (0:1:nt-1)*par.tint*par.morfac/3600;
Tsim = t(end)*3600;
rhos  = 2650;

Lmix = 2*s.R.*par.Tp./max(0.01,(par.rho.*s.c.*s.Hrms));
Lmix = min(Lmix,h+0.5*s.Hrms);
Lmix2 = sqrt(2*s.R.*par.Tp./max(0.01,(par.rho.*s.c)));
Lmix2 = min(Lmix2,h+0.5*s.Hrms);

% transform Su to Sz
Stemp = zeros(size(s.Su));
Stemp(1,:) = s.Su(1,:);
Stemp(2:end,:) = 0.5*(s.Su(1:end-1,:)+s.Su(2:end,:));
dz = s.zb(:,end)-s.zb(:,1);
s.Sdzbtot = (1-par.np)*flipud(cumsum(flipud(dz)))*(xw(2,2)-xw(1,2))/Tsim;
s.Sdzavtot = (1-par.np)*flipud(cumsum(flipud(s.dzav(:,end))))*(xw(2,2)-xw(1,2))/Tsim;
s.Scutot = mean(Stemp,2);

% sediment transports from undertow:
s.ut  = s.ue-s.u;
Slongwaves = s.u.*s.cc.*h;
Sundertow = s.ut.*s.cc.*h;
Swaveassym = s.ua.*s.cc.*h;

hours = [1 2 3 4 5 6 7 8 9 12 18];
hours2 = [0 1 2 3 4 5 6 7 8 9 12 18];
indh = []; indh(1) = 1;
indh(2:length(hours2)) = round(hours*3600/par.morfac/par.tint);

figure(9); set(gca,'FontSize',16);
plot(xw,s.Cm*par.rhos,'k-','LineWidth',1.5); hold on; plot(m.xm2,m.cm,'ks','LineWidth',1.5); 
plot(xw,s.Ceqm*rhos,'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('c [g/l] & c_{eq} [g/l]'); axis([0 max(max(xw)) 0 6]);
print sed_conc_x.png -dpng

figure(10); set(gca,'FontSize',16);
plot(xw(:,2),mean((s.DR/par.rho).^(2/3).*par.Tp./s.Tbore,2),'k--','LineWidth',1.5); hold on;
plot(xw(:,2),mean(s.kb,2),'k','LineWidth',1.5);
xlabel('x [m]'); ylabel('k [m^2/s^2] & k_b [m^2/s^2]'); axis([0 max(max(xw)) 0 1.2]);
print turbulence_x.png -dpng

figure(11); set(gca,'FontSize',16);
plot(xw(:,2),mean(Lmix2./h,2),'k','LineWidth',1.5); hold on;
plot(xw(:,2),mean(min(1,(exp(h./Lmix2)-1).^-1),2),'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('L_{mix}/h [-] & (exp(h/L_{mix})-1)^-1 [-]'); axis([0 max(max(xw)) 0 2]);
print mixing_length_x.png -dpng

figure(12); set(gca,'FontSize',16);
plot(xw(:,2),mean(s.Beta,2),'k','LineWidth',1.5); hold on; 
plot(xw(:,2),mean(tan(asin(s.Beta)),2),'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('\beta [-] & (\partial\eta/\partialx)^-1 [-]'); axis([0 max(max(xw)) 0 1.5]);
print surfaceslope_x.png -dpng

figure(13); set(gca,'FontSize',16);
plot(xw(:,2),mean(s.Sdzavtot,2),'k-','LineWidth',1); hold on;
plot(xw(:,2),mean(s.Sdzbtot,2),'k-','LineWidth',1.5); hold on;
plot(xw(:,2),mean(Slongwaves,2),'k--','LineWidth',1.0);
plot(xw(:,2),mean(Sundertow,2),'k-.','LineWidth',1.0);
plot(xw(:,2),mean(Swaveassym,2),'k:','LineWidth',1.0);
% plot(m.x{end},m.Sdzb{end},'k--','LineWidth',1.5);
% plot(xw(:,2),mean(Swaveassym,2)+mean(Sundertow,2),'r','LineWidth',1.0);
xlabel('x [m]'); ylabel('S [m^3/m/s]'); axis([100 max(max(xw)) -0.00005 0.00002]);
print sediment_transports_x.png -dpng

% transport intervals
figure(14);
for i = 1:length(indh)-1
    subplot(3,4,i);
    plot(xw(:,2),mean(Stemp(:,indh(i):indh(i+1)),2),'k-','LineWidth',2.5); hold on;
    plot(xw(:,2),mean(Slongwaves(:,indh(i):indh(i+1)),2),'k--','LineWidth',1.5); plot(xw(:,2),mean(Sundertow(:,indh(i):indh(i+1)),2),'k-.','LineWidth',1.5); plot(xw(:,2),mean(Swaveassym(:,indh(i):indh(i+1)),2),'k:','LineWidth',1.5); 
    %plot(xw(:,2),mean(Slongwaves(:,indh(i):indh(i+1))+Sundertow(:,indh(i):indh(i+1))+Swaveassym(:,indh(i):indh(i+1)),2),'k--','LineWidth',2.5); 
    xlabel('x [m]'); ylabel('S [m^3/m/s]'); axis([100 220 -0.00015 0.00005]);
end
print decomposed_sediment_transport_intervals.png -dpng


% figure(4)
% subplot(211);
% plot(xw(:,2),s.Sdzbtot,'k--','LineWidth',1.5); hold on; plot(m.x{1},m.Sdzb{end},'k-','LineWidth',1.5);
% plot(xw(:,2),s.Sdzavtot,'r--','LineWidth',1.5); plot(xw(:,2),s.Scutot,'r-','LineWidth',1.5); 
% xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]'); axis([0 max(max(xw)) -3E-5 1E-5]);
% subplot(212);
% plot(xw(:,2),s.Scutot,'k-','LineWidth',1.5); hold on;
% plot(xw(:,2),mean(Slongwaves,2),'g-','LineWidth',1.5); 
% plot(xw(:,2),mean(Sundertow,2),'r-','LineWidth',1.5); 
% plot(xw(:,2),mean(Swaveassym,2),'c-','LineWidth',1.5); 
% plot(xw(:,2),mean(Slongwaves+Sundertow,2),'k--','LineWidth',1.5); 
% plot(xw(:,2),mean(Slongwaves+Sundertow-Swaveassym,2),'k:','LineWidth',1.5); 
% xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]'); axis([0 max(max(xw)) -3E-5 1E-5]);


