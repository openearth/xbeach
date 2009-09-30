function pst(s,m,par,nt,nx,ny,xw,yw);

x = 0:1:230;

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

% corrlation long waves and (near bed) turbulence
for i = 1:nx
     [R,P] = corrcoef(s.DR(i,:)',detrend(s.u(i,:))');
     rho(i) = R(1,2); pval(i) = P(1,2);
     [R,P] = corrcoef(s.cc(i,:)',detrend(s.u(i,:))');
     rho2(i) = R(1,2); pval2(i) = P(1,2);
     [R,P] = corrcoef(s.kb(i,:)',detrend(s.u(i,:))');
     rho3(i) = R(1,2); pval3(i) = P(1,2);
end

% what are simulated sed.conc. at times of mobile frame measurements 
MF = 1;
hours = [0 0.1 0.3 1.0 2.04 6.0];
sec = hours*3600;
if sum(isnan(m.tMF))==0
for i = 1:length(m.tMF);    
    seccor(i) = sec(m.interval(i));
end
tMF = m.tMF+seccor;
% now link with simulations...
windowSize = 90; % ten seconds window 
fs = fspecial('average',windowSize);
cma = [];
for i=1:length(xw)
    cma(i,:) = imfilter(s.cc(i,:),fs,'symmetric');
end
ind1 = round(tMF/par.morfac/par.tint);
for i = 1:length(m.xswf)
    ind2(i) = find(m.xswf(i)==xw(:,2));
end
else
    MF = 0;
end
%% sediment concentrations
figure(11); set(gca,'FontSize',16);
plot(xw,s.Cm*par.rhos,'k-','LineWidth',1.5); hold on; plot(m.xum,m.Cm,'ks','LineWidth',1.5); 
% plot(xw,s.Ceqm*rhos,'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('c [g/l]'); axis([0 max(max(xw)) 0 2*rhos*max(max(s.Cm))]);
print sed_conc_x.png -dpng;

figure(12); set(gca,'FontSize',16);
plot(xw,s.Cm*par.rhos,'k-','LineWidth',1.5); hold on; plot(m.xum,m.Cm,'ks','LineWidth',1.5); 
% plot(xw,s.Ceqm*rhos,'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('c [g/l] & c_{eq} [g/l]'); axis([160 max(max(xw)) 0 25]);
print sed_conc_x_zoom.png -dpng;

xt = [170 175 180 185 190 195 200 205];
for i = 1:length(xt)
    ind(i) = find(xw(:,2)==xt(i));
end
fact = nanmean(par.rhos*mean(s.cc(ind,:),2)./m.Cm')

if MF==1
% scatter plot simulated and measured sediment concentrations
figure(13); set(gca,'FontSize',16);
for i = 1:length(ind1)
    loglog(m.Cmt(i),rhos*cma(ind2(i),ind1(i)),'ks','LineWidth',1.5); hold on;
end
plot([0.1 50],[0.1 50],'k','LineWidth',1.5);
plot([0.1 50],0.5*[0.1 50],'k--','LineWidth',1.5);
plot([0.1 50],2*[0.1 50],'k--','LineWidth',1.5);
xlabel('c_{measured} [g/l]'); ylabel('c_{simulated} [g/l]'); axis image
axis([0.1 50 0.1 50]); grid on;
print sed_conc_scatter.png -dpng;
end

% wave averaged, bore averaged and near bed turbulence
figure(14); set(gca,'FontSize',16);
plot(xw(:,2),mean((s.DR/par.rho).^(2/3),2),'k:','LineWidth',1.5); hold on;
plot(xw(:,2),mean((s.DR/par.rho).^(2/3).*par.Tp./s.Tbore,2),'k--','LineWidth',1.5); hold on;
plot(xw(:,2),mean(s.kb,2),'k','LineWidth',1.5);
xlabel('x [m]'); ylabel('k [m^2/s^2]'); axis([41 max(max(xw)) 0 1.2*max(max(mean((s.DR/par.rho).^(2/3).*par.Tp./s.Tbore,2)))]);
print turbulence.png -dpng;

% mixing length
figure(15); set(gca,'FontSize',16);
plot(xw(:,2),mean(Lmix./h,2),'k','LineWidth',1.5); hold on; plot(xw(:,2),mean(Lmix2./h,2),'k--','LineWidth',1.5); plot(xw(:,2),mean(0.9*s.Hrms./h,2),'k:','LineWidth',1.5);
plot(xw(:,2),mean(min(1,(exp(h./Lmix)-1).^-1),2),'r-','LineWidth',1.5);
xlabel('x [m]'); ylabel('L_{mix}/h [-] & (exp(h/L_{mix})-1)^-1 [-]'); axis([41 max(max(xw)) 0 1.5]);
print mixinglength.png -dpng

% surface slope
figure(16); set(gca,'FontSize',16);
plot(xw(:,2),mean(s.Beta,2),'k','LineWidth',1.5); hold on; 
plot(xw(:,2),mean(tan(asin(s.Beta)),2),'k--','LineWidth',1.5);
xlabel('x [m]'); ylabel('\beta & \partial\eta/\partialx [-]'); axis([41 max(max(xw)) 0 1.2*max(mean(tan(asin(s.Beta)),2))]);
print surfaceslope.png -dpng

% correlations:
% 1) roller dissipation and long wave flows
% 2) sediment concentration and long wave flows
% 3) near bed turbulence and long wave flows
figure(17); set(gca,'FontSize',16);
plot(xw,rho,'k-','LineWidth',1.5); hold on; plot(xw,rho2,'k--','LineWidth',1.5); hold on; plot(xw,rho3,'k:','LineWidth',1.5); hold on;
grid on; xlabel('x [m]'); ylabel('\rho_{(D_R,u)} & \rho_{(c,u)} & \rho_{(k_b,u)} [-]');
axis([41 max(max(xw)) -1 1]);

% sediment transports
figure(18); set(gca,'FontSize',16);
plot(xw(:,2),s.Sdzbtot,'k-','LineWidth',2.5); hold on; plot(x,m.Sdzb,'k--','LineWidth',2.5);
plot(xw(:,2),s.Sdzavtot,'k-.','LineWidth',1.5); plot(xw(:,2),s.Scutot,'k:','LineWidth',1.5); 
xlabel('x [m]'); ylabel('S [m^3/m/s]'); axis([41 max(max(xw)) 1.2*min(min(s.Sdzbtot),min(m.Sdzb)) 2*max(max(s.Sdzbtot),max(m.Sdzb))]);
print sediment_transport.png -dpng

% decompose flow-based sediment transports
figure(19); set(gca,'FontSize',16);
plot(xw(:,2),s.Scutot,'k-','LineWidth',2.5); hold on;
plot(xw(:,2),mean(Slongwaves,2),'k--','LineWidth',1.5); plot(xw(:,2),mean(Sundertow,2),'k-.','LineWidth',1.5); plot(xw(:,2),mean(Swaveassym,2),'k:','LineWidth',1.5); 
%plot(xw(:,2),mean(Slongwaves+Sundertow+Swaveassym,2),'k--','LineWidth',2.5); 
xlabel('x [m]'); ylabel('S [m^3/m/s]'); axis([41 max(max(xw)) 1.2*min(min(s.Sdzbtot),min(m.Sdzb)) 1.2*max(max(mean(Swaveassym,2)),max(m.Sdzb))]);
print decomposed_sediment_transport.png -dpng


