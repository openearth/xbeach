function pu(s,m,par,nt,nx,ny,xw,yw);

h = s.zs-s.zb;
s.Urms_hfm = sqrt(mean(s.urms.^2,2));
s.Urms_lfm = std(s.u')';
s.Urms_m = sqrt(s.Urms_hfm.^2 + s.Urms_lfm.^2);
s.Uon_m = 1/sqrt(2)*mean(s.uon,2);
s.Uoff_m = 1/sqrt(2)*mean(s.uoff,2);
s.Um = mean(s.ue,2);
s.Ua = mean(s.ua,2);
s.Skm = mean(s.Sk,2);
s.Asm = mean(s.As,2);
ee = 0.125*par.rho*par.g*s.Hrms.^2;
ussw = ee./max(s.c,sqrt(par.g*par.hmin))./par.rho./max(s.hh,0.5*par.hmin);
usr = 2*s.R./max(s.c,sqrt(par.g*par.hmin))./par.rho./max(s.hh,0.5*par.hmin);

hours = [1 2 3 4 5 6 7 8 9 12 18];
hours2 = [0 1 2 3 4 5 6 7 8 9 12 18];
indh = []; indh(1) = 1;
indh(2:length(hours2)) = round(hours*3600/par.morfac/par.tint);

figure(4); set(gca,'FontSize',16);
plot(xw,s.Um,'k-','LineWidth',1.5); hold on; plot(m.xm2,m.um,'kv','LineWidth',1.5);
plot(xw,s.Ua,'k:','LineWidth',1.5); plot(m.xps,0.1*mean(m.ua),'k^','LineWidth',1.5);
plot(xw,mean(s.u,2),'k--','LineWidth',1.5); 
xlabel('x [m]'); ylabel('U_{m} & U_{A} [m/s]'); axis([min(min(xw)) max(max(xw)) -0.3 0.2]);
print('mean_flow.png','-dpng');

figure(5); set(gca,'FontSize',16);
plot(xw,s.Urms_m,'k-','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.urms1.^2)),'ks','LineWidth',1.5);  
plot(xw,s.Urms_hfm,'k--','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.urms_hf1.^2)),'k^','LineWidth',1.5); 
plot(xw,s.Urms_lfm,'k-.','LineWidth',1.5); hold on; plot(m.xps,sqrt(mean(m.urms_lf1.^2)),'kv','LineWidth',1.5); 
xlabel('x [m]'); ylabel('U_{rms} [m]'); axis([min(min(xw)) max(max(xw)) 0 1.0]);
print('orbital_flows.png','-dpng');

figure(6); set(gca,'FontSize',16);
Beta = atan(s.Asm./s.Skm); B = sqrt(s.Asm.^2+s.Skm.^2);
BetaFW = mean(atan(m.As1./m.Sk1)); BFW = mean(sqrt(m.As1.^2+m.Sk1.^2));
plot(xw(:,2),s.Skm,'k-','LineWidth',1.5); hold on; plot(m.xps,m.Sk1,'k^','LineWidth',1.5); 
plot(xw(:,2),s.Asm,'k--','LineWidth',1.5); plot(m.xps,m.As1,'kv','LineWidth',1.5); 
xlabel('x [m]'); ylabel('S_K & A_S [-]'); axis([min(min(xw)) max(max(xw)) -2 2]);
print('Skewness&Asummtery.png','-dpng');

figure(7); set(gca,'FontSize',16);
plot(xw(:,2),Beta,'k--','LineWidth',1.5); hold on; plot(m.xps,BetaFW,'k^','LineWidth',1.5);
plot(xw(:,2),B,'k-','LineWidth',1.5); plot(m.xps,BFW,'kv','LineWidth',1.5);
xlabel('x [m]'); ylabel('\beta [-] & B [-]'); axis([min(min(xw)) max(max(xw)) -2 2]);
print('non_linearity&phase.png','-dpng'); 

% flows decomposed for intervals
figure(8); set(gca,'FontSize',16);
for i = 1:length(indh)-1
    subplot(3,4,i);
    plot(xw(:,2),mean(s.ue(:,indh(i):indh(i+1)),2),'k','LineWidth',1.5); hold on;
    plot(xw(:,2),mean(s.u(:,indh(i):indh(i+1)),2),'k--','LineWidth',1.5);
    plot(xw(:,2),mean(-ussw(:,indh(i):indh(i+1)),2),'k-.','LineWidth',1.5);
    plot(xw(:,2),mean(-usr(:,indh(i):indh(i+1)),2),'k:','LineWidth',1.5);
    xlabel('x [m]'); ylabel('U_{m} [m]'); axis([160 max(max(xw)) -0.5 0]);
end
print mean_flow_decomp_intervals.png -dpng

