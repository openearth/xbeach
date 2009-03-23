function pu(s,m,par,nt,nx,ny,xw,yw);

h = s.zs-s.zb;
s.Urms_hfm = sqrt(mean(s.urms.^2,2));
s.Urms_lfm = std(s.u')';
s.Urms_m = sqrt(s.Urms_hfm.^2 + s.Urms_lfm.^2);
s.Uon_m = 1/sqrt(2)*mean(s.uon,2);
s.Uoff_m = 1/sqrt(2)*mean(s.uoff,2);
s.Um = mean(s.ue,2);
s.Ua = mean(s.ua,2);

ee = 0.125*par.rho*par.g*s.H.^2;
ussw = ee./max(s.c,sqrt(par.g*par.hmin))./par.rho./max(s.hh,par.hmin);
usr = 2*s.R./max(s.c,sqrt(par.g*par.hmin))./par.rho./max(s.hh,par.hmin);

fac = 0.15;

% return flow
figure(8); set(gca,'FontSize',16);
plot(xw,s.Um,'k-','LineWidth',1.5); hold on; plot(m.xum,m.Um,'kv','LineWidth',1.5); % plot(m.xm2,m.um2,'rv','LineWidth',1.5);
plot(xw,s.Ua,'k:','LineWidth',1.5); plot(m.xhrms,fac*m.ua,'k^','LineWidth',1.5); % plot(m.xurms,1/fac*m.uaMF,'k.','LineWidth',1.5);
% plot(xw,s.Um+s.Ua,'k--','LineWidth',1.5); 
plot(xw,mean(s.u,2),'k-.','LineWidth',1.5); 
xlabel('x [m]'); ylabel('U_{m} & U_{A} [m/s]'); axis([min(min(xw)) max(max(xw)) 1.5*min(min(s.Um)) 1.5*max(max(s.Ua))]);
print mean_flow.png -dpng

figure(9); set(gca,'FontSize',16);
plot(xw,s.Urms_m,'k-','LineWidth',1.5); hold on; plot(m.xurms,m.Urms,'ks','LineWidth',1.5);  
plot(xw,s.Urms_hfm,'k--','LineWidth',1.5); hold on; plot(m.xurms,m.Urms_hf,'k^','LineWidth',1.5); 
plot(xw,s.Urms_lfm,'k-.','LineWidth',1.5); hold on; plot(m.xurms,m.Urms_lf,'kv','LineWidth',1.5); 
xlabel('x [m]'); ylabel('U_{rms} [m/s]'); axis([min(min(xw)) max(max(xw)) 0 1.5*max(s.Urms_m)]);
print orbital_flow.png -dpng

figure(10); set(gca,'FontSize',16);
plot(xw(:,2),mean(s.ue,2),'k','LineWidth',1.5); hold on;
plot(xw(:,2),mean(s.u,2),'k--','LineWidth',1.5);
plot(xw(:,2),mean(-ussw,2),'k-.','LineWidth',1.5);
plot(xw(:,2),mean(-usr,2),'k:','LineWidth',1.5);
% plot(xw(:,2),mean(-usr-ussw+s.u,2),'r','LineWidth',1.5);
plot(m.xum,m.Um,'ks','LineWidth',1.5);
xlabel('x [m]'); ylabel('U_{m} [m]'); axis([160 max(max(xw)) -0.5 0]);
print mean_flow_decomp.png -dpng

% figure(11);
% x = [170:5:205];
% for i = 1:length(x)
%     ind = find(x(i)==m.xswf);
%     subplot(2,4,i);
%     for j = 1:length(ind)
%         semilogx(m.cm{ind(j)},m.zcm{ind(j)},'ks'); hold on;
%     end
%     ind2 = find(x(i)==xw(:,2));
%     htemp = (1+mean(s.sig(i,:,:),3))*mean(s.hh(ind2,:),2);
%     semilogx(par.rhos*mean(s.cuq3d(i,:,:),3),htemp,'k'); hold on;
% end
% 
% figure(12);
% for i = 1:length(x)
%     ind = find(x(i)==m.xswf);
%     subplot(2,4,i);
%     for j = 1:length(ind)
%         plot(m.um{ind(j)},m.zum{ind(j)},'ks','MarkerSize',2); hold on;
%     end
%     % plot(m.um{ind},m.zum{ind},'ks','MarkerSize',2);
%     ind2 = find(xw(:,2)==x(i));
%     htemp = (1+mean(s.sig(i,:,:),3))*mean(s.hh(ind2,:),2);
%     plot(mean(s.veloc(i,:,1,:),4),htemp,'k'); hold on;
%     title(['x = ',num2str(x(i)),' [m]']);
% end
% 
% figure(12);
% for i = 1:length(x)
%     ind = find(x(i)==m.xswf);
%     subplot(2,4,i);
%     % plot(m.um{ind},m.zum{ind},'ks','MarkerSize',2);
%     ind = find(xw(:,2)==x(i));
%     htemp = (1+mean(s.sig(i,:,:),3))*mean(s.hh(ind,:),2);
%     plot(mean(s.veloc(i,:,3,:),4),htemp,'k'); hold on;
%     title(['x = ',num2str(x(i)),' [m]']);
% end

