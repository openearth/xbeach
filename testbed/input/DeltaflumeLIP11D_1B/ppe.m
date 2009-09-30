function ppe(s,m,par,nt,nx,ny,xw,yw);
figure(15); set(gca,'FontSize',16);
subplot(211);
plot(m.x{1},m.z{1},'k'); hold on;
plot(xw(:,2),s.zb(:,end),'k--','LineWidth',1.5); 
plot(m.x{end},m.z{end},'k','LineWidth',1.5);
plot([0 xw(end,2)],[0 0],'k:');
axis([0 max(max(xw)) -4.5 1.5]); xlabel('x [m]'); ylabel('z_{b} [m]'); 
print profile_evolution.png -dpng

figure(16); set(gca,'FontSize',16);
subplot(211);
plot(m.x{1},m.z{1},'k');hold on;
plot(xw(:,2),s.zb(:,end),'k--','LineWidth',1.5); % +0.055
plot(m.x{end},m.z{end},'k','LineWidth',1.5);
plot([0 xw(end,2)],[0 0],'k:')
axis([100 max(max(xw)) -2.0 1.5]); xlabel('x [m]'); ylabel('z_{b} [m]'); 
print profile_evolution_detailed.png -dpng

