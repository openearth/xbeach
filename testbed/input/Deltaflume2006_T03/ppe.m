function ppe(s,m,par,nt,nx,ny,xw,yw);

x= 0:1:230; dx = 1;
hours = [0.1 0.3 1.0 2.04 6.0];
hours2 = [0.0 0.1 0.3 1.0 2.04 6.0];
indh = []; indh(1) = 1;
indh(2:6) = round(hours*3600/par.morfac/par.tint);

% compute erosion volumes above maximum still water level
Am = zeros(1,6); As = Am;
for i = 2:6;
    Am(i) = sum(max(0,m.z(1,:))-max(0,m.z(i,:)))*dx;
    As(i) = sum(max(0,s.zb(:,1))-max(0,s.zb(:,indh(i))))*dx;
end
Am1 = [0 1.15 2.48 5.31 7.10 9.85];
Am2 = [0 1.04 2.17 5.24 7.28 10.0];

figure(20); subplot(211); set(gca,'FontSize',10);
for i = 1:6
    plot(x,m.z(i,:),'k'); hold on;
    plot(xw(:,2),s.zb(:,indh(i)),'k--','LineWidth',1.5);
end
plot([0 xw(end,2)],[0 0],'k')
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([0 max(max(xw)) -4.5 2]);
print profiles.png -dpng

figure(21); subplot(211); set(gca,'FontSize',10);
for i = 1:6
    plot(x,m.z(i,:),'k'); hold on;
    plot(xw(:,2),s.zb(:,indh(i)),'k--','LineWidth',1.5);
end
plot([0 xw(end,2)],[0 0],'k');
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([160 max(max(xw)) -2.0 2]);
print detailed_profiles.png -dpng

figure(22); set(gca,'FontSize',16);
plot(hours2,Am(1,:),'k-s','LineWidth',1.5,'MarkerFaceColor','k'); hold on;
plot(hours2,As(1,:),'k--s','LineWidth',1.5,'MarkerSize',12);
plot(hours2,Am1,'k-o','LineWidth',1.5,'MarkerFaceColor','k');
plot(hours2,Am2,'k-o','LineWidth',1.5,'MarkerFaceColor','k');

xlabel('t [hours]'); ylabel('A [m^3/m]');
print t_dune_erosion.png -dpng


