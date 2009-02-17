function ppe(s,m,par,nt,nx,ny,xw,yw);

x= 0:1:230;
hours = [0.1 0.3 1.0 2.04 6.0];
indh = []; indh(1) = 1;
indh(2:6) = round(hours*3600/par.morfac/par.tint);

figure(19); subplot(211); set(gca,'FontSize',10);
for i = 1:6
    plot(x,m.z(i,:),'k'); hold on;
    plot(xw(:,2),s.zb(:,indh(i)),'k--','LineWidth',1.5);
end
plot([0 xw(end,2)],[0 0],'k')
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([0 max(max(xw)) -4.5 2]);
print profiles.png -dpng

figure(20); subplot(211); set(gca,'FontSize',10);
for i = 1:6
    plot(x,m.z(i,:),'k'); hold on;
    plot(xw(:,2),s.zb(:,indh(i)),'k--','LineWidth',1.5);
end
plot([0 xw(end,2)],[0 0],'k');
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([160 max(max(xw)) -2.0 2]);
print detailed_profiles.png -dpng

