function plotdata()
clear s;close all
[runid,testid,datadir]=testinfo

%% simulation results
nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'ccg'};{'Sug'};{'ua'};];

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

%% figures;
% load post storm profile

figure(2); set(gca,'FontSize',16);
plot(xw,s.zb(:,1),'k--','LineWidth',1.5); hold on;
plot(xw,s.zb(:,end),'r','LineWidth',1.5);
plot(xw,mean(s.H,2)+max(s.zs(1,:)),'b','LineWidth',1.5);
plot(xw,mean(s.ue,2)+max(s.zs(1,:)),'r','LineWidth',1.5);
plot(xw,mean(s.ua,2)+max(s.zs(1,:)),'r--','LineWidth',1.5);
plot(xw,mean(s.ccg,2)*2650+max(s.zs(1,:)),'g','LineWidth',1.5);
plot(xw,repmat(max(s.zs(1,:)),1,nx),'k:');
axis([0 max(xw(:,2)) -20 s.zb(end,1)+2.5]);
xlabel('x [m]'); ylabel('z_b [m]');
A0 = sum(max(0,s.zb(:,1)-max(s.zs(1,:))))*5;
AE = sum(max(0,s.zb(:,end)-max(s.zs(1,:))))*5;
AA = A0-AE;
title(['total erosion above max surge level = ',num2str(AA),' [m^3/m ']);
print('profile_evolution.png','-dpng');

movie = 1;
if movie == 1

figure;
for i = nt
    if mod(i,10)==0
    plot(xw,s.zb(:,1),'k--'); hold on;
    plot(xw,s.zs(:,i),'b');
    plot(xw,s.zb(:,i),'k');
    plot(xw,s.ue(:,i)+max(s.zs(1,:)),'r');
    plot(xw,s.ccg(:,i)*2650+max(s.zs(1,:)),'g');
    pause(0.1); hold off;
    end
end

end
pname = ['..\..\report\',testid '_' runid '_fig1' '.jpg'];
eval(['print -djpeg ' pname]);



