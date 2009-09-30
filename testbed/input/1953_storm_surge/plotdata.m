clear s

morfac = 1;
tint   = 10;
%% simulation results
nam = [{'zb'};{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'cc'};{'Su'};{'ua'};];

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

% % setup
% s.setup = mean(s.zs-max(0,s.zb),2);
% 
% % long wave height
% h = s.zs-s.zb;
% windowSize = 60; % ten seconds window 
% fs = fspecial('average',windowSize);
% hm = [];
% for i=1:length(xw)
%     hm(i,:) = imfilter(h(i,:),fs,'symmetric');
% end
% zs = h-hm;
% 
% s.Hrms_lfm = sqrt(8)*std(zs')';
% s.Hrms_lfm2 = sqrt(8)*std(h')';
% s.Hrms_hfm = sqrt(mean(s.Hrms.^2,2));
% % s.Hrms_lfm = sqrt(8)*std(detrend((s.zs-s.zb)'))';
% s.Hrms_m = sqrt(s.Hrms_hfm.^2 + s.Hrms_lfm.^2);
% s.Urms_hfm = sqrt(mean(s.urms.^2,2));
% s.Urms_lfm = std(s.u')';
% s.Urms_m = sqrt(s.Urms_hfm.^2 + s.Urms_lfm.^2);
% s.Uon_m = 1/sqrt(2)*mean(s.uon,2);
% s.Uoff_m = 1/sqrt(2)*mean(s.uoff,2);
% s.Um = mean(s.ue,2);
% s.Ua = mean(s.ua,2);
% s.Skm = mean(s.Sk,2);
% s.Asm = mean(s.As,2);
% s.Cm = mean(s.cc,2);
% s.Ceqm = mean(s.ceq,2);
% dz = s.zb(:,end)-s.zb(:,1);
% s.usR = mean(2*s.R./(par.rho.*max(h,0.1).*s.c),2);
% s.usW = mean(0.125*par.g*s.Hrms.^2./(max(h,0.1).*s.c),2);

%% figures;
% load post storm profile
% addpath E:\jsmthieldevries\delft3d\toolbox\
% wlsettings('E:\jsmthieldevries\delft3d\toolbox\');
% bath = tekal('open','E:\jsmthieldevries\data\duinafslag\M1263-T4\zm_post.tek');
% xz2 = bath.Field.Data;
% xz2(:,1) = xz2(:,1)*4.56;
% xz2(:,2) = xz2(:,2)*3.3;

% to NAP..
zsmax = 3.9;
zsmin = 0.93;
dz = 13.86-zsmax;

figure; 
subplot(211);
plot(xw,s.zb(:,1)-dz,'k--','LineWidth',1.5); hold on;
plot(xw,s.zb(:,end)-dz,'k','LineWidth',1.5);
plot(xw,repmat(zsmax,1,nx),'k:');
plot(xw,repmat(zsmin,1,nx),'k:');
axis([700 1100 -4 12]);
xlabel('x [m]'); ylabel('z_b [m]');
print profile_evolution.png -dpng

figure; 
subplot(211);
plot(xw,mean(s.Hrms,2)+max(s.zs(1,:))-dz,'r','LineWidth',1.5);
plot(xw,mean(s.ue,2)+max(s.zs(1,:))-dz,'r--','LineWidth',1.5);
% plot(xw,mean(s.ua,2)+max(s.zs(1,:))-dz,'r--','LineWidth',1.5);
plot(xw,mean(s.cc,2)*2650+max(s.zs(1,:))-dz,'r-.','LineWidth',1.5);

% erosion volume above max storm surge level
A0 = sum(max(s.zb(:,1)-dz,zsmax))*4.56;
AE = sum(max(s.zb(:,end)-dz,zsmax))*4.56;
Atot = A0-AE;

% total erosion volume
Atot2 = sum(max(s.zb(:,1)-s.zb(:,end),0))*4.56;

% title(['total erosion above max surge level = ',num2str(AA),' [m^3/m] ']);

t = [0:1:nt-1]/3600*tint*morfac;
t1 = [0.000	 9150.38 21581.09 26164.13 33227.03 42502.98 52312.57 54274.48 65402.48 78476.70 85021.65 91550.91 98095.87 104640.83 111170.09]/3600;
t2 = cumsum([0.000 18316.461 0 22224.601 0 19619.174 0 24861.418 0 26148.435 ])/3600;
tp = 1./[0.1135 0.1135 0.1135 0.1135 0.10944 0.10944 0.0983 0.0983 0.0983 0.0983];
hrms =[4.9516 4.9516 5.4463 5.4463 6.039 6.039 6.4357 6.4357 6.1043 6.1043]/sqrt(2);
zs = [11.385 12.210 12.144 12.474 12.210 12.870 13.860 13.860 12.375 11.550 11.550 11.880 12.375 11.979 10.890];
figure; 
subplot(311); plot(t,s.zs(1,:)-dz,'k'); hold on; plot(t1,zs-dz,'color',[0.6 0.6 0.6],'LineWidth',3);
xlabel('t [hours]'); ylabel('\eta [m]'); axis([0 t(end) 0 5]); 
subplot(312); plot(t,s.Hrms(1,:),'k'); hold on; plot(t2,hrms,'color',[0.6 0.6 0.6],'LineWidth',3);
xlabel('t [hours]'); ylabel('H_{rms} [m]'); axis([0 t(end) 0 12]); 
subplot(313); plot(t2,tp,'color',[0.6 0.6 0.6],'LineWidth',3);
xlabel('t [hours]'); ylabel('T_p [s]'); axis([0 t(end) 8 12]); 
print hydrodynamics.png -dpng


% % movieplot
% figure;
% for i = nt-1000:nt
%     plot(xw,s.zb(:,1),'k--'); hold on;
%     plot(xw,s.zs(:,i),'b');
%     plot(xw,s.zb(:,i),'k');
%     plot(xw,s.ue(:,i)+max(s.zs(1,:)),'r');
%     plot(xw,s.cc(:,i)*2650+max(s.zs(1,:)),'g');
%     axis([0 max(xw(:,2)) 0 s.zb(end,1)+2.5]);
%     pause(0.1); hold off;
% end

% fclose all; close all;




