clear s m

%% coefficients and settings
fid = fopen('..\data\datanam.txt');
Tnam = fgetl(fid);
nvar = str2num(fgetl(fid));
for i = 1:nvar-1
    line = fgetl(fid);
    [nam,remain] = strtok(line,' ');
    [dummy,remain] = strtok(remain,'=');
    par.(nam) = str2num(remain(2:end));
end
line = fgetl(fid);
[nam,remain] = strtok(line,' ');
[dummy,remain] = strtok(remain,'=');
par.(nam) = remain(3:end);
fclose(fid);
%% measurements
% load(['..\data\',Tnam]);
% m.Hrms = Hrms; m.Hrms_hf = Hrms_hf; m.Hrms_lf = Hrms_lf; m.Hrms_in = Hrms_in; m.Hrms_out = Hrms_out;
% m.setup = setup; m.gf = gf; m.rho = rho;
% m.Urms = Urms; m.Urms_hf = Urms_hf; m.Urms_lf = Urms_lf; m.Urms_hft = Urms_hft; m.Urms_lft = Urms_lft;
% m.Um = Um; m.Umt = Umt; m.ua = ua; m.uaMF = uaMF;
% m.Asu = Asu; m.Sku = Sku; m.AsuMF = AsuMF; m.SkuMF = SkuMF;
% m.Cm = Cm; m.Cmt = Cmt;
% m.Sdzb = Sdzb;
% m.Cm4cm = Cm4cm;
% m.z = z;
% m.xhrms = xhrms; m.xurms = xurms; m.xum = xum; m.xswf = xswf;
% m.interval = interval; m.tMF = tMF;
% clear s;

addpath D:\d3d\toolbox\
wlsettings('D:\d3d\toolbox\');
dir = ['D:\documents\deltagoot\data_CD\FlumeWall\'];
profiles = [{[dir,par.test,'\',par.test,'_zt0000.tek']}; 
            {[dir,par.test,'\',par.test,'_zt0010.tek']}; 
            {[dir,par.test,'\',par.test,'_zt0020.tek']}; 
            {[dir,par.test,'\',par.test,'_zt0060.tek']};];
x = 0:1:220; z  =zeros(length(profiles),length(x));
for i = 1:length(profiles)
    bath = tekal('open',profiles{i});
    xz1 = bath.Field.Data;
    z(i,:) = interp1(xz1(:,1),xz1(:,2),x)-4.5;
end


%% simulation results
nam = [{'zb'};]; %{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'cc'};{'ceq'};{'Su'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'Tbore'};{'Beta'};{'R'};{'DR'};{'D'};{'kb'};{'c'};];

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

%% waveheights
% pwh(s,m,par,nt,nx,ny,xw,yw);

%% flows
% pu(s,m,par,nt,nx,ny,xw,yw);

%% sediment
% pst(s,m,par,nt,nx,ny,xw,yw);

%% profile evolution
% ppe(s,m,par,nt,nx,ny,xw,yw);
dx = 1;
hours = [1.0 2.04 6.0];
hours2 = [0.0 1.0 2.04 6.0];
indh = []; indh(1) = 1;
indh(2:4) = round(hours*3600/par.morfac/par.tint);

% compute erosion volumes above maximum still water level
Am = zeros(1,4); As = Am;
for i = 2:4;
    Am(i) = sum(max(0,z(1,:))-max(0,z(i,:)))*dx;
    As(i) = sum(max(0,s.zb(:,1))-max(0,s.zb(:,indh(i))))*dx;
end

figure(1); subplot(211); set(gca,'FontSize',10);
for i = 1:4
    plot(x,z(i,:),'k'); hold on;
    plot(xw(:,2),s.zb(:,indh(i)),'k--','LineWidth',1.5);
end
plot([0 xw(end,2)],[0 0],'k')
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([0 max(max(xw)) -4.5 2]);
print profiles.png -dpng

figure(2); subplot(211); set(gca,'FontSize',10);
for i = 1:4
    plot(x,z(i,:),'k'); hold on;
    plot(xw(:,2),s.zb(:,indh(i)),'k--','LineWidth',1.5);
end
plot([0 xw(end,2)],[0 0],'k');
xlabel('x [m]'); ylabel('z_{b} [m]'); axis([160 max(max(xw)) -2.0 2]);
print detailed_profiles.png -dpng

figure(3); set(gca,'FontSize',16);
plot(hours2,Am(1,:),'k-s','LineWidth',1.5,'MarkerFaceColor','k'); hold on;
plot(hours2,As(1,:),'k--s','LineWidth',1.5,'MarkerSize',12);
xlabel('t [hours]'); ylabel('A [m^3/m]');
print t_dune_erosion.png -dpng

%% movieplot
% figure;
% for i = nt-10:nt
%     plot(xw,s.zb(:,1),'k--'); hold on;
%     plot(xw,s.zb(:,i),'k');
%     plot(xw,s.zs(:,i),'b');
%     plot(xw,s.Hrms(:,i),'g');
%     plot(xw,s.ue(:,i),'r');
%     axis([0 230 -4.5 2]);
%     title(num2str(i));
%     pause; hold off;
% end

%% use output for other stuff

save([par.test,'.mat'],'par','s','m','xw');