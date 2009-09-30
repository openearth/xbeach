clear s m

%% coefficients and settings
fid = fopen('..\data\datanam.txt');
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
addpath D:\d3d\toolbox\
wlsettings('D:\d3d\toolbox\');

x = 0:1:220;
if strcmp(par.test,'T1')==1
    prof = {'zt0000.tek','zt0006.tek','zt0013.tek','zt0035.tek','zt0095.tek'};
else
    prof = {'zt0000.tek','zt0013.tek','zt0035.tek','zt0095.tek'};
end
for i = 1:length(prof)
    location = ['E:\data\testbank\H0298\H0298-Test-',par.test(2),'\',prof{i}];    
    bath = tekal('open',location);
    xz = bath.Field.Data;
    z{i} = interp1(xz(:,1),xz(:,2),x);
    z{i}(isnan(z{i})) = xz(end,2);
end

% load hard layer
hardlayer = load('hardlayer_T2.dep');

%% simulation results
nam = [{'zb'};{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'cc'};{'Su'};{'dzav'};];

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

zhard = s.zb(:,1)-hardlayer(2,:)';

%% profile evolution
ind = round([0 0.6 1.3 3.5 9.6]*3600/par.morfac/par.tint); ind(1) = 1;
indh = find(zhard==max(zhard));
figure; set(gca,'FontSize',16)
plot(xw(:,2),s.zb(:,1),'k-','LineWidth',1.5); hold on; 
for i = 1:length(z)
    plot(xw(:,2),s.zb(:,ind(i)),'k--','LineWidth',1.5); hold on; plot(x,z{i}-par.zs0,'k');
end
plot(xw(1:indh,2),zhard(1:indh),'r-','LineWidth',6);
plot([0 210],[0 0],'k:');
axis([150 210 -2 2.5]); xlabel('x [m]'); ylabel('z_b [m]');
print profile_evolution.png -dpng

t = (0:1:nt-1)*par.tint*par.morfac/3600;
Tsim = t(end)*3600;
rhos  = 2650;
dz = s.zb(:,end)-s.zb(:,1);
s.Sdzbtot = (1-par.np)*flipud(cumsum(flipud(dz)))*(xw(2,2)-xw(1,2))/Tsim;
s.Sdzavtot = (1-par.np)*flipud(cumsum(flipud(s.dzav(:,end))))*(xw(2,2)-xw(1,2))/Tsim;

% figure; 
% subplot(211);
% plot(xw,mean(s.Su,2),'k-.'); hold on;
% plot(xw,mean(s.Sdzavtot,2),'k--');
% plot(xw,mean(s.Sdzbtot,2),'k-','LineWidth',1.5);
% axis([140 210 -5.0e-5 0.5e-5]);
% subplot(212);
% plot(xw,s.zb(:,1),'k-'); hold on;
% plot(xw,hardlayer','r-','LineWidth',1.5); hold on;
% axis([140 210 -2.5 2.2]);

% % movieplot
% figure;
% for i = 1:nt
%    plot(xw,s.zb(:,1),'k--'); hold on;
%    plot(xw,s.zb(:,i),'k-s');
%    plot(xw,s.zs(:,i),'b');
%    plot(xw,s.Hrms(:,i),'g');
%    plot(xw,s.ue(:,i),'r');
%    plot(xw,par.rhos*s.cc(:,i),'m','LineWidth',1.5);
%    plot([0 230],[par.topr par.topr],'k--');
%    axis([0 230 -4.5 2]);
%    title(num2str(i));
%    pause(0.05); hold off;
% end
