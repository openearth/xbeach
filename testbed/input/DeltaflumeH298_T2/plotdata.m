function plotdata();

clear s m

%% coefficients and settings
[runid,testid,datadir]=testinfo
fid = fopen([datadir 'datanam.txt']);
nvar = str2num(fgetl(fid));
for i = 1:nvar
    line = fgetl(fid);
    [nam,remain] = strtok(line,' ');
    [dummy,remain] = strtok(remain,'=');
    par.(nam) = str2num(remain(2:end));
end
par.(nam) = remain(end-1:end);
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
    location = ['F:\Data\H0298\H0298-Test-',par.test(end),'\',prof{i}];    
    bath = tekal('open',location);
    xz = bath.Field.Data;
    z{i} = interp1(xz(:,1),xz(:,2),x);
    z{i}(isnan(z{i})) = xz(end,2);
end

% load hard layer
hardlayer = load('hardlayer_T2.dep');

%% simulation results
nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'dzav'};];
nam2 = [{'ccg'};{'Susg'};{'Subg'};];

% dimensions
fid = fopen('dims.dat','r'); 
temp = fread(fid,[7,1],'double'); 
nt = temp(1);
nx = temp(2)+1;
ny = temp(3)+1;
ntheta = temp(4);
kmax = temp(5);
ngd = temp(6);
nd = temp(7);
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

for j = 1:length(nam2)
    temp = zeros(nx,ngd,nt);
    fid = fopen([nam2{j},'.dat'],'r');
    for i=1:nt
        for ii=1:ngd
            ttemp = fread(fid,[nx,ny],'double');
            temp(:,ii,i) = ttemp(:,2);               
        end
    end
    fclose(fid);
    s.(nam2{j}) = temp;
end


zhard = s.zb(:,1)-hardlayer(2,:)';

%% profile evolution
ind = round([0 0.6 1.3 3.5 9.6]*3600/par.morfac/par.tint); ind(1) = 1;
indh = find(zhard==max(zhard));
figure; set(gca,'FontSize',16)
plot(xw(:,2),s.zb(:,1),'k-','LineWidth',1.5); hold on;
for i = 2:length(z)
    plot(xw(:,2),s.zb(:,ind(i)),'k--','LineWidth',1.5); hold on; plot(x,z{i}-par.zs0,'k-','LineWidth',1.5);
end
plot(xw(1:indh,2),zhard(1:indh),'r','LineWidth',4);
plot(xw(:,2),max(s.zs'),'b--','LineWidth',1.5);
plot(xw(:,2),min(s.zs'),'b--','LineWidth',1.5);
plot([0 210],[0 0],'k:','LineWidth',2);
axis([150 210 -2 2.5]); xlabel('x [m]'); ylabel('z_b [m]');
pname = ['..\..\report\',testid '_' runid '_profile_evolution' '.jpg'];
eval(['print -djpeg ' pname]);

h = s.zs-s.zb;
t = (0:1:nt-1)*par.tint*par.morfac/3600;
Tsim = t(end)*3600;
rhos  = 2650;
% transform Su to Sz
Ssus1 = zeros(size(s.Susg));
Ssus1(1,:,:) = s.Susg(1,:,:);
Ssus1(2:end,:,:) = s.Susg(1:end-1,:,:); %0.5*(s.Su(1:end-1,:)+s.Su(2:end,:));
Sbed1 = zeros(size(s.Subg));
Sbed1(1,:,:) = s.Subg(1,:,:);
Sbed1(2:end,:,:) = s.Subg(1:end-1,:,:); %0.5*(s.Su(1:end-1,:)+s.Su(2:end,:));
dz = s.zb(:,end)-s.zb(:,1);
dx = zeros(1,length(xw));
dx(2:end) = diff(xw(:,1));
s.Sdzbtot = (1-par.np)*flipud(cumsum(flipud(dz))).*dx'/Tsim;
s.Sdzavtot = (1-par.np)*flipud(cumsum(flipud(s.dzav(:,end)))).*dx'/Tsim;
s.Sustot = mean(Ssus1,3);
s.Subtot = mean(Sbed1,3);
s.Sutot = s.Subtot+s.Sustot

% sediment transports from undertow:
% s.ut  = s.ue-s.u;
% Slongwaves = s.u.*s.cc.*h;
% Sundertow = s.ut.*s.cc.*h;
% Swaveassym = s.ua.*s.cc.*h;

figure; set(gca,'FontSize',16);
plot(xw(:,2),mean(s.ue,2),'k','LineWidth',2);
axis([150 210 -0.75 0.05]);
xlabel('x [m]'); ylabel('U_m [m/s]');

figure; set(gca,'FontSize',16);
plot(xw(:,2),mean(squeeze(s.ccg(:,1,:))*rhos,2),'k','LineWidth',2);
axis([150 210 0 5]);
xlabel('x [m]'); ylabel('c [g/l]');

figure;
plot(xw,s.Sdzbtot,'k--','LineWidth',1.5); hold on; 
plot(xw,s.Sdzavtot,'r--','LineWidth',1.5); plot(xw,s.Subtot(:,1),'r:','LineWidth',1.5); plot(xw,s.Sustot(:,1),'r-.','LineWidth',1.5); plot(xw,s.Sdzavtot+s.Sutot(:,1),'r');
xlabel('x [m]'); ylabel('S_{mean} [m^3/m/s]');

% axis([165 max(max(xw)) -4E-4 1E-4]);



%movieplot
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
%    pause(0.1); hold off;
% end

end