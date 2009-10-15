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
hardlayer = load('hardlayer_T3.dep');

%% simulation results
nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'v'};{'ve'};{'urms'};{'dzav'};];
nam2 = [{'ccg'};{'Susg'};{'Subg'};{'ceqbg'};{'ceqsg'};];

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

%% Movieplot
ind{1} = 1:1:nx;
nam3 = {'dzbed'};       dim3 = [nd; length(ind{1});];
nam4 = {'pbbed'};       dim4 = [nd; ngd; length(ind{1});];

for j = 1:length(nam3)
    temp = zeros(dim3(2,j),dim3(1,j),nt);
    fid = fopen([nam3{j},'.dat'],'r');
    for i=1:nt
        for ii=1:dim3(1,j)
            ttemp = fread(fid,[nx,ny],'double');
            temp(:,ii,i) = ttemp(ind{j},2);               
        end
    end
    fclose(fid);
    s.(nam3{j}) = temp;
end

for j = 1:length(nam4)
    temp = zeros(dim4(3,j),dim4(1,j),dim4(2,j),nt);
    fid = fopen([nam4{j},'.dat'],'r');
    for i=1:nt
        for ii=1:dim4(2,j)
            for jj=1:dim4(1,j)
                ttemp = fread(fid,[nx,ny],'double');
                temp(:,jj,ii,i) = ttemp(ind{j},2);    
            end
        end
    end
    fclose(fid);
    s.(nam4{j}) = temp;
end

figure;
ca=[-0.2 1.2];
cmap=colormap;
c1r=interp1(ca,[1 length(cmap)],1);
c1=cmap(round(c1r),:);
c0r=interp1(ca,[1 length(cmap)],0);
c0=cmap(round(c0r),:);
i=1;
temp = zeros(length(xw),nd+1);
temp(:,2:end) = cumsum(s.dzbed(:,:,i),2);
zlayer = repmat(s.zb(:,i),1,nd+1)-temp;
xt = repmat(xw(:,2),1,nd+1);
temp(:,2:end) = s.pbbed(:,:,1,i);
temp(:,1) = temp(:,2);
p1=    pcolor(xt,zlayer,temp); view(2); colorbar; caxis(ca); hold on; 
p2=    plot(xw,s.zb(:,1),'k--','LineWidth',1.5); 
p3=    plot(xw,s.zb(:,i),'k-','LineWidth',1.5);
p4=    plot(xw,s.zs(:,i),'b','LineWidth',1.5);
p5=    plot(xw,2650*s.ccg(:,1,i),'color',c1,'LineWidth',1.5);
p6=    plot(xw,2650*s.ccg(:,2,i),'color',c0,'LineWidth',1.5);
p7=    plot(xw,s.H(:,i),'c','LineWidth',1.5);
p8=    plot(xw,s.ue(:,i),'g','LineWidth',1.5);
axis([100 220 -4 7]);
for i =250:nt
    temp = zeros(length(xw),nd+1);
    temp(:,2:end) = cumsum(s.dzbed(:,:,i),2);
    zlayer = repmat(s.zb(:,i),1,nd+1)-temp;
    temp(:,2:end) = s.pbbed(:,:,1,i);
    temp(:,1) = temp(:,2);
    set(p1,'ydata',zlayer,'cdata',temp);
    % set(p1,'ydata',zlayer(:),'cdata',temp(:));
    set(p3,'ydata',s.zb(:,i));
    set(p4,'ydata',s.zs(:,i));
    set(p5,'ydata',2650*s.ccg(:,1,i));
    set(p6,'ydata',2650*s.ccg(:,2,i));
    set(p7,'ydata',s.H(:,i));
    set(p8,'ydata',s.ue(:,i));
    title(num2str(i));
    pause(0.1); hold off;    
end


%% continuity: RETE STOER!
morfac = 10;
morstart = 1000;
tint = 100;
ncorr = morstart/tint+1;
dt = 1;
ind = 10;
xu(1:nx-1) = 0.5d0*(xw(1:nx-1,2)+xw(2:nx,2));
dx = xu(ind:end)-xu(ind-1:end-1);
Vstart = sum(dx'.*(s.zb(ind:end-1,1)-min(s.zb(ind:end-1,1))));
% step 1: integrate profiles over flume length as function of time....
DVbed = sum(repmat(dx,nt,1)'.*s.zb(ind:end-1,:),1)-sum(dx'.*s.zb(ind:end-1,1),1);
plot(ts,DVbed/Vstart); hold on;
% step 2: computre sediment volume in water column over time....
Vsus = -sum(repmat(dx,nt,1)'.*squeeze(sum(s.ccg(ind:end-1,:,:),2)).*s.hh(ind:end-1,:)/0.6,1);
plot(ts,(Vsus-Vsus(1))*morfac/Vstart,'r');
% step 3: compute sediment volume that leaves the model through the
% offshore boundary....
DVbound = tint*squeeze(cumsum((sum(s.Susg(ind,:,:)+s.Subg(ind,:,:),2))/0.6,3));
% step 4: correct for sediemnt suspension before t=morstart
Vsuscorr = (Vsus-Vsus(ncorr))*morfac;
Vsuscorr(1:ncorr) = 0;
plot(ts,DVbound/Vstart/0.6,'g--');
% step 5: check continuity
plot(ts,(DVbed-Vsuscorr-DVbound')/Vstart,'k','LineWidth',1.5);

%%
figure;
borders = [0 100 150 221];
subplot(311); pcolor(squeeze(s.ccg(:,1,:))); axis(borders); caxis([0 0.01]); colorbar;
subplot(312); pcolor(squeeze(s.ceqsg(:,1,:))); axis(borders); caxis([0 0.01]); colorbar;
subplot(313); pcolor(squeeze(s.ceqbg(:,1,:))); axis(borders); caxis([0 0.01]); colorbar;

figure;
subplot(211); pcolor(squeeze(s.ue)); axis(borders); caxis([-4 4]); colorbar;
subplot(212); pcolor(squeeze(s.ve)); axis(borders); caxis([-4 4]); colorbar;

figure;
subplot(211); pcolor(squeeze(s.u)); axis(borders); caxis([-6 6]); colorbar;
subplot(212); pcolor(squeeze(s.v)); axis(borders); caxis([-4 4]); colorbar;

figure
DZ = diff(s.zb');
pcolor(DZ'); caxis([-1 1]); colorbar; axis(borders); %


%% profile evolution
ind = round([0 0.6 1.3 3.5 9.6]*3600/par.morfac/par.tint); ind(1) = 1;
indh = find(zhard==max(zhard));
figure; set(gca,'FontSize',16)
plot(xw(:,2),s.zb(:,1),'k-','LineWidth',1.5); hold on;
for i = 2:length(z)
    plot(xw(:,2),s.zb(:,ind(i+1)),'k--','LineWidth',1.5); hold on; plot(x,z{i}-par.zs0,'k-','LineWidth',1.5);
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


end