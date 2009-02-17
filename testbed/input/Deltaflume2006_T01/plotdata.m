clear s m
[runid,testid,datadir]=testinfo

%% coefficients and settings
fid = fopen([datadir 'datanam.txt']);
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
load([datadir,Tnam]);
m.Hrms = Hrms; m.Hrms_hf = Hrms_hf; m.Hrms_lf = Hrms_lf;
m.setup = setup; m.gf = gf; m.rho = rho;
m.Urms = Urms; m.Urms_hf = Urms_hf; m.Urms_lf = Urms_lf; m.Urms_hft = Urms_hft; m.Urms_lft = Urms_lft;
m.Um = Um; m.Umt = Umt; m.ua = ua; m.uaMF = uaMF;
m.Asu = Asu; m.Sku = Sku; m.AsuMF = AsuMF; m.SkuMF = SkuMF;
m.Cm = Cm; m.Cmt = Cmt;
m.Sdzb = Sdzb;
m.Cm4cm = Cm4cm;
m.z = z;
m.xhrms = xhrms; m.xurms = xurms; m.xum = xum; m.xswf = xswf;
m.interval = interval; m.tMF = tMF;
clear s;

%% simulation results
nam = [{'zb'};{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'cc'};{'ceq'};{'Su'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'Tbore'};{'Beta'};{'R'};{'DR'};{'D'};{'kb'};{'c'};];

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
pwh(s,m,par,nt,nx,ny,xw,yw);

%% flows
pu(s,m,par,nt,nx,ny,xw,yw);

%% sediment
pst(s,m,par,nt,nx,ny,xw,yw);

%% profile evolution
ppe(s,m,par,nt,nx,ny,xw,yw);

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