function plotdata()
clear s m

%% coefficients and settings
[runid,testid,datadir]=testinfo
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
m.Hrms = Hrms; m.Hrms_hf = Hrms_hf; m.Hrms_lf = Hrms_lf; m.Hrms_in = Hrms_in; m.Hrms_out = Hrms_out;
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
m.um = um; m.zum = zum; m.cm = cm; m.zcm = zcm;
clear s;

%% simulation results
% nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'ccg'};{'ceqg'};{'Sug'};{'Fx'};{'dzav'};{'Fx'};{'urep'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'BR'};{'Tbore'};{'R'};{'DR'};{'D'};{'kb'};{'c'};]; %
nam = [{'zb'};{'H'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'ccg'};{'ceqg'};{'Sug'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'BR'};{'Tbore'};{'R'};{'DR'};{'D'};{'kb'};{'c'};]; %
% dimensions
fid = fopen('dims.dat','r'); 
temp = fread(fid,[4,1],'double'); 
nt = temp(1);
nx = temp(2)+1;
ny = temp(3)+1;
kmax = temp(4);
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

if kmax>1
    
% read 3d arrays....
x = [170:5:205];
for i = 1:length(x)
    ind(i) = find(xw(:,2)==x(i));
end
nam2 = {'sig','cuq3d'};
for j = 1:length(nam2)
    temp = zeros(length(ind),kmax,nt);
    fid = fopen([nam2{j},'.dat'],'r');
    for i=1:nt
        for ii=1:kmax
            ttemp = fread(fid,[nx,ny],'double');
            temp(:,ii,i) = ttemp(ind,2);               
        end
    end
    fclose(fid);
    s.(nam2{j}) = temp;
end

% read 4d arrays
nam3 = {'veloc'};
for j = 1:length(nam3)
    temp = zeros(length(ind),kmax,3,nt);
    fid = fopen([nam3{j},'.dat'],'r');
    for i=1:nt
        for ii=1:3
            for iii=1:kmax
                ttemp = fread(fid,[nx,ny],'double');
                temp(:,iii,ii,i) = ttemp(ind,2);
            end
        end
    end
    fclose(fid);
    s.(nam3{j}) =temp;
end

figure;
for i = 1:nt
    for j = 1:length(ind)
        subplot(2,4,j);
        htemp = (1+s.sig(j,:,i))*s.hh(ind(j),i);
        %plot(squeeze(s.veloc(j,:,1,i)),htemp,'b'); hold on;
        %plot(squeeze(s.veloc(j,:,3,i)),htemp,'b'); hold on;
        plot(par.rhos*squeeze(s.cuq3d(j,:,i)),htemp,'g');
        hold off;
    end
    pause(0.1); 
end

end

det = 1;
if det == 1;
%% waveheights
pwh(s,m,par,nt,nx,ny,xw,yw);

%% flows
pu(s,m,par,nt,nx,ny,xw,yw);

%% sediment
pst(s,m,par,nt,nx,ny,xw,yw);

%% profile evolution
ppe(s,m,par,nt,nx,ny,xw,yw);

else
% movieplot
figure;
for i = 1:nt
    plot(xw,s.zb(:,1),'k--'); hold on;
    plot(xw,s.zb(:,i),'k-');
    plot(xw,2650*s.ccg(:,i),'c-.');
    plot(xw,s.zs(:,i),'b');
    plot(xw,s.H(:,i),'g');
    plot(xw,s.ue(:,i),'r');
    axis([0 230 -4.5 2]);
    title(num2str(i));
    pause(0.1); hold off;
end

end
%% use output for other stuff

save([par.test,'.mat'],'par','s','m','xw');