function plotdata()
clear s m

%% coefficients and settings

%% measurements

%% profiles
addpath f:\d3d\toolbox\
wlsettings('f:\d3d\toolbox\');

prof = {'zt0000.tek','zt0018.tek','zt0040.tek','zt0193.tek'};

x = 0:1:230;
z = zeros(length(prof),length(x));
temp = ['F:\subversion_xbeach\trunk\testbed\data\Deltaflume_M1263_T3\'];
for i = 1:length(prof)
    bed = tekal('read',[temp,prof{i}]);
    xt = bed.Field.Data(:,1)';
    xt = [xt,230];
    zt = bed.Field.Data(:,2)';
    zt = [zt,zt(end)];
    z(i,:) = interp1(xt,zt,x)'-4.5;
end

%% simulation results
nam = [{'zb'};{'H'};{'zs'};];
% dimensions
fid = fopen('dims.dat','r'); 
temp = fread(fid,[4,1],'double'); 
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

%% profile evolution
ppe(s,m,par,nt,nx,ny,xw,yw);

else
% movieplot
figure;
for i = 1:nt
    plot(x, z(4,:)+4.5,'k-','LineWidth',1.5); hold on;
    plot(xw,s.zb(:,1),'k--'); 
    plot(xw,s.zb(:,i),'k-');
    plot(xw,s.zs(:,i),'b');
    plot(xw,s.H(:,i),'g');
    axis([0 230 0 7]);
    title(num2str(i));
    pause(0.2); hold off;
end

zsi = load('waterlevels.wls');
figure; 
plot(150+ts*2310,s.zs(1,:),'b--'); hold on;
plot(zsi(:,1),zsi(:,2),'b','LineWidth',1.5);


end
%% use output for other stuff

save([par.test,'.mat'],'par','s','m','xw');