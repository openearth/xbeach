clear s m

%% coefficients and settings
fid = fopen('..\data\datanam.txt');
Tnam = fgetl(fid);
nvar = str2num(fgetl(fid));
for i = 1:nvar
    line = fgetl(fid);
    [nam,remain] = strtok(line,' ');
    [dummy,remain] = strtok(remain,'=');
    par.(nam) = str2num(remain(2:end));
end
fclose(fid);
%% measurements
load(['..\data\',Tnam]);
m=s;
clear s;

%% simulation results
nam = [{'zb'};{'Hrms'};{'zs'};{'hh'};{'u'};{'ue'};{'urms'};{'cc'};{'ceq'};{'Su'};{'Fx'};{'dzav'};{'Fx'};{'ua'};{'uon'};{'uoff'};{'As'};{'Sk'};{'Tbore'};{'Beta'};{'R'};{'DR'};{'kb'};{'c'};];

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
% for i = 1:nt
%     plot(xw,s.zb(:,1),'k--'); hold on;
%     plot(xw,s.zb(:,i),'k');
%     plot(xw,s.zs(:,i),'b');
%     plot(xw,s.Hrms(:,i),'g');
%     plot(xw,s.ue(:,i),'r');
%     plot(m.x{end},m.z{end},'k-.');
%     axis([0 210 -4.5 1.5]);
%     title(num2str(i));
%     pause(0.01); hold off;
% end
