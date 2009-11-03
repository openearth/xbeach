function plotdata()
clear all;close all
[runid,testid,datadir]=testinfo
%% simulation results
nam = [{'zb'};{'zs'};{'u'};{'ue'};{'H'};]; %{'dzav'};

% dimensions
fid = fopen('dims.dat','r'); 
temp = fread(fid,[7,1],'double'); 
nt = temp(1);
nx = temp(2)+1;
ny = temp(3)+1;
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

% load hard layer new approach...
dzg = 0.05;
load('gdist1.inp');
gds = gdist1(2:3:end-1,:);
dz = sum(gds,1)*dzg;

zhl = s.zb(:,1)-dz';

%% movieplot
figure;
for i =1:nt
    plot(xw,s.zb(:,1),'k--'); hold on;
    plot(xw,zhl,'r','LineWidth',4);
    % plot(xgr,zgr.*bgr-3.8,'ro');%c{bgr+1});
    plot(xw,s.zb(:,i),'k','LineWidth',1.5);
    plot(xw,s.zs(:,i),'b');
    plot(xw,s.H(:,i),'r');
    %plot(xw,2650*s.Su(:,i),'g');
    plot(xw(1:end-1,:),s.ue(1:end-1,i),'r');
    % plot(xw,2650*s.cc(:,i),'g');
    axis([0 50 -2 2]);
    title(num2str(i));
    name = ['avi\hard_',num2str(i,'%04.0f'),'.png'];
    print('-dpng','-r300',name);
    pause(0.01); hold off;
end

print endprofile.png -dpng
pname = ['..\..\report\',testid '_' runid '_fig1' '.jpg'];
eval(['print -djpeg ' pname]);
