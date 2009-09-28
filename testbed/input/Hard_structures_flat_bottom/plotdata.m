function plotdata()
clear all;close all
[runid,testid,datadir]=testinfo
%% simulation results
nam = [{'zb'};{'zs'};{'u'};{'ue'};{'H'};]; %{'dzav'};

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

% get hard layer...
% bgr = load('graindist.grd');
% zgr = zeros(size(bgr));
% for i =1:20;
%     zgr(i,:) = s.zb(:,1)+(i-1 )*0.2;
% end  
% xgr = repmat(xw(:,2),1,20)';

% load hard layer new approach...
dzhl = load('hardlayer.dep');
zhl = s.zb(:,1)-dzhl';

%% movieplot
figure;
for i =nt:nt
    plot(xw,s.zb(:,1),'k--'); hold on;
    plot(xw,zhl,'r','LineWidth',4);
    % plot(xgr,zgr.*bgr-3.8,'ro');%c{bgr+1});
    plot(xw,s.zb(:,i),'k');
    plot(xw,s.zs(:,i),'b');
    plot(xw,s.H(:,i),'r');
    %plot(xw,2650*s.Su(:,i),'g');
    plot(xw(1:end-1,:),s.ue(1:end-1,i),'r');
    % plot(xw,2650*s.cc(:,i),'g');
    axis([0 50 -2 2]);
    title(num2str(i));
    pause(0.1); hold off;
end

print endprofile.png -dpng
pname = ['..\..\report\',testid '_' runid '_fig1' '.jpg'];
eval(['print -djpeg ' pname]);
