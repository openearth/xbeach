close all; clear all;

load hardlayer_T2.dep
load bed_T2.dep

dzg = 0.1;
dx = 1;
x = [0:1:length(bed_T2)-1]*dx;

temp = hardlayer_T2(1,:);
temp = min(temp,2);
ndsand = floor(temp/dzg);
gdist1 = zeros(20,length(ndsand));
for i = 1:length(ndsand); 
    gdist1(1:ndsand(i),i)=1; 
end;
gdist2 = 1-gdist1;

zl = fliplr([-19:1:0]*dzg);
zlm = repmat(zl,length(x),1)' + repmat(bed_T2(1,:),20,1);
xm = repmat(x,20,1);


figure; 
subplot(211);
plot(x,bed_T2(1,:),'r-','LineWidth',3); hold on; pcolor(xm,zlm,gdist1); colorbar;
axis([160 220 -4 3]);
subplot(212);
plot(x,bed_T2(1,:),'r-','LineWidth',3); hold on; pcolor(xm,zlm,gdist2); colorbar;
axis([160 220 -4 3]);

%% make XBeach file
gdist1XB = zeros(3*20,length(ndsand));
for i =1:20
   gdist1XB(i*3-2:i*3,:) = repmat(gdist1(i,:),3,1); 
   gdist2XB(i*3-2:i*3,:) = repmat(gdist2(i,:),3,1); 
end
    

save('gdist1.inp','gdist1XB','-ascii');
save('gdist2.inp','gdist2XB','-ascii');