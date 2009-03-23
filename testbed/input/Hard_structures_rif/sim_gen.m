close all; clear all;

% Set up simple simulation to test effect rif

% profile withou rif:
zref = [-20 -3 0 3 15 15];
dzref = [17 3 3 12];
sref = [180 70 20 3];
xref = [0 cumsum(dzref.*sref) sum(dzref.*sref)+100];
 
dx = 1;
x = 0:dx:xref(end);
z1 = interp1(xref,zref,x);

% profile with rif:
xrifb = x(find(z1==-10));
zrifb = 10-2; % keep rif 2 meters below mean sea level so it is expected to remain below 
rifs = 5;
xrif = [xrifb-15-8*rifs xrifb-15 xrifb+15 xrifb+15+8*rifs];
zrif = [z1(find(x==xrif(1))) z1(find(x==xrif(2)))+8 z1(find(x==xrif(3)))+8 z1(find(x==xrif(4)))];
xt = xrif(1):1:xrif(end);
z2t = interp1(xrif,zrif,xt);
ind1 = find(x==xrif(1));
ind2 = find(x==xrif(end));

z2 = z1;
z2(ind1:ind2) = z2t;

sandlayer = z2*0.0;
sandlayer(:) = 20;
sandlayer(ind1-round(50/dx):ind2+round(50/dx)) = 0; 
zsl = z2-20;
zsl(ind1-round(50/dx):ind2+round(50/dx)) = z2(ind1-round(50/dx):ind2+round(50/dx));

% interpolate all on grid with dx = 5m
xi = 0:5:x(end);
z1i = interp1(x,z1,xi);
z2i = interp1(x,z2,xi);
sli = interp1(x,sandlayer,xi);
zsli = interp1(x,zsl,xi);

figure; 
subplot(2,1,1); plot(xref,zref,'k'); hold on; plot(x,z2,'b','LineWidth',1.5); %plot(x,z1,'b--','LineWidth',1.5);
plot(x,zsl,'r','LineWidth',1.5);
subplot(2,1,2); hold on; plot(xi,z2i,'b','LineWidth',1.5); %plot(xi,z1i,'b--','LineWidth',1.5);
plot(xi,zsli,'r','LineWidth',1.5);

z1ip = repmat(z1i,3,1);
z2ip = repmat(z2i,3,1);
slip = repmat(sli,3,1);

save('bed.dep','z2ip','-ascii');
save('bed_no_rif.dep','z1ip','-ascii');
save('hardlayer.dep','slip','-ascii');

