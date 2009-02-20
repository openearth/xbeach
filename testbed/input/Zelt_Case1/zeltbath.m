% zelt bathymetry from Tuba, 1/20/97

figure(1)
clear
clg
orient portrait
%l=20.;
l=8;
g=9.81;
Ls = l;
dx=.125;

x=(-3.92.*l/pi:dx:4.66*l/pi);
%x=(0:.5:3*l);
y=(0:dx:2*l);
[X,Y]=meshgrid(x,y);

h0=0.4*l/pi;
z0=cos(pi*y/l);
z1=zeros(size(y))+h0;
xt=0.;
nx=length(x);
ny=length(y);

T=sqrt(g*h0)/l


for i=1:nx
if x(i) >= xt, h(:,i)=h0-(0.4*x(i))./(3.-z0(:));
else, h(:,i)=z1(:);
end
end



mesh(x,y,-h)
view(-135,30)
axis([-5.5*l/pi 4.5*l/pi 0 2.2*l -1.*h0 2*h0])
hold on



xn=3*l/pi-l/pi*z0;
xn(ny+1)=x(1);
xn(ny+2)=x(1);
xn(ny+3)=3*l/pi-l/pi*z0(1);
yn=y;
yn(ny+1)=y(ny);
yn(ny+2)=y(1);
yn(ny+3)=y(1);
h1=zeros(1,ny+3);

plot3(xn,yn,h1+0.01,'linewidth',1)

hold on

% x-arrow
%plot3([-10 -7],[0 0],[0.01 0.01],'linewidth',1)

%y-arrow
%plot3([-10 -10],[0 3],[0.01 0.01],'linewidth',1)

%z-arrow
plot3([-10 -10],[0 0],[0 0.4],'linewidth',1)

% arrow heads
plot3([-7 -7.5],[0  0.5],[0.01 0.01],'linewidth',1)
plot3([-7 -7.5],[0 -0.5],[0.01 0.01],'linewidth',1)

plot3([-10 -9.5],[3 2.5],[0.01 0.01],'linewidth',1)
plot3([-10 -10.5],[3 2.5],[0.01 0.01],'linewidth',1)

plot3([-10 -10.5],[0 0],[.4 0.3],'linewidth',1)
plot3([-10 -9.5],[0 0],[0.4 0.3],'linewidth',1)

text(-5,-1,0,'x')
text(-11,4,.4,'y')
text(-10,0,.6,'z')

%measuring stations
plot3(2*l/pi,0,0.01,'.','markersize',20);
plot3(l/pi*(3-.5*sqrt(2)),.25*l,0.01,'.','markersize',20);
plot3(l/pi*3,.5*l,0.01,'.','markersize',20);
plot3(l/pi*(3+.5*sqrt(2)),.75*l,0.01,'.','markersize',20);
plot3(l/pi*4,l,0.01,'.','markersize',20);

measstaty = floor( [ 0 .25*l .5*l .75*l l]/dx+1)
for ii=1:length(measstaty)
    [jj]=max(find(h(measstaty(ii),:)>0.0))
end


text(2*l/pi+1,0,.5,'0')
text(l/pi*(3-.5*sqrt(2))+1,.25*l,.5,'0.25')
text(l/pi*3+1,.5*l,.5,'0.5')
text(l/pi*(3+.5*sqrt(2))+1,.75*l,.5,'0.75')
text(l/pi*4+1,l,.5,'1')
text(l/pi*4+1,l,1,'y/Ly = ')


xlabel('cross shore (m)');
ylabel('along shore (m)');
zlabel('depth (m)');

% Ls-arrow
plot3([-10 0],[17 17],[-1 -1],'linewidth',1)
% arrowheads
plot3([0 -.5],[17  17.5],[-1 -1],'linewidth',1)
plot3([0 -.5],[17  16.5], [-1 -1],'linewidth',1)
plot3([-10 -9.5],[17  17.5],[-1 -1],'linewidth',1)
plot3([-10 -9.5],[17  16.5], [-1 -1],'linewidth',1)
text(-5,19,-1,'Ls');

% hs-arrow
plot3([-10 -10],[4 4],[0 -1],'linewidth',1)
% arrowheads
plot3([-10 -10],[4 4.5],[-1 -.9],'linewidth',1)
plot3([-10 -10],[4 3.5], [-1 -.9],'linewidth',1)
plot3([-10 -10],[4 4.5],[0 -.1],'linewidth',1)
plot3([-10 -10],[4 3.5], [0 -.1],'linewidth',1)
text(-10,3,-.4,'hs');

% Ly-arrow
plot3([-11 -11],[0 16],[-1 -1],'linewidth',1)
% arrowheads
plot3([-11 -11.5],[0 .5],[-1 -1],'linewidth',1)
plot3([-11 -10.5],[0 .5], [-1 -1],'linewidth',1)
plot3([-11 -11.5],[16  15.5],[-1 -1],'linewidth',1)
plot3([-11 -10.5],[16  15.5], [-1 -1],'linewidth',1)
text(-12.5,8,-1,'2 Ly');



axis('off')



