function plotdata()
%clear variables
close all
[runid,testid,datadir]=testinfo

fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fip=fopen('point001.dat');
for i=1:nt*10
    data1(i,:)=fread(fip,[2 1],'double');
end
fclose(fip)
fip=fopen('point002.dat');
for i=1:nt*10
    data2(i,:)=fread(fip,[2 1],'double');
end
fclose(fip)
fip=fopen('point003.dat');
for i=1:nt*10
    data3(i,:)=fread(fip,[2 1],'double');
end
fclose(fip)
fip=fopen('point004.dat');
for i=1:nt*10
    data4(i,:)=fread(fip,[2 1],'double');
end
fclose(fip)
fip=fopen('point005.dat');
for i=1:nt*10
    data5(i,:)=fread(fip,[2 1],'double');
end
fclose(fip)

% tuba's values
l = 20;
g = 9.81;
ho = 0.1*(4./pi)*20;
H = .02*ho;
T = l/sqrt(g*ho);

load ([datadir 'etat1.dat'])
load ([datadir 'etat2.dat'])
load ([datadir 'etat3.dat'])
load ([datadir 'etat4.dat'])
load ([datadir 'etat5.dat'])
load ([datadir 'etamin.dat'])
load ([datadir 'etamax.dat'])
t=0:0.04:(length(etat1)-1)*0.04;
[i,j] = max(etat1);
corr = t(j)/T;
figure(100)
clf
subplot(121)
plot(t/T-corr,etat1/H+20,'r--')
hold on
plot(t/T-corr,etat2/H+15,'r--')
hold on
plot(t/T-corr,etat3/H+10,'r--')
hold on
plot(t/T-corr,etat4/H+5,'r--')
hold on
plot(t/T-corr,etat5/H,'r--')

grid on


subplot(122)

plot([1:length(etamax)]/length(etamax)*2,etamax/H,'r--')
hold on
plot([1:length(etamax)]/length(etamax)*2,etamin/H,'r--')

% mijn waarden
l = 8;
g = 9.81;
ho = 0.1*(4./pi)*l;
H = .02*ho;
T = l/sqrt(g*ho);
[i,j] = max(data5(:,2));
corr = data5(j)/T;
subplot(121)
plot(data5(:,1)/T-corr, data5(:,2)/H+20)
hold on

plot(data4(:,1)/T-corr, data4(:,2)/H+15)

plot(data3(:,1)/T-corr, data3(:,2)/H+10)

plot(data2(:,1)/T-corr, data2(:,2)/H+5)

plot(data1(:,1)/T-corr, data1(:,2)/H)

set(gca,'Box','on');

set(gca,'ytick',(-5:5:25))
vi=['-5';' 0';' 0';' 0';' 0';' 0';' 5'];
set(gca,'yticklabels',vi)

set(gca,'xtick',(-5:5:20))
vi=['-5';' 0';' 5';'10';'15';'20'];
set(gca,'xticklabels',vi)

xlabel('t/T');
ylabel('\zeta/H')



fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fid,[nx+1,ny+1],'double');
y=fread(fid,[nx+1,ny+1],'double');
fclose(fixy)
fiz=fopen('maxzs.dat','r');
fiy=fopen('minzs.dat','r');
fib=fopen('zb.dat','r');
for i=1:nt;
    maxzs=fread(fiz,[nx+1,ny+1],'double');
    minzs=fread(fiy,[nx+1,ny+1],'double');
    bathy=fread(fib,[nx+1,ny+1],'double');
end

for j=1:ny;
  zsmax(j)=maxzs(max(find(maxzs(:,j)>bathy(:,j)+0.005)),j);
  zsmin(j)=minzs(max(find(minzs(:,j)>bathy(:,j)+0.005)),j);
end

subplot(122)
hold on
plot([1:ny]/ny*2,zsmax/H)
plot([1:ny]/ny*2,zsmin/H)

xlabel('y/L_y');
ylabel('\zeta/H');
PN = fliplr(pwd);
[runid,R] = strtok(PN,'\'); runid = fliplr(runid);
[testid,R] = strtok(R,'\'); testid = fliplr(testid);
[dataid,R] = strtok(R,'\'); dataid = fliplr(dataid);
pname = [dataid '_' testid '_' runid '_fig' '1' '.jpg'];
eval(['print -djpeg ' pname]);

