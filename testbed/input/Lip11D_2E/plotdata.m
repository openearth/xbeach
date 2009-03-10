function plotdata()
close all
fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double')
nx=fread(fid,[1],'double')
ny=fread(fid,[1],'double')
fclose(fid)
fixy=fopen('xy.dat','r');
x=fread(fid,[nx+1,ny+1],'double');
y=fread(fid,[nx+1,ny+1],'double');
fclose(fixy)
fizs=fopen('zs.dat','r');
fizb=fopen('zb.dat','r');
fihrms=fopen('H.dat','r');
fiurms=fopen('urms.dat','r');
fiu=fopen('u.dat','r');
fiue=fopen('ue.dat','r');
it=0
zbt=zeros(nx+1,nt-300);
zst=zeros(nx+1,nt-300);
ut=zeros(nx+1,nt-300);
uet=zeros(nx+1,nt-300);
hrmst=zeros(nx+1,nt-300);
urmst=zeros(nx+1,nt-300);
for i=1:nt;
    zb=fread(fizb,[nx+1,ny+1],'double');
    zs=fread(fizs,[nx+1,ny+1],'double');
    hrms=fread(fihrms,[nx+1,ny+1],'double');
    urms=fread(fiurms,[nx+1,ny+1],'double');
    u=fread(fiu,[nx+1,ny+1],'double');
    ue=fread(fiue,[nx+1,ny+1],'double');
    if i==1
        z0=zb;
    end
    
    if mod(i,1)==0&i>300
        it=it+1;
        zbt(:,it)=zb(:,2);
        zst(:,it)=zs(:,2);
        ut(:,it)=u(:,2);
        uet(:,it)=ue(:,2);
        hrmst(:,it)=hrms(:,2);
        urmst(:,it)=urms(:,2);
    end
end;
fclose(fid)
xm=x(:,2);
zsm=mean(zst,2);
Hrms_hi=sqrt(mean(hrmst.^2,2));
t=[0:size(zbt,2)-1];nt=length(t);
f_cutoff=0.1;df=1/t(end);av_yes=0;

[wav_high , wav_low] = filterhjb(zst', f_cutoff , df , av_yes);
Hrms_lo=sqrt(8)*(std(wav_low));
Hrms_hi=sqrt(Hrms_hi.^2+8*(std(wav_high))'.^2);
Urms=sqrt(mean(urmst.^2,2));
%Urms_lo=(std(uet'))';
[wav_high , wav_low] = filterhjb(uet', f_cutoff , df , av_yes);
Urms_lo=(std(wav_low));
uemean=(mean(uet'))';
guls=3*mean(ut.*urmst.^2,2);
out=zeros(180,2)
out(:,1)=xm(1:180);
out(:,2)=zsm(1:180);
tekal('write','eta.tek',out);
out(:,2)=Hrms_hi(1:180);
tekal('write','Hrms_h_i.tek',out);
out(:,2)=Hrms_lo(1:180);
tekal('write','Hrms_l_o.tek',out);
out(:,2)=Urms(1:180);
tekal('write','urms.tek',out);
out(:,2)=Urms_lo(1:180);
tekal('write','urms_l_o.tek',out);
out(:,2)=guls(1:180);
tekal('write','guls.tek',out);
out(:,2)=uemean(1:180);
tekal('write','u_m_e_a_n.tek',out);
tek(:,1)=xm(:);
tek(:,2)=zbt(:,nt)-zbt(:,1);
tekal('write','sedero.tek',tek);
showbank
plotpro
x1=x(:,2);
figure(3);
subplot(411);
pcolor(t,x1,hrmst);shading interp;caxis([0 2])
subplot(412);
pcolor(t,x1,zst-mean(mean(zst)));shading interp;caxis([-1 1 ])
subplot(413);
pcolor(t,x1,uet);shading interp;caxis([-1 1])
subplot(414);
pcolor(t,x1,urmst);shading interp;caxis([-1 1])
