a=load('fig1\ms1.xyz');
t1=a(:,1)/1000;w1=a(:,2)/10000;
for i=1:12;w1(2:end-1)=.25*w1(1:end-2)+.5*w1(2:end-1)+.25*w1(3:end);end
a=load('fig1\ms2.xyz');
t2=a(:,1)/1000;w2=a(:,2)/10000;
a=load('fig1\ms3.xyz');
t3=a(:,1)/1000;w3=a(:,2)/10000;
a=load('fig1\ms4.xyz');
t4=a(:,1)/1000;w4=a(:,2)/10000;
a=load('fig1\ms5.xyz');
t5=a(:,1)/1000;w5=a(:,2)/10000;
figure(4);
subplot(211)
plot(t1,w1,t2,w2,t3,w3,t4,w4,t5,w5,'linewidth',2)
legend('MS1','MS2','MS3','MS4','MS5')
a=load('fig2\ms2.xyz');
t2=a(:,1)/1000;v2=a(:,2)/10000;
a=load('fig2\ms3.xyz');
t3=a(:,1)/1000;v3=a(:,2)/10000;
a=load('fig2\ms4.xyz');
t4=a(:,1)/1000;v4=a(:,2)/10000;
a=load('fig2\ms5.xyz');
t5=a(:,1)/1000;v5=a(:,2)/10000;
subplot(212)
plot(t2,v2,t3,v3,t4,v4,t5,v5,'linewidth',2)
legend('MS2','MS3','MS4','MS5');hold on
fid = fopen('zs0input.dat','wt');
for i=1:length(t1);
    fprintf(fid,'%8.2f  %8.2f %8.2f\n',t1(i)*60,w1(i),0.7);
end
fclose(fid);
plot(t,zst(35,:),t,zst(39,:),t,zst(103,:),t,zst(109,:),'linewidth',2)
%legend('MP2','MP3','MP4','MP6')
