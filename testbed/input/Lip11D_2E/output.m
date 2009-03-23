function [it,s]=output(it,s,par);
if  mod(par.t,par.tint)==0;
    if it==0;
        load bot2e;
        s.xp=xp;s.zp=zp;
    end
    figure(1);
    %     pcolor(s.x,s.y,s.zs);shading interp;caxis([-.5 2]);colorbar;axis equal;
    %    hold on;quiver(s.x,s.y,s.u,s.v);drawnow;hold off
    %subplot(211);
    plot(s.x,s.zb,s.x,s.zs,'o-',s.x,s.H,s.x,(s.ceq.*max(s.h,0)),'linewidth',1);%axis([-125 50 -.5 .5])
    hold;
    plot(s.x,s.zb,'o-','linewidth',2);%axis([-125 50 -.5 .5]);
%    plot(s.x,s.zb-s.sedero,'g--','linewidth',1);%axis([-125 50 -.5 .5])
    plot(s.x,s.zb0,'g--','linewidth',2);%axis([-125 50 -.5 .5])
    plot(s.xp,s.zp(:,2),'k:','linewidth',2);
    plot(s.x,s.ueu,'r','linewidth',2);
    hold;
    title(num2str(par.t));
    axis([150 200 -1 6]);
    if 1;par.t>par.tstart; 
       it=it+1
       s.zbt(:,it)=s.zb(:,2);
       s.ht(:,it)=s.h(:,2);
       s.zst(:,it)=s.zs(:,2);
       s.ut(:,it)=s.u(:,2);
       s.uet(:,it)=s.u(:,2)-s.ust(:,2);
       s.Ht(:,it)=s.H(:,2);
       s.urmst(:,it)=s.urms(:,2);
       s.Fxt(:,it)=s.Fx(:,2);
       s.Sxxt(:,it)=s.D(:,2);
	   s.Sxyt(:,it)=s.Sxy(:,2);
    end
    %subplot(212);plot(s.x,s.uu,'linewidth',1);%axis([100 150 -1 1])
    %hold on;
    %subplot(313);plot(s.x,s.vv,'o-','linewidth',1);%axis([100 150 -.1 .1])
    %hold on;
    drawnow;
end
if par.t==par.tstop;
    xm=s.x(:,2);
    zsm=mean(s.zst,2);
    Hrms_hi=sqrt(mean(s.Ht.^2,2));
    Hrms_lo=sqrt(8)*(std(s.zst'))';   
    Urms=sqrt(mean(s.urmst.^2,2));
    Urms_lo=(std(s.uet'))';
    guls=3*(mean(s.uet.*s.urmst.^2,2)-mean(s.uet,2).*Urms.^2);
    out(:,1)=xm(1:180);
    out(:,2)=zsm(1:180);
    tekal('write','eta.tek',out);
    out(:,2)=Hrms_hi(1:180);
    tekal('write','Hrms_hi.tek',out);
    out(:,2)=Hrms_lo(1:180);
    tekal('write','Hrms_lo.tek',out);
    out(:,2)=Urms(1:180);
    tekal('write','urms.tek',out);
    out(:,2)=Urms_lo(1:180);
    tekal('write','urms_lo.tek',out);
    out(:,2)=guls(1:180);
    tekal('write','guls.tek',out);
end;
end