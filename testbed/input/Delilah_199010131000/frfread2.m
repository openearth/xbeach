% read FRF-8m array spectra 
% and interpolate to hourly values for longshore
% current computations
%
%
clear all
close all
ic = 0
for jd = 13:13
   for jh = [16]
      ic = ic + 1
 filein = 'LA9010131600.ASC';
 if jd < 10
  filein(8:8) = num2str(jd);
else
  filein(7:8) = num2str(jd);
end
if jh < 10
  filein(10:10) = num2str(jh);
else
  filein(9:10) = num2str(jh);
end

 fid=fopen(filein);

dn  = 10;
nk  = 29;   % number of frequency bins
ni  = 12;
nj  = 8;

for l = 1:14
   line = fgetl(fid);
   if isempty(line), break, end
   if ~isstr(line), break, end

   if l == 1
      Hm0_frf(ic) = str2num(line(11:20));
      fp_frf(ic)  = str2num(line(21:30));
      th_frf(ic)  = str2num(line(31:40));
   end
   if l == 3
      wl_frf(ic) = str2num(line(21:30));
      ws_frf(ic) = str2num(line(51:60));
   end
   if l == 4
      wd_frf(ic) = str2num(line(1:10));
   end
   
end
for k = 1:nk
   line = fgetl(fid);
   if isempty(line), break, end
   if ~isstr(line), break, end
   f(k)= str2num(line(11:20));
   A(k)= str2num(line(25:40));
   for i = 1:ni
       line = fgetl(fid);
       if isempty(line), break, end
       if ~isstr(line), break, end
       if i == ni
        for j = 1:3
         s(k,j+(i-1)*nj)   = str2num(line(1+(j-1)*dn:j*dn))*A(k);           
        end
       else
        for j = 1:nj
         s(k,j+(i-1)*nj)   = str2num(line(1+(j-1)*dn:j*dn))*A(k);           
        end
       end
     end
  end
                        
fclose(fid);

% plot results
theta = [-90:2:90];
mesh(theta,f,s)
xlabel('\theta^{o}')
ylabel('f (Hz)')
pause(1)

theta = theta*pi/180;
S_d     = s*180/pi;

% get mean characteristics
 
df = f(2)-f(1);
dthet= theta(2)-theta(1);
w = 2*pi*f;
% select short waves only
f_lo = 0.05;
f_hi = 0.3;
[kf] = find(f<f_hi & f>f_lo);
% only incident  wave directions
[jm] = find(theta<pi/2 & theta> -pi/2);
% mean wave height
Hm0(ic)   = 4*sqrt(sum(sum(S_d(kf,jm)))*dthet*df);
% mean wave direction
theta_0(ic) = sum(S_d(kf,jm)*theta(jm)')/sum(sum(S_d(kf,jm)));
% peak period
[Sp,jp] = max(sum(S_d(kf,jm)'));
wp = 2*pi*f(kf(jp));
fp(ic) = wp/2/pi;
tm(ic) = jd + jh/24 + ((2+16/60)*.5)/24;  % time at the middle of the data collection 

end
end

