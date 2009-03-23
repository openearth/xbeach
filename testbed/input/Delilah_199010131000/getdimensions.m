function [nt,x,y]=getdimensions

fid=fopen('dims.dat','r');
nt=fread(fid,[1],'double');
nx=fread(fid,[1],'double');
ny=fread(fid,[1],'double');
fclose(fid);

fidxy=fopen('xy.dat','r');
x=fread(fidxy,[nx+1,ny+1],'double');
y=fread(fidxy,[nx+1,ny+1],'double');
fclose(fidxy);