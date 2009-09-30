function XBdims=getdimensions
% XBdims=getdimensions
% 
% Must start in XBeach output directory
% Output structure XBdims
% XBdims.nt       = number of regular spatial output timesteps
% XBdims.nx       = number of grid cells in x-direcion
% XBdims.ny       = number of grid cells in y-direcion
% XBdims.ntheta   = number of theta bins in wave module
% XBdims.kmax     = number of sigma layers in Quasi-3D model
% XBdims.ngd      = number of sediment classes
% XBdims.nd       = number of sediment class layers
% XBdims.ntp      = number of point output timesteps
% XBdims.ntc      = number of transect output timesteps
% XBdims.ntm      = number of time-average output timesteps
% XBdims.version  = Revision number of XBeach program
% XBdims.tsglobal = times at which regular output is given
% XBdims.tspoints = times at which point output is given
% XBdims.tsmean   = times at which time-averaged output is given
% XBdims.x        = world x-coordinates grid
% XBdims.y        = world y-coordinates grid
% 

% Created 05-06-2009 : XBeach-group Delft

fid=fopen('dims.dat','r');

% version test
temp=fread(fid,[11],'double'); 
frewind(fid);
if temp(11)>0
    version='old';
else
    version='eid';
end

switch version
    case ('old');
    XBdims.nt=fread(fid,[1],'double'); 
    XBdims.nx=fread(fid,[1],'double');         
    XBdims.ny=fread(fid,[1],'double'); 
    XBdims.ngd=fread(fid,[1],'double');       
    XBdims.nd=fread(fid,[1],'double');        
    XBdims.ntp=fread(fid,[1],'double');
    XBdims.ntm=fread(fid,[1],'double');  
    XBdims.tsglobal=fread(fid,[XBdims.nt],'double');
    XBdims.tspoints=fread(fid,[XBdims.ntp],'double');
    XBdims.tsmean=fread(fid,[XBdims.ntm],'double');
    fclose(fid);       
    case('eid')        
    XBdims.nt=fread(fid,[1],'double'); 
    XBdims.nx=fread(fid,[1],'double');         
    XBdims.ny=fread(fid,[1],'double');  
    XBdims.ntheta=fread(fid,[1],'double');  
    XBdims.kmax=fread(fid,[1],'double'); 
    XBdims.ngd=fread(fid,[1],'double');       
    XBdims.nd=fread(fid,[1],'double');        
    XBdims.ntp=fread(fid,[1],'double');
    XBdims.ntc=fread(fid,[1],'double');
    XBdims.ntm=fread(fid,[1],'double');  
    XBdims.version=-1*fread(fid,[1],'double');
    XBdims.tsglobal=fread(fid,[XBdims.nt],'double');
    XBdims.tspoints=fread(fid,[XBdims.ntp],'double');
    XBdims.tscross=fread(fid,[XBdims.ntc],'double');
    XBdims.tsmean=fread(fid,[XBdims.ntm],'double');
    fclose(fid);
end

fidxy=fopen('xy.dat','r');
XBdims.x=fread(fidxy,[XBdims.nx+1,XBdims.ny+1],'double');   
XBdims.y=fread(fidxy,[XBdims.nx+1,XBdims.ny+1],'double');   
fclose(fidxy);