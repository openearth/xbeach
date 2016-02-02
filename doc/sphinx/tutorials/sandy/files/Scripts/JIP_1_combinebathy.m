%% XBeach JIP tutorial: how to set-up a storm-impact model?
% v0.1 Nov-15
clear all, close all, clc

%% 0. Define
destin = 'e:\Projecten\TKI\Data\';

% Determine area of interest (not boundaries from XBeach grid, but close)
output_lon_w        = 579100; 
output_lon_e        = 589100;  
output_lat_n        = 4439200;  
output_lat_s        = 4428850;  

% Determine size of grid cells
output_delta_x      = 2;
output_delta_y      = 2;
[output_X,output_Y] = meshgrid(output_lon_w:output_delta_x:output_lon_e, output_lat_n:-output_delta_y:output_lat_s);

%% 1. Load CRM
cd(destin);
[input_Z_CRM R_CRM]            = arcgridread_v2('CRM.asc'); % 

% Get grid information from input dataset     
[input_n_y input_n_x]           = size(input_Z_CRM);
input_delta_x                   = R_CRM(2,1);
input_delta_y                   = R_CRM(1,2);
input_lon_w                     = R_CRM(3,1);
input_lon_e                     = R_CRM(3,1)+input_delta_x*(input_n_x-1);
input_lat_n                     = R_CRM(3,2);
input_lat_s                     = R_CRM(3,2)+input_delta_y*(input_n_y-1);

% Convert coordaintes
[input_X_CRM, input_Y_CRM]          = meshgrid(input_lon_w:R_CRM(2,1):input_lon_e,input_lat_n:R_CRM(1,2):input_lat_s);
[input_X_CRM_x,input_Y_CRM_y,logs]  = convertCoordinates(input_X_CRM,input_Y_CRM,'CS1.name','WGS 84','CS1.type','geo','CS2.name','WGS 84 / UTM zone 18N','CS2.type','xy')

% Griddata to output grid
output_Z_CRM                     = griddata(input_X_CRM_x,input_Y_CRM_y,input_Z_CRM, output_X, output_Y);

%% 2. Load LIDAR pre
cd(destin);
[input_Z_pre R_pre]            = arcgridread_v2('2012_USGS_preSandy_NJ.asc'); % 

% Get grid information from input dataset     
[input_n_y input_n_x]           = size(input_Z_pre);
input_delta_x                   = R_pre(2,1);
input_delta_y                   = R_pre(1,2);
input_lon_w                     = R_pre(3,1);
input_lon_e                     = R_pre(3,1)+input_delta_x*(input_n_x-1);
input_lat_n                     = R_pre(3,2);
input_lat_s                     = R_pre(3,2)+input_delta_y*(input_n_y-1);

% Convert coordaintes
[input_X_pre, input_Y_pre]      = meshgrid(input_lon_w:R_pre(2,1):input_lon_e,input_lat_n:R_pre(1,2):input_lat_s);

% Interpolate to grid
output_Z_pre                     = interp2(input_X_pre,input_Y_pre,input_Z_pre, output_X, output_Y);

%% 3. Load LIDAR post
cd(destin);
[input_Z_post R_post]            = arcgridread_v2('2012_USGS_postSandy_NJ.asc'); % 

% Get grid information from input dataset     
[input_n_y input_n_x]           = size(input_Z_post);
input_delta_x                   = R_post(2,1);
input_delta_y                   = R_post(1,2);
input_lon_w                     = R_post(3,1);
input_lon_e                     = R_post(3,1)+input_delta_x*(input_n_x-1);
input_lat_n                     = R_post(3,2);
input_lat_s                     = R_post(3,2)+input_delta_y*(input_n_y-1);

% Convert coordaintes
[input_X_post, input_Y_post]      = meshgrid(input_lon_w:R_post(2,1):input_lon_e,input_lat_n:R_post(1,2):input_lat_s);

% Interpolate to grid
output_Z_post                   = interp2(input_X_post,input_Y_post,input_Z_post, output_X, output_Y);

%% 4. Combine data. Apply CRM where we don't have LiDAR
idlidar = (isnan(output_Z_pre) | isnan(output_Z_post));
output_Z_pre(idlidar) = output_Z_CRM(idlidar);
output_Z_post(idlidar) = output_Z_CRM(idlidar);

%% 5. Smoothing or other things
%

%% 6. Result and save
figure; pcolor(output_X, output_Y, output_Z_pre); shading flat; caxis([-20 10])
figure; pcolor(output_X, output_Y, output_Z_post); shading flat; caxis([-20 10])
figure; pcolor(output_X, output_Y, output_Z_post-output_Z_pre); shading flat; caxis([-1 1]); colormap(erosed)