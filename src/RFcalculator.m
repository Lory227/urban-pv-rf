% Radiative forcing calculation based on albedo change
% This code reads the GSA albedo results, and calculates the monthly
% radiative forcing against the albedo change

% author: Yilong Zhou
% Year: 2024

% References:
% Bright & O'Halloran (2019), doi: 10.5194/gmd-12-3975-2019
%  Wohlfahrt, Tomelleri & Hammerle, doi: https://doi.org/10.5281/zenodo.4432576

clear; 

Currentdir = pwd;
lastidx = max(strfind(Currentdir,'\'));
HOMEdir = Currentdir(1:lastidx);
INPUTdir = [HOMEdir,'INPUTS\']; addpath(INPUTdir);
addpath(HOMEdir+"Supplementary functions");

%% Load CACK 
% Code retrieved from https://doi.org/10.5281/zenodo.4432576
time_start = tic;
Currentdir = pwd;
filename_cack = [Currentdir,'\CACKv1.0\CACKv1.0.nc']; % CACK filename
disp('Loading CACK ...')
info = ncinfo(filename_cack);
lon_cack = ncread(info.Filename,'Longitude');
lat_cack = (ncread(info.Filename,'Latitude')); % inverted
% mo_cack = (ncread(info.Filename,'Month'));
yy_cack = (ncread(info.Filename,'Year'));
yy_cack = yy_cack + 2000; 
cack = ncread(info.Filename,'CACK');
cack_sig = ncread(info.Filename,'Sigma_total');
disp('Done loading CACK')

%% Determine the longitude and latitude of pixels, and retrieve the kernal of the location

% **************** The tile of DSM simulated *****************
DSM_name = "C_37EN1";

% Load the DSM
[Tiff.A,Tiff.info] = readgeoraster(INPUTdir+"DSMs\"+DSM_name+".tif");          % load the digital surface model
bbox = [Tiff.info.XWorldLimits', Tiff.info.YWorldLimits'];                % bouding box
H = Tiff.A; info = Tiff.info;
[~,~,X,Y] = mapinpolygon(H,bbox,bbox(:,1),bbox(:,2),1,true);

% Create a 10 by 13 grid
x_idx = 250:500:5000; y_idx = 250:500:6250; 
x_grid = X(1,x_idx); y_grid = Y(y_idx,1);
[x,y] = meshgrid(x_grid,y_grid);
[Latitude,Longitude] = projinv(info.ProjectedCRS,x,y);

clearvars DSM_name HOMEdir Tiff H x_idx y_idx X Y

% Extract kernel for each pixel in the DSM
cellSize = size(Latitude,1);
Latitude = Latitude(:); Longitude = Longitude(:);
cack_p_1 = cell(cellSize,1); cack_p_sig_1 = cell(cellSize,1);

for i = 1:size(Latitude,1)
    [cack_p_1{i}, cack_p_sig_1{i}] = fct_CACK_extract(lat_cack, lon_cack, ...
        yy_cack, cack, cack_sig, ...
        Latitude(i), Longitude(i), ...
        2001, 2016);
end
cack_p_1 = reshape(cack_p_1,13,10);
cack_p_sig_1 = reshape(cack_p_sig_1,13,10);
Latitude1 = Latitude; Longitude1 = Longitude;

% ************** The tile of DSM simulated *****************
DSM_name = "C_37EN2";

% Load the DSM
[Tiff.A,Tiff.info] = readgeoraster(INPUTdir+"DSMs\"+DSM_name+".tif");          % load the digital surface model
bbox = [Tiff.info.XWorldLimits', Tiff.info.YWorldLimits'];                % bouding box
H = Tiff.A; info = Tiff.info;
[~,~,X,Y] = mapinpolygon(H,bbox,bbox(:,1),bbox(:,2),1,true);

% Create a 10 by 13 grid
x_idx = 250:500:5000; y_idx = 250:500:6250; 
x_grid = X(1,x_idx); y_grid = Y(y_idx,1);
[x,y] = meshgrid(x_grid,y_grid);
[Latitude,Longitude] = projinv(info.ProjectedCRS,x,y);

clearvars DSM_name HOMEdir Tiff H x_idx y_idx X Y

% Extract kernel for each pixel in the DSM
cellSize = size(Latitude,1);
Latitude = Latitude(:); Longitude = Longitude(:);
cack_p_2 = cell(cellSize,1); cack_p_sig_2 = cell(cellSize,1);

for i = 1:size(Latitude,1)
    [cack_p_2{i}, cack_p_sig_2{i}] = fct_CACK_extract(lat_cack, lon_cack, ...
        yy_cack, cack, cack_sig, ...
        Latitude(i), Longitude(i), ...
        2001, 2016);
end
cack_p_2 = reshape(cack_p_2,13,10);                             % monthly average climatological kernel at chosen point
cack_p_sig_2 = reshape(cack_p_sig_2,13,10);              % uncertainty of monthly average climatological kernel at chosen point
Latitude2 = Latitude; Longitude2 = Longitude;

% ************** The tile of DSM simulated *****************
DSM_name = "C_37EZ1";

% Load the DSM
[Tiff.A,Tiff.info] = readgeoraster(INPUTdir+"DSMs\"+DSM_name+".tif");          % load the digital surface model
bbox = [Tiff.info.XWorldLimits', Tiff.info.YWorldLimits'];                % bouding box
H = Tiff.A; info = Tiff.info;
[~,~,X,Y] = mapinpolygon(H,bbox,bbox(:,1),bbox(:,2),1,true);

% Create a 10 by 3 grid
x_idx = 250:500:5000; y_idx = 500:500:1500; 
x_grid = X(1,x_idx); y_grid = Y(y_idx,1);
[x,y] = meshgrid(x_grid,y_grid);
[Latitude,Longitude] = projinv(info.ProjectedCRS,x,y);
projectedCRS = info.ProjectedCRS;

clearvars DSM_name HOMEdir Tiff H x_idx y_idx X Y

% Extract kernel for each pixel in the DSM
cellSize = size(Latitude,1);
Latitude = Latitude(:); Longitude = Longitude(:);
cack_p_3 = cell(cellSize,1); cack_p_sig_3 = cell(cellSize,1);

for i = 1:size(Latitude,1)
    [cack_p_3{i}, cack_p_sig_3{i}] = fct_CACK_extract(lat_cack, lon_cack, ...
        yy_cack, cack, cack_sig, ...
        Latitude(i), Longitude(i), ...
        2001, 2016);
end
cack_p_3 = reshape(cack_p_3,3,10);                             % monthly average climatological kernel at chosen point
cack_p_sig_3 = reshape(cack_p_sig_3,3,10);              % uncertainty of monthly average climatological kernel at chosen point
Latitude3 = Latitude; Longitude3 = Longitude;

% ************** The tile of DSM simulated *****************
DSM_name = "C_37EZ2";

% Load the DSM
[Tiff.A,Tiff.info] = readgeoraster(INPUTdir+"DSMs\"+DSM_name+".tif");          % load the digital surface model
bbox = [Tiff.info.XWorldLimits', Tiff.info.YWorldLimits'];                % bouding box
H = Tiff.A; info = Tiff.info;
[~,~,X,Y] = mapinpolygon(H,bbox,bbox(:,1),bbox(:,2),1,true);

% Create a 10 by 3 grid
x_idx = 250:500:5000; y_idx = 500:500:1500; 
x_grid = X(1,x_idx); y_grid = Y(y_idx,1);
[x,y] = meshgrid(x_grid,y_grid);
[Latitude,Longitude] = projinv(projectedCRS,x,y);

clearvars DSM_name HOMEdir Tiff H x_idx y_idx X Y

% Extract kernel for each pixel in the DSM
cellSize = size(Latitude,1);
Latitude = Latitude(:); Longitude = Longitude(:);
cack_p_4 = cell(cellSize,1); cack_p_sig_4 = cell(cellSize,1);

for i = 1:size(Latitude,1)
    [cack_p_4{i}, cack_p_sig_4{i}] = fct_CACK_extract(lat_cack, lon_cack, ...
        yy_cack, cack, cack_sig, ...
        Latitude(i), Longitude(i), ...
        2001, 2016);
end
cack_p_4 = reshape(cack_p_4,3,10);                             % monthly average climatological kernel at chosen point
cack_p_sig_4 = reshape(cack_p_sig_4,3,10);              % uncertainty of monthly average climatological kernel at chosen point
Latitude4 = Latitude; Longitude4 = Longitude;

Latitude = [reshape(Latitude1,13,10),reshape(Latitude2,13,10);reshape(Latitude3,3,10),reshape(Latitude4,3,10)];
Longitude= [reshape(Longitude1,13,10),reshape(Longitude2,13,10);reshape(Longitude3,3,10),reshape(Longitude4,3,10)];
cack_p = [cack_p_1, cack_p_2; cack_p_3, cack_p_4];
cack_p_sig = [cack_p_sig_1, cack_p_sig_2; cack_p_sig_3, cack_p_sig_4];
clear Latitude1 Latitude2 Latitude3 Latitude4 Longitude1 Longitude2 Longitude3 Longitude4

%% Load the albedo and PV results
load(INPUTdir+"albedoResults.mat");
load(INPUTdir+"PVresults.mat");

% The roof tilt angles are extracted for surface area calculation
Tilts = arrayfun(@(x) arrayfun(@(y) y.Tilt, x.Planes, 'UniformOutput', false), BP, 'UniformOutput', false);
Yield = arrayfun(@(x) x.Yield, BP,  'UniformOutput', false); 
Yield_filt = arrayfun(@(x) x.Yield_filt, BP, 'UniformOutput',false);  

Ppv = 0;    % PV potential       (kWh)
Ppv_filt = 0;      % PV potential (>650kWh/kWp)
Apv = 0;     % surface area of roof with integrated PV 
Apv_filt = sum(cell2mat(arrayfun(@(x) x.Area_filt,BP,'UniformOutput',false)))/1e6;    % surface area of roof with integrated PV (>650kWh/kWp)
Ae = 510.1e6;               % surface area of the earth (km2)

% The DSM used for PV calculation has cell size of 0.5 by 0.5 m2
for i = 1:length(BP)
    if ~isempty(Yield{i})
        currentY = Yield{i};
        for j = 1:length(currentY)
            currentpv = double(sum(currentY{j}) * 0.25 / cosd(Tilts{i}{j}));    % 0.25 is the cell area
            Ppv = Ppv + currentpv;
            Apv = Apv + 0.25 * length(currentY{j}) / cosd(Tilts{i}{j});
        end
    end
    if ~isempty(Yield_filt{i})
        currentY_filt = Yield_filt{i};
        for j = 1:length(currentY_filt)
            if ~isempty(currentY_filt{j})
                currentpv_filt = sum(currentY_filt{j}) * 0.25 / cosd(Tilts{i}{j});
                Ppv_filt = Ppv_filt + currentpv_filt;
            end
        end
    end
end
Apv = Apv / 1e6;                                                % km2

clear Yield Yield_filt BP currentY currentY_filt currentpv currentpv_filt

%% Calculate the radiative forcing

% ******************** Potivie RF *******************

rf_alb_unfilteredPV = 0; rf_alb_filteredPV = 0; 
rf_alb_sig_unfilteredPV = 0; rf_alb_sig_filteredPV = 0;
pixelsize = 0.25;

for i = 1:length(albedo_monthly(:))
    rf_alb_unfilteredPV = rf_alb_unfilteredPV + mean((albedo_monthly{i} - albedo_monthly_unfilteredPV{i})' / 100 .* cack_p{i}) * pixelsize/Ae;                                    
    rf_alb_filteredPV = rf_alb_filteredPV + mean((albedo_monthly{i} - albedo_monthly_filteredPV{i})' / 100 .* cack_p{i}) * pixelsize/Ae;                                    
    rf_alb_sig_unfilteredPV = rf_alb_sig_unfilteredPV + mean((albedo_monthly{i} - albedo_monthly_unfilteredPV{i})' /100.* cack_p_sig{i}) * pixelsize/Ae;
    rf_alb_sig_filteredPV = rf_alb_sig_filteredPV + mean((albedo_monthly{i} - albedo_monthly_filteredPV{i})' /100 .* cack_p_sig{i}) * pixelsize/Ae;
end

% *************** Negative RF ****************
% Carbon emission intensity for the Netherlands (kg CO2 kWh-1)
% \https://www.oecd-ilibrary.org/energy/co2-emissions-from-fuel-combustion_22199446\
y_intensity = [0.421, 0.446, 0.445, 0.472, 0.493, 0.464, 0.437];            % 2010 & 2012 - 2017
kwh_co2_NL = mean(y_intensity);
%kwh_co2_NL = 0.013;                  % Case for highly decarbonized energy grid: 0.013 to 0.03 

% Following based on fct_RFco2_unc.mat in Wohlfahrt et al. (2021)
% logspace time
ts = -2; te = 3; tn = 1000;     % log-spcaced time time vector for BET integration
t = [0 logspace(ts, te, tn)]; % log-spaced time vector (y)
dt = diff(t); % time difference vector (y)

% Constants 
a_irf = [0.2173, 0.2240, 0.2824, 0.2763]; % Table 5 in Joos et al. (2013)
tau_irf = [394.4, 36.54, 4.304]; % Table 5 in Joos et al. (2013)
kco2 = 1.76 * 10^-15; % radiative efficiency of CO2 (W m-2 kg-1) (Bright et al., 2016)

% Sanity test
% rf_alb = 150*0.1*100/Ae;
% rf_alb_sig = 0; kwh_co2_NL = 0.4; Ppv = 1000*1e6;

% random vectors 
%rv_rf = normrnd(0, rf_alb_sig, [1, tn]);

% impulse response function for CO2
irf = a_irf(1) + a_irf(2) .* exp(-t./tau_irf(1)) + a_irf(3) .* exp(-t./tau_irf(2)) + a_irf(4) .* exp(-t./tau_irf(3));
irf(end) = [];

% instantaneous RF from CO2 source/sink; Eq. 11 from Bright et al. (2016)
sss_unfilteredPV = ones(numel(t)-1,1) .* Ppv .* kwh_co2_NL ; % source/sink strength (kg CO2/y)         
sss_filteredPV = ones(numel(t)-1,1) .* Ppv_filt .* kwh_co2_NL ; % source/sink strength (kg CO2/y)          
RF_co2_unfilteredPV = cumsum(sss_unfilteredPV .* irf' .* dt') .* kco2; % convolve source/sink strength with IRF and integrate over time
RF_co2_filteredPV = cumsum(sss_filteredPV .* irf' .* dt') .* kco2; % convolve source/sink strength with IRF and integrate over time

[row, col] = find(RF_co2_unfilteredPV > repmat(rf_alb_unfilteredPV + rf_alb_sig_unfilteredPV,numel(t)-1,1)); % find occurences of RF_co2 > RF_alb
idx = grpstats(row, col, 'min'); % row index in each column - minimum is the first occurence of RF_co2 > RF_alb
bet_n = numel(idx); % number of columns (= random draws) where a minimum was found - should equal n
intm = t(idx); % break-even times
bet_avg_unfilteredPV = nanmean(intm); % average

[row, col] = find(RF_co2_filteredPV > repmat(rf_alb_filteredPV + rf_alb_sig_filteredPV,numel(t)-1,1)); % find occurences of RF_co2 > RF_alb
idx = grpstats(row, col, 'min'); % row index in each column - minimum is the first occurence of RF_co2 > RF_alb
bet_n = numel(idx); % number of columns (= random draws) where a minimum was found - should equal n
intm = t(idx); % break-even times
bet_avg_filteredPV = nanmean(intm); % average

% Negative RF @ 5 ,10, 15, 20, 25, 30 years (541,601,637,662,681,697)
nRF_y_unfilt= RF_co2_unfilteredPV([541,601,637,662,681,697]);
nRF_y_filt = RF_co2_filteredPV([541,601,637,662,681,697]);

display("The positive radiative forcing for scenario one is found to be " + num2str(rf_alb_unfilteredPV) + "W/m^2");
display("The breakeven time for scenario one is calcualted to be " + num2str(bet_avg_unfilteredPV*365) + "days");
display("The positive radiative forcing for scenario two is found to be " + num2str(rf_alb_filteredPV) + "W/m^2");
display("The breakeven time for scenario two is calcualted to be " + num2str(bet_avg_filteredPV*365) + "days");
display("The PV integrated roof area for scenario one is " + num2str(Apv) + "m^2, and the annual PV energy yield is " + num2str(Ppv/1e6) + "GWh/year." );
display("The PV integrated roof area for scenario two is " + num2str(Apv_filt) + "m^2, and the annual PV energy yield is " + num2str(Ppv_filt/1e6) + "GWh/year." );


time_end = toc(time_start);

