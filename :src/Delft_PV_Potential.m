% This script calculates the PV potential of the building roofs in Delft
% based on simplified skyline based model. 

% References: 1) https://doi.org/10.1002/solr.202100478
%                       2) https://doi.org/10.1016/j.rser.2023.113885
%                       3) https://doi.org/10.1038/s41560-018-0318-6
%                       

% Author: Yilong Zhou
% Date: March, 2024

clear; clc;

Currentdir = pwd;
lastidx = max(strfind(Currentdir,'\'));
HOMEdir = Currentdir(1:lastidx); 
INPUTdir = [HOMEdir,'INPUTS\']; 
addpath(INPUTdir);
addpath(HOMEdir+"Supplementary functions"); 
addpath(Currentdir+"\Coefficients");

%% Load the DSM of Delft
INPUTdir = [Currentdir,'\Inputs\'];
fileList = dir([INPUTdir,'*.tif']);

keyword = 'Gebied';                    
idx = find(contains({fileList.name},keyword));

[H,info] = readgeoraster(fullfile(INPUTdir,fileList(idx).name));
bbox = [info.XWorldLimits',info.YWorldLimits'];
H(H>10000) = NaN; % Remove noise (water reflect) points
% ptCloud = heightmaptopointcloud(H,bbox,0.5);      % Display the point cloud
% pcshow(ptCloud.Location);                                   % Display the point cloud

% Initiate solar map for imagesc
solarmap = zeros(size(H)); solarmap_filt = zeros(size(H));

%% Load the shapefiles for buildings in Delft
[S,A] = shaperead([INPUTdir,'pand_delft_gebied.shp']);
logicalIndex = arrayfun(@(x) x.aantal_ver~=0, A);
BP= S(logicalIndex); A = A(logicalIndex);
clear S logicalIndex

%% Crop out the DSM context for each building footprint (+/- raduis) for roof plane detection
% Initializing the parameters
radius = 10;
H_Context = cell(size(BP,1),1);
bbox_Context = cell(size(BP,1),1);

% Specify the number of workers you want to use
numworkers = 32; 

% Create a new parallel pool with the specified number of workers
parpool('local', numworkers);

% Run through each building footprint
parfor n = 1:length(BP)
    [H_Context{n}, ~, ~, bbox_Context{n}] = buildingHmap(BP(n).BoundingBox, H, bbox, radius, 0.5);
end
delete(gcp('nocreate'));
pause(10);

% ptcloud = heightmaptopointcloud(H_Context{24000},bbox_Context{24000},0.5);
% pcshow(ptcloud.Location);

%% Loop through each building for roof plane dection

simTime = struct('planeExtraction',[],'solarCal',[]);

% Start a stopwatch
t_initStart = tic;

% Initialize parameters
if ~exist('MinArea','var');         MinArea = 10;               end

% Specify the number of workers you want to use
numworkers = 56; 

% Create a new parallel pool with the specified number of workers
parpool('local', numworkers);

% Define the ParforProgMon object
ppm = ParforProgMon('Roof plane detection: ', size(BP,1));

% ********* START PLANE EXTRACTION LOOP *********

% Set steps to run the loop for.k
steps = 1:size(BP,1); errorInfo = cell(size(steps));

parfor n = steps
    warning('off','all');
    rGnd = [];rTop = []; nrSect = []; rType = []; rPts = [];
    rAzim = []; rTilt = []; rNormal = []; rArea = []; rBnd = [];
    
    % Retrieve the x,y and bounding box of the footprint
    x = cell2mat({BP(n).X}'); y = cell2mat({BP(n).Y}');
    x = x(1:end-1); y = y(1:end-1);

    % Reduce points grid to points inside building bbox
    [H_b,bbox_b,~,~] = mapinpolygon(H_Context{n},bbox_Context{n},x,y,0.5);
    
    try
    % Extract all building faces (roof and walls) larger than minimum area
    [rGnd,rTop,nrSect,rType,rAzim,rTilt,rNormal,rArea,rBnd,rPts] = roofplanedetection(H_b,bbox_b,[x',y'],MinArea,1.5,0.5);
    
    catch ME
        errorInfo{n} = sprintf('Error at n = %d: %s', n, ME.message);    
    end

    % WRITE TO RESULTS
    planeStruct = struct('Type',rType, 'Azimuth',rAzim, 'Tilt',rTilt,...
        'Area',rArea,'fN',rNormal, 'Boundary',rBnd, 'Points',rPts);
    BP(n).Planes = planeStruct;  BP(n).Area = sum(cell2mat(rArea));
    BP(n).Top = rTop; BP(n).Ground = rGnd;

    %clear H_b bbox_b rType rAzim rTilt rArea rNormal rBnd rPts planeStruct rGnd rTop
    ppm.increment()
end

delete(gcp('nocreate'));

% Find out which buildings are not plane detected 
count = 0; restpes = [];
hold on                  
for i = 1:length(errorInfo)
    if ~isempty(errorInfo{i})
        %disp(errorInfo{i}); 
        count = count + 1;
        resteps = [restpes,i];
        %mapshow(BP(i).X,BP(i).Y,'DisplayType','Polygon');
    end
end

if ~isempty(resteps)

% Specify the number of workers you want to use
numworkers = 56; 

% Create a new parallel pool with the specified number of workers
parpool('local', numworkers);

% ********* START PLANE EXTRACTION LOOP *********

% Rerun the plane extraction for the undetected surfaces in the previous step
errorInfo2 = cell(size(resteps));

parfor n = resteps
    warning('off','all');
    rGnd = [];rTop = []; nrSect = []; rType = []; rPts = [];
    rAzim = []; rTilt = []; rNormal = []; rArea = []; rBnd = [];
    
    % Retrieve the x,y and bounding box of the footprint
    x = cell2mat({BP(n).X}'); y = cell2mat({BP(n).Y}');
    x = x(1:end-1); y = y(1:end-1);

    % Reduce points grid to points inside building bbox
    [H_b,bbox_b,~,~] = mapinpolygon(H_Context{n},bbox_Context{n},x,y,0.5);  % bbox_b(:,1),bbox_b(:,2)
    
    try
    % Extract all building faces (roof and walls) larger than minimum area
    [rGnd,rTop,nrSect,rType,rAzim,rTilt,rNormal,rArea,rBnd,rPts] = roofplanedetection(H_b,bbox_b,[x',y'],MinArea,1.3,0.5);
    
    catch ME
        errorInfo2{n} = sprintf('Error at n = %d: %s', n, ME.message);    
        %errorInfo{n} = n; % Simple version
    end

    % WRITE TO RESULTS
    planeStruct = struct('Type',rType, 'Azimuth',rAzim, 'Tilt',rTilt,...
        'Area',rArea,'fN',rNormal, 'Boundary',rBnd, 'Points',rPts);
    BP(n).Planes = planeStruct;  BP(n).Area = sum(cell2mat(rArea));
    BP(n).Top = rTop; BP(n).Ground = rGnd;

    %clear H_b bbox_b rType rAzim rTilt rArea rNormal rBnd rPts planeStruct rGnd rTop
end

delete(gcp('nocreate'));

end

time_end = toc(t_initStart);
simTime.planeExtraction = time_end;

%% Crop out the DSM context for each building footprint (+/- raduis) for horizon scan
% Initializing the parameters
radius = 100;
H_Context = cell(size(BP,1),1); bbox_Context = cell(size(BP,1),1);

% Specify the number of workers you want to use
numworkers = 32; 

% Create a new parallel pool with the specified number of workers
parpool('local', numworkers); 

% Run through each building footprint
parfor n = 1:length(BP)
    [H_Context{n}, ~, ~, bbox_Context{n}] = buildingHmap(BP(n).BoundingBox, H, bbox, radius, 0.5);
end
delete(gcp('nocreate'));
pause(10);

% ptcloud = heightmaptopointcloud(H_Context{24000},Context_bbox{24000},0.5);
% pcshow(ptcloud.Location);

%% Calculate the solar PV potential

% Start the stopwatch
t_start = tic;

% ******* Initialize parameters *******
REGIONstring = 'Delft';
altN = 90;                        % Number of steps of sky dome in 90 degrees
azimN = 360;                  % Number of steps of sky dome in 360 degrees
horiN = 360/5;               % Number of steps of horizon scanner azimuth
timeInt = 30;                  % minutes between sunhour points

% *******  PV module characteristics *******
M.R = 0.1;                                    % module reflectivity index
M.eps_top = 0.2;                        % module top emissivity
M.eps_back = 0.89;                    % module back emissivity
M.tracking_type = 'fixed';
M.albedo_gnd = 0.2;

M.efficiency =  0.215; % STC overall module efficiency
M.deff_dT = -0.0029;    % deta/dT [%/Celcius]
NOCT=45;                % [Celcius]
M.INOCT = NOCT + 18;    % TNOCT+18(direct mount) in Celsius
M.L = 1.558;            % module length [m] (vertical)
M.W = 1.046;            % module width [m] horizontal
M.Voc = 68.2;          % module Voc at STC
M.Ns = 96;              % number of cells in series
M.n = 1.2;              % ideality factor (tech dependent)
M.Pmax = 0.350;                       % Peak power [kWp]

% Some initial cells
mapidx = cell(size(BP,1),1);
mapfiltidx = cell(size(BP,1),1);
dim = size(solarmap);
errorInfo2 = cell(size(steps));
S = cell(size(BP,1),1);

% ******* Get Airmass correction matrix *******
height_km = nanmean(H,'all')/1000;
f_AM =  AMfactorMat(height_km,altN,azimN);

% ******* Loop through each building surface *******

% Specify the number of workers you want to use
numworkers = 48; % Change this number to the desired number of workers

% Create a new parallel pool with the specified number of workers
parpool('local', numworkers);


parfor n = steps

    if isempty(BP(n).Planes(1).Points); continue; end
    currentBP = BP(n).Planes;
    currentContext = H_Context{n}; currentbbox = bbox_Context{n};
    buildingfaces = length(currentBP);

    % ******* Find maximum solar dome time factor for clear sky *******
    [LON,LAT] = rd2wgs(nanmean(currentbbox(:,1)),nanmean(currentbbox(:,2)));   
    coord = [LAT,LON];                                                      % Delft coordinates = [52.01,4.357];

    % Get solar time [h] for each azimuth and elevation angle of the sun
    sunMat = genSunPositionMat(coord,altN,azimN,timeInt);                   % annalemma of the sun
    f_aoi_max = maxAOIfactor(sunMat);

    % Initialize some empty items
    bd_idx = cell(buildingfaces,1); bdfilt_idx = cell(buildingfaces,1);
    Area_filt = 0; boundarys = cell(buildingfaces,1);

    for r = 1 : buildingfaces

        % Get roof orientation and points
        face_tilt = currentBP(r).Tilt;                               % Tilt of surface mesh
        face_azim = currentBP(r).Azimuth;                            % Orientation of surface mesh
        face_pts_sim = currentBP(r).Points + currentBP(r).fN * 0.5;  % Points on surface mesh lifted by 50cm along normal
        face_pts = currentBP(r).Points;
        face_Npts = size(currentBP(r).Points,1);                     % Number of points on surface mesh

        % Get the horizon elevation for each azimuth angle of all face points
        try
            
            shapes = Horizonscanner(currentContext, currentbbox, face_pts_sim, face_azim, face_tilt, radius, azimN);
        

        % Sky sector view matrix following A.Calcabrini
        [sun_if,sky_vf] = genSkyDomeMat(sunMat,f_AM,f_aoi_max,face_azim,face_tilt);

        % Get skyline matrix for all points on the building face (SCF and SVF calculation)
        SCF = ones(face_Npts,1); SVF = zeros(face_Npts,1);
        irr = zeros(face_Npts,1); yld = zeros(face_Npts,1);
        for p = 1 : face_Npts
            skyline_raw = shapes{p};
            [SCF(p),SVF(p)] = SCFnSVFCalc_Revised_Fast(skyline_raw, sun_if, sky_vf);
        end

        % Irradiance and Yield Calculation
        [irr,yld] = IrrYieldCalc_Revised(REGIONstring,M,face_azim,face_tilt,SVF,SCF,1,1);
        Yld = yld/(M.L*M.W);                                  % Yield results

        % Apply the specific yield threshold
        yld_filt = yld / M.Pmax;
        yld_idx = find(yld_filt > 650);
        Y_filt = double(yld(yld_idx)/(M.L*M.W));
        P_filt = face_pts(yld_idx,:);

        % Find the boundary of filtered points
        [B_idx,~] = boundary(P_filt(:,1),P_filt(:,2),0.5);
        polyX = P_filt(B_idx,1); polyY = P_filt(B_idx,2);

        % Add to results (Entire building surfaces)
        BP(n).Irradiation{r,1} = irr;
        BP(n).Yield{r,1} = Yld; 
        boundarys{r,1} = [polyX,polyY];
        if ~isempty(yld_idx)
            BP(n).Yield_filt{r,1} = Y_filt;
            BP(n).Points_filt{r,1} = P_filt;
        end
        
        colIndices = round((currentBP(r).Points(:,1) - bbox(1,1) - 0.25) / 0.5);
        rowIndices = round((bbox(2,2) - 0.25 - currentBP(r).Points(:,2)) / 0.5);
        bd_idx{r} = sub2ind(dim, rowIndices, colIndices);
        %solarmap(sub2ind(size(solarmap), rowIndices, colIndices)) = Y{r};

        % Replace for the filtered ones
        colIndices = round((P_filt(:,1) - bbox(1,1) - 0.25) / 0.5);
        rowIndices = round((bbox(2,2) - 0.25 - P_filt(:,2)) / 0.5);
        bdfilt_idx{r} = sub2ind(dim, rowIndices, colIndices);
        %solarmap_filt(sub2ind(size(solarmap_filt), rowIndices, colIndices)) = Y_filt{r};
        Area_filt = Area_filt + size(Y_filt,1)*0.25/cosd(face_tilt);

        catch ME

            errorInfo2{n} = sprintf('Error at n = %d: %s', n, ME.message);

        end
    end
    BP(n).Area_filt = Area_filt;
    mapidx{n} = bd_idx; mapfiltidx{n} = bdfilt_idx;
    S{n} = boundarys;
    
end
delete(gcp('nocreate'));
pause(10);


time_end = toc(t_start);
simTime.solarCal = time_end;

%% Save to the shapefile
% Creating a structure for the shapefile

% for k = 1:length(polygons)
%     S(k).Geometry = 'Polygon';
%     S(k).BoundingBox = [min(polygons{k}(1,:)), min(polygons{k}(2,:)); max(polygons{k}(1,:)), max(polygons{k}(2,:))];
%     S(k).X = polygons{k}(1,:);
%     S(k).Y = polygons{k}(2,:);
% end


%% Assign the yield to each DSM pixel
for i = steps
    yieldidx = cell2mat(mapidx{i});
    yieldnum = cell2mat(BP(i).Yield);
    if isempty(yieldidx); i 
    else
        solarmap(yieldidx) = yieldnum;

    end
end

for i = steps
    yieldidx = cell2mat(mapfiltidx{i});
    try
    yieldnum = cell2mat(BP(i).Yield_filt);
    if isempty(yieldidx); i
    else
    solarmap_filt(yieldidx) = yieldnum;
    end
    catch ME
        infomap{i} = sprintf('Error at i = %d: %s', i, ME.message);
    end
    
end

%% Create the structure for roof plane polygons
p = 0;
buildings_poly_filt = struct('Geometry', {}, 'BoundingBox', {}, 'X', {}, 'Y', {});
for i = steps
    ksteps = length(S{i});
    for k = 1:ksteps
        if ~isempty(S{i}{k})
            buildings_poly_filt(p+1).Geometry = 'Polygon';
            buildings_poly_filt(p+1).BoundingBox = [min(S{i}{k}(:,1)), min(S{i}{k}(:,2)); max(S{i}{k}(:,1)), max(S{i}{k}(:,2))];
            buildings_poly_filt(p+1).X = S{i}{k}(:,1);
            buildings_poly_filt(p+1).Y = S{i}{k}(:,2);
            p = p + 1;
        end
    end
end

% hold on
% for p = 1:length(buildings_poly_filt)
%     mapshow(buildings_poly_filt(p).X,buildings_poly_filt(p).Y,'DisplayType','polygon');
% end

p = 0;
buildings_poly = struct('Geometry', {}, 'BoundingBox', {}, 'X', {}, 'Y', {});
for i = steps
    ksteps = length(BP(i).Planes);
    for k = 1:ksteps
        if ~isempty(BP(i).Planes(k).Boundary)
            bdry = BP(i).Planes(k).Boundary;
            buildings_poly(p+1).Geometry = 'Polygon';
            buildings_poly(p+1).BoundingBox = [min(bdry(:,1)), min(bdry(:,2)); max(bdry(:,1)), max(bdry(:,2))];
            buildings_poly(p+1).X = bdry(:,1);
            buildings_poly(p+1).Y = bdry(:,2);
            p = p + 1;
        end
    end
end

% hold on
% for p = 1:length(buildings_poly)
%     mapshow(buildings_poly(p).X,buildings_poly(p).Y,'DisplayType','polygon');
% end

