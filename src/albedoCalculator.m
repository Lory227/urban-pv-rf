%% This fucntional tool calculates local albedo of the given area using LiDAR data
%  The underlying albedo calcualtion algorithm is based on the geometric spectral albedo mdoel
%  developed by H.Ziar et.al. (2019) https://doi.org/10.1016/j.apenergy.2019.113867.
%
%  Inputs:        DSM from LiDAR data
%                       Weather data for the location of interest (Meteonorm & SMARTS)               
%                       Spectral reflectances of surface materials (ASTER library)       
%                       Ecosystem Unit Map
%
%
%  Outputs:     F_A1:                          View factor for illuminated and visible segments 
%                       F_A2:                          View facotr for shaded and visible segements
%                       H_albedo:                  H value for albedo caculation
%                       Reflectivity:               Reflectivity of the area of interest
%                       TIME:                          Simulation time for each step 
%                       albedo:                       Hourly albedo results 
%                       albedo_weighted:    Irradiance-weighted albedo
%                       binary_DSM:             Material map of the area of interest 
%                       roughness:                 Roughness of the area of interest 
%                       shadingProfile:          Visibility map from the albedometer over the area
%
%      This project uses ParforProgMon to monitor the progress of parallel computing. 
%      It was originally developed by Dylan Muir, Willem-Jan de Goeij, and The MathWorks, Inc., under the BSD 3-Clause License.
%      See license.txt under ParforMon folder for details.
%       
%      ** The simualted albedo results can be found in the folder INPUTS\albedoResults.mat


%  Author: Yilong Zhou 
%  Date: January, 2025

clear; clc;

Currentdir = pwd;
lastidx = max(strfind(Currentdir,'\'));
HOMEdir = Currentdir(1:lastidx);
INPUTdir = [HOMEdir,'INPUTS']; addpath(INPUTdir);
addpath(HOMEdir+"Supplementary functions");

%% Initialization
% Select the DSM tile to be simulated
list = {'C_37EN1', 'C_37EN2', 'C_37EZ1', 'C_37EZ2'};       % Four DSM tiles containing Delft
[indx,tf] = listdlg('ListString',list,'SelectionMode','single',...
    'PromptString','Select one DSM tile','ListSize', [200 150]);      
DSM_name = list{indx};               % DSM tile to be simulated
if str2double(DSM_name(end)) == 1; Tile_name = "Delft1_9T"; else; Tile_name = "Delft_9T"; end

% Determine if roofs are integrated with PV
choice = questdlg('Do you want to put PV on roofs?','Confirmation', 'Yes', 'No', 'No');
if strcmp(choice,'No'); boolean = 0; else; boolean = 1; end

% Select one scenario to simulate the roof PV integration
% Scenario 1:       Roof PV is filtered with an economic threshold of 650 kWh/kWp (filt)
% Scenario 2:       Roof PV is unfiltered (nonfilt)
if boolean == 1
    choice = questdlg('Do you want to filter roof PV with 650 kWh/kWp threshold?','Confirmation', 'Yes', 'No', 'No');
    if strcmp(choice,'No') 
        pvpolygon = 'buildings_poly.shp'; 
    else            
        pvpolygon = 'buildings_poly_filt.shp'; 
    end
end

%% Load the DSM

[Tiff.A,Tiff.info] = readgeoraster(string(INPUTdir)+"\DSMs\"+DSM_name+".tif");          % Load the digital surface model
bbox = [Tiff.info.XWorldLimits', Tiff.info.YWorldLimits'];                % Find the bouding box (bbox)
H = Tiff.A; info = Tiff.info;                                                                     % Extract the DSM height data and info
[~,~,X,Y] = mapinpolygon(H,bbox,bbox(:,1),bbox(:,2),1,true);      % Map the height data based on the bbox

% Gnerate a grid of points for albedo simulation
% C_37EN_1, C_37EN2:    grid size of 10 by 13
% C_37EZ1, C_37EZ2:        grid size of 10 by 3
x_idx = 250:500:5000; 
if strcmp(DSM_name,"C_37EN1") || strcmp(DSM_name,"C_37EN2") == 1
    y_idx = 250:500:6250;       
    gridsteps = 1:13;           % 10 by 13 grid
else
    y_idx = 500:500:6250;       
    gridsteps = 1:3;             % 10 by 3 grid

end
x_grid = X(1,x_idx); y_grid = Y(y_idx,1);

clearvars -except x_grid y_grid HOMEdir Currentdir INPUTdir Tile_name a b gridsteps boolean pvpolygon

%% Loop through each albedometer position

for a = 1:10
    for b = gridsteps
        
        TIME = struct();
        
        % Crop the DSM given the location of the albedometer
        [Tiff.A,Tiff.info] = readgeoraster(string(INPUTdir)+'\DSMs\'+Tile_name+'.tif');          % Load the digital surface model of 9 tiles
        bbox_H = [Tiff.info.XWorldLimits', Tiff.info.YWorldLimits'];                            % Bouding box
        H = Tiff.A; info = Tiff.info;
        clear Tiff
        albedoMeter = [x_grid(a),y_grid(b)];
        bbox_grid = [albedoMeter(1)-2500,albedoMeter(2)-3125;albedoMeter(1)+2501,albedoMeter(2)+3126];
        [H,bbox,X,Y] = mapinpolygon(H,bbox_H,bbox_grid(:,1),bbox_grid(:,2),1,true);             % The DSM for the albedometer location is extracted

        %% Calculate the view factor of the extracted DSM

        time_start = tic;

        % Define the height of the albedometer (500m above the mean height of the DSM)
        %  doi: 10.5067/DOC/CEOSWGCV/LPV/ALBEDO.001 
        Ha = 500 + mean(H(:));                      

        % Calculate the viewing angle and the view factor of each pixel
        [vf,theta_v,proj_albedo_x,proj_albedo_y] = calculate_view_factors(size(H,2),size(H,1),1,Ha,H);
        VF = sum(vf(:));

        if VF < 0.95
            warning('The view factor is too low! Place your albedometer lower.');
        end

        % The view factor calculation time
        TIME.ViewFactor = toc(time_start);

        %% Load the TIFF files for buildings and vegetation derived from original LiDAR data

        time_start = tic;

        % Load the DSM for buildings of the tile
        [binary_building,~] = readgeoraster(string(INPUTdir)+'\DSMs\'+ Tile_name + "_buildings.tif");
        binary_building(binary_building~=0) = 1;
        [binary_building,~,~,~] = mapinpolygon(binary_building,bbox_H,bbox_grid(:,1),bbox_grid(:,2),1,true);

        % Load the DSM for vegetations
        [binary_vegetation,~] = readgeoraster(string(INPUTdir)+'\DSMs\'+ Tile_name +"_vegetation.tif");
        binary_vegetation(binary_vegetation~=0) = 1;
        [binary_vegetation,~,~,~] = mapinpolygon(binary_vegetation,bbox_H,bbox_grid(:,1),bbox_grid(:,2),1,true);

        % TheTIFF seperation time
        TIME.TIFFseperation = toc(time_start);

        %% Further classify the DSM into subgroups (preprocessed from Ecosystem Unit Map)
        % The Ecosystem Unit Map was first processed in ArcGIS Pro, and
        % each material class was exported as an individual shapefile
        % https://www.cbs.nl/en-gb/background/2017/12/ecosystem-unit-map

        % Water:            Water (code:2)
        % Ground:          Sand/Meadow/GreyAlsphat&Roads/Grassland (code:3/4/5/6)
        % Vegetation:    Coniferous Forest/Deciduous Forest/Mixed Forest (code:7/8/9)
        % Building:         Residential/Greenhouses/Industrial (code:10/11/12)

        time_start = tic;

        % Initiate a binary matrix for the DSM material
        binary_DSM = single(zeros(size(H,1),size(H,2)));

        % *********Classify the DSM into Water subgroup*********
        Water = shaperead(string(INPUTdir)+'\Shapefiles\Water.shp','BoundingBox',bbox);

        % Loop through the shapes in the shapefile
        binary_water = single(zeros(size(H,1),size(H,2)));

        % Shapefile of Water (code 2)
        Mask_Water = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Water)
            for i = 1:size(Water,1)
                % Get the current shape
                shape = Water(i);
                % Update the mask
                mask_water = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Water = Mask_Water | mask_water;
            end
            Mask_Water = Mask_Water .* 2;
            binary_water = binary_water + Mask_Water;
            binary_DSM = binary_water;
        end

        clear Water Mask_Water mask_water

        % *********Classify the DSM into Ground subgroup*********
        Sand = shaperead(string(INPUTdir)+'\Shapefiles\Ground_Sand.shp','BoundingBox',bbox);
        Meadow = shaperead(string(INPUTdir)+'\Shapefiles\Ground_Meadow.shp','BoundingBox',bbox);
        GreyAlsphat = shaperead(string(INPUTdir)+'\Shapefiles\Ground_GreyAlsphat.shp','BoundingBox',bbox);
        Roads = shaperead(string(INPUTdir)+'\Shapefiles\Ground_Roads.shp','BoundingBox',bbox);
        Grassland = shaperead(string(INPUTdir)+'\Shapefiles\Ground_Grassland.shp','BoundingBox',bbox);

        % Loop through the shapes in the shapefiles
        binary_ground = single(zeros(size(H,1),size(H,2)));

        % Shapefile of Sand (code 3)
        Mask_Ground = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Sand)
            for i = 1 : size(Sand,1)
                % Get the current shape
                shape = Sand(i);
                % Update the mask
                mask_sand = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Ground = Mask_Ground | mask_sand;
            end
            binary_ground(Mask_Ground==1) = 3;                   % Assign the material code number
            binary_DSM(Mask_Ground==1) = 3;
        end
        clear Sand mask_sand

        % Shapefile of Meadow (code 4)
        Mask_Ground = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Meadow)
            for i = 1 : size(Meadow,1)
                % Get the current shape
                shape = Meadow(i);
                % Update the mask
                mask_meadow = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Ground = Mask_Ground | mask_meadow;
            end
            binary_ground(Mask_Ground==1) = 4;                   % Assign the material code number
            binary_DSM(Mask_Ground==1) = 4;
        end
        clear Meadow mask_meadow

        % Shapefile of Alsphat (code 5)
        Mask_Ground = single(zeros(size(H,1),size(H,2)));
        if ~isempty(GreyAlsphat)
            for i = 1 : size(GreyAlsphat,1)
                % Get the current shape
                shape = GreyAlsphat(i);
                % Update the mask
                mask_greyalsphat = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Ground = Mask_Ground | mask_greyalsphat;
            end
            binary_ground(Mask_Ground==1) = 5;                   % Assign the material code number
            binary_DSM(Mask_Ground==1) = 5;
        end
        clear GreyAlsphat mask_greyalsphat

        % Shapefile of Roads (code 5)
        Mask_Ground = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Roads)
            for i = 1 : size(Roads,1)
                % Get the current shape
                shape = Roads(i);
                % Update the mask
                mask_roads = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Ground = Mask_Ground | mask_roads;
            end
            binary_ground(Mask_Ground==1) = 5;                   % Assign the material code number
            binary_DSM(Mask_Ground==1) = 5;
        end
        clear Roads mask_roads

        % Shapefile of Grassland (code 6)
        Mask_Ground = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Grassland)
            for i = 1 : size(Grassland,1)
                % Get the current shape
                shape = Grassland(i);
                % Update the mask
                mask_grassland = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Ground = Mask_Ground | mask_grassland;
            end
            binary_ground(Mask_Ground==1) = 6;                   % Assign the material code number
            binary_DSM(Mask_Ground==1) = 6;
        end
        clear Grassland mask_grassland Mask_Ground

        % *********Classify the DSM into Vegetation subgroup*********
        Coniferous = shaperead(string(INPUTdir)+'\Shapefiles\Vegetation_Coniferous_Forest.shp','BoundingBox',bbox);
        Deciduous = shaperead(string(INPUTdir)+'\Shapefiles\Vegetation_Deciduous_Forest.shp','BoundingBox',bbox);
        Mixed = shaperead(string(INPUTdir)+'\Shapefiles\Vegetation_MixedForest.shp','BoundingBox',bbox);

        % Shapefile of Coniferous (code 7)
        Mask_Coniferous = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Coniferous)
            for i = 1 : size(Coniferous,1)
                % Get the current shape
                shape = Coniferous(i);
                % Update the mask
                mask_coniferous = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Coniferous = Mask_Coniferous | mask_coniferous;
            end
            binary_coniferous = Mask_Coniferous & binary_vegetation;
            binary_vegetation(binary_coniferous) = 7;
            binary_DSM(binary_coniferous) = 7;
        end
        clear binary_coniferous mask_coniferous Mask_Coniferous Coniferous

        % Shapefile of Deciduous (code 8)
        Mask_Deciduous = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Deciduous)
            for i = 1 : size(Deciduous,1)
                % Get the current shape
                shape = Deciduous(i);
                % Update the mask
                mask_deciduous = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Deciduous = Mask_Deciduous | mask_deciduous;
            end
            binary_deciduous = Mask_Deciduous & binary_vegetation;
            binary_vegetation(binary_deciduous) = 8;
            binary_DSM(binary_deciduous) = 8;
        end
        clear binary_deciduous mask_deciduous Mask_Deciduous Deciduous

        % Shapefile of Mixed (code 9)
        Mask_Mixed = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Mixed)
            for i = 1 : size(Mixed,1)
                % Get the current shape
                shape = Mixed(i);
                % Update the mask
                mask_mixed = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Mixed = Mask_Mixed | mask_mixed;
            end
            binary_mixed = Mask_Mixed & binary_vegetation;
            binary_vegetation(binary_mixed) = 9;
            binary_DSM(binary_mixed) = 9;
        end

        % Assign the remaining vegetation points to mixed forest cluster
        binary_DSM(binary_vegetation==1) = 9;
        clear binary_mixed mask_mixed Mask_Mixed Mixed

        % *********Classify the DSM into Building subgroup*********
        Greenhouse = shaperead(string(INPUTdir)+'\Shapefiles\Building_Greenhouse.shp','BoundingBox',bbox);
        Residential = shaperead(string(INPUTdir)+'\Shapefiles\Residential_Area.shp','BoundingBox',bbox);

        % Shapefile of Residential Buildings (code 10)
        Mask_Residential = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Residential)
            for i = 1 : size(Residential,1)
                % Get the current shape
                shape = Residential(i);
                % Update the mask
                mask_residential = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Residential = Mask_Residential | mask_residential;
            end
            binary_residential = Mask_Residential & binary_building;
            binary_building(binary_residential) = 10;
            binary_DSM(binary_residential) = 10;
        end
        clear binary_residential mask_residential Mask_Residential Residential

        % Shapefile of Greenhouse (code 11)
        Mask_Greenhouse = single(zeros(size(H,1),size(H,2)));
        if ~isempty(Greenhouse)
            for i = 1 : size(Greenhouse,1)
                % Get the current shape
                shape = Greenhouse(i);
                % Update the mask
                mask_greenhouse = maskgenerator(bbox,shape.X,shape.Y,1);
                Mask_Greenhouse = Mask_Greenhouse | mask_greenhouse;
            end
            binary_greenhouse = Mask_Greenhouse & binary_building;
            binary_building(binary_greenhouse) = 11;
            binary_DSM(binary_greenhouse) = 11;
        end
        clear binary_greenhouse mask_greenhouse Mask_Greenhouse Greenhouse

        % For the rest of the buildings, assume they are industrial buildings (code 12)
        binary_DSM(binary_building==1) = 12;
        binary_building(binary_building==1) = 12;

        % For pixels not assigned to any subgroup, classify them as ground.
        % Determine the dominant ground code in this region and assign it to these pixels.
        countmax = 0;
        for i = 3:6
            count = sum(ismember(binary_DSM, i),'all');  % Count occurrences of ground subgroups
            if count > countmax
                countmax = count;
                codeidx = i;
            end
        end
        binary_DSM(binary_DSM==0) = codeidx;
        clear countmax count codeidx binary_ground binary_building binary_vegetation binary_water

        % The subgroup classification time
        TIME.SubgroupClassification = toc(time_start);
        pause(10);

        %% Shading Profile Calculation 
        % This part computes if the DSM cell is visible to the albedometer with LOS assessment
        % For each albedometer position, this shading profile is unique and
        % only needs to be computed once.

        time_start = tic;

        % Converting coordinates in geographical format
        info.XWorldLimits = [bbox(1,1),bbox(2,1)];
        info.YWorldLimits = [bbox(1,2),bbox(2,2)];
        info.RasterSize = [size(H,1),size(H,2)];
        [YLimits,XLimits] = projinv(info.ProjectedCRS,info.XWorldLimits,info.YWorldLimits);
        R = georefcells(YLimits,XLimits,size(H),'ColumnsStartFrom','north');

        % Converting world coordinates to instrinsic, then to geographic coordiantees
        [SurIntrinsicX,SurIntrinsicY] = worldToIntrinsic(info,albedoMeter(1),albedoMeter(2));
        [LatSur,LonSur] = intrinsicToGeographic(R,SurIntrinsicX,SurIntrinsicY);
        [ObsIntrinsicX,ObsIntrinsicY] = worldToIntrinsic(info,X,Y);
        [LatObs,LonObs] = intrinsicToGeographic(R,ObsIntrinsicX,ObsIntrinsicY);

        % Line-of-sight (LOS) assessment with los2
        LonSur = ones(size(X,1),size(X,2)) * single(LonSur);
        LatSur = ones(size(X,1),size(X,2)) * single(LatSur);
        shadingProfile = single(zeros(size(X,1),size(X,2)));

        % Specify the number of workers you want to use
        numworkers = 60; % Change this number to the desired number of workers (60)

        % Create a new parallel pool with the specified number of workers
        parpool('local', numworkers);

        % Define the ParforProgMon object
        ppm = ParforProgMon('Processing: ', size(X,1));

        parfor j = 1:size(X,1)
            shadingProfile(j,:) = los2(H,R,LatSur(j,:),LonSur(j,:),LatObs(j,:),LonObs(j,:),Ha);
            ppm.increment();                % increment the progress bar
        end

        % Plot the shading profile
        % figure(2)
        % [r, c] = size(shadingProfile);
        % imagesc((1:c)+0.5, (1:r)+0.5, shadingProfile);
        % colormap(gray);
        % close all

        TIME.ShadingProfile = toc(time_start);
        clear LatObs LonObs ObsIntrinsicY ObsIntrinsicX ans
        delete(gcp('nocreate'));

        %% Scale up/down the SMARTS data based on Meteonorm data, and calculate the effective view factor in overcase conditions

        % Load the meteorological data
        load('Brightness4.mat');
        sun_alt = single(sun_alt); sun_azim = single(sun_azim);

        % Calculate the roughness of the DSM
        [roughness,~] = DSM_roughness(H); 

        % Calculate the effective view factor in overcast sky condition
        clearH = ~overcastH;            % Indices for clear sky conditions
        F_overcast = vf .* shadingProfile; F_overcast_sum = sum(F_overcast,'all');
        
        % Retrieve the spectrum wavelength range
        wavelength = readmatrix("astmg173.xls",'Sheet','SMARTS2','Range','A3:A2004');
        wavelength = single(wavelength(241:1992));              % From 400nm to 3950nm given the wavelength range of materials

        % Rigid scale of spectrum data based on Meteonorm and SMARTS GHIs
        GHI_spectrum = single(GHI_spectrum(241:1992,:));
        GHI_tot = trapz(wavelength,GHI_spectrum)';          % GHI of SMARTS data
        GHI_norm = GHI_Meteo./GHI_tot;                            % Calculate the weight between meteonorm GHI and SMARTS GHI
        GHI_norm = repmat(GHI_norm',[1752,1]);              % Reshape the array
        GHI_spectrum = GHI_spectrum .* GHI_norm;        % Scale the SMARTS GHI rigidly based on Meteonorm data

        clear GHI_norm P_vis

        % Load the factor H (H = DNI/DHI*cosd(theta)), which is
        % precalcualted with Meteonorm data while generating the Brightness4.mat
        H_albedo = single(H_albedo);
        H_albedo(isnan(H_albedo)) = 0;

        %% Load the spectral reflectance of materials and calculate the hourly reflectivity

        load([INPUTdir,'\Materials.mat']);           % All preprossed spectral reflectance data from ASTER library
        GHI_tot = trapz(wavelength,GHI_spectrum);
        
        % Reflectivity of all the materials
        Ref.Water = trapz(wavelength,GHI_spectrum .* Materials.Water(:,2))./GHI_tot;
        Ref.Tiles = trapz(wavelength,GHI_spectrum .* Materials.Tiles(:,2))./GHI_tot;
        Ref.Sand = trapz(wavelength,GHI_spectrum .* Materials.Sand(:,2))./GHI_tot;
        Ref.RoofAlsphat = trapz(wavelength,GHI_spectrum .* Materials.RoofAlsphat(:,2))./GHI_tot;
        Ref.MixedForest{1} = trapz(wavelength,GHI_spectrum .* Materials.MixedForest{1}(:,2))./GHI_tot;
        Ref.MixedForest{2} = trapz(wavelength,GHI_spectrum .* Materials.MixedForest{2}(:,2))./GHI_tot;
        Ref.MixedForest{3} = trapz(wavelength,GHI_spectrum .* Materials.MixedForest{3}(:,2))./GHI_tot;
        Ref.MixedForest{4} = trapz(wavelength,GHI_spectrum .* Materials.MixedForest{4}(:,2))./GHI_tot;
        Ref.Meadow = trapz(wavelength,GHI_spectrum .* Materials.Meadow(:,2))./GHI_tot;
        Ref.GreyAlsphat = trapz(wavelength,GHI_spectrum .* Materials.GreyAlsphat(:,2))./GHI_tot;
        Ref.Greenhouse = trapz(wavelength,GHI_spectrum .* Materials.Greenhouse(:,2))./GHI_tot;
        Ref.Grass= trapz(wavelength,GHI_spectrum .* Materials.Grass(:,2))./GHI_tot;
        Ref.DeciduousForest{1} = trapz(wavelength,GHI_spectrum .* Materials.DeciduousForest{1}(:,2))./GHI_tot;
        Ref.DeciduousForest{2} = trapz(wavelength,GHI_spectrum .* Materials.DeciduousForest{2}(:,2))./GHI_tot;
        Ref.DeciduousForest{3} = trapz(wavelength,GHI_spectrum .* Materials.DeciduousForest{3}(:,2))./GHI_tot;
        Ref.DeciduousForest{4} = trapz(wavelength,GHI_spectrum .* Materials.DeciduousForest{4}(:,2))./GHI_tot;
        Ref.ConiferousForest{1} = trapz(wavelength,GHI_spectrum .* Materials.ConiferousForest{1}(:,2))./GHI_tot;
        Ref.ConiferousForest{2} = trapz(wavelength,GHI_spectrum .* Materials.ConiferousForest{2}(:,2))./GHI_tot;
        Ref.ConiferousForest{3} = trapz(wavelength,GHI_spectrum .* Materials.ConiferousForest{3}(:,2))./GHI_tot;
        Ref.ConiferousForest{4} = trapz(wavelength,GHI_spectrum .* Materials.ConiferousForest{4}(:,2))./GHI_tot;

        % Create the reflectivity arrays for Spring, Summer, Autumn and Winter
        % Each row represents one material and the row number corresponds to the
        % material index defined in the previous section

        % Spring
        R_spring = [Ref.Water; Ref.Sand; Ref.Meadow; Ref.GreyAlsphat; Ref.Grass; Ref.ConiferousForest{1};...
            Ref.DeciduousForest{1}; Ref.MixedForest{1}; Ref.Tiles; Ref.Greenhouse; Ref.GreyAlsphat];

        % Summer
        R_summer = [Ref.Water; Ref.Sand; Ref.Meadow; Ref.GreyAlsphat; Ref.Grass; Ref.ConiferousForest{2};...
            Ref.DeciduousForest{2}; Ref.MixedForest{2}; Ref.Tiles; Ref.Greenhouse; Ref.GreyAlsphat];

        % Fall
        R_fall = [Ref.Water; Ref.Sand; Ref.Meadow; Ref.GreyAlsphat; Ref.Grass; Ref.ConiferousForest{3};...
            Ref.DeciduousForest{3}; Ref.MixedForest{3}; Ref.Tiles; Ref.Greenhouse; Ref.GreyAlsphat];

        % Winter
        R_winter = [Ref.Water; Ref.Sand; Ref.Meadow; Ref.GreyAlsphat; Ref.Grass; Ref.ConiferousForest{4};...
            Ref.DeciduousForest{4}; Ref.MixedForest{4}; Ref.Tiles; Ref.Greenhouse; Ref.GreyAlsphat];
        
        % Reflectivity of all daylight hours for one year
        % Spring:           Marth, April, May; 
        % Summer:       June, July, August;
        % Autumn:        September, October, November; 
        % Winter:          December, January, February
        R_year = single([R_winter(:,1:522), R_spring(:,523:1769), ...
            R_summer(:,1770:3202), R_fall(:,3203:3902), R_winter(:,3903:4387)]);
        
        % Replace the roof material with PV material if boolean is 1
        if boolean == 1
            % Load the PV material
            load('PVlaminate.mat');
            PV_material = trapz(wavelength,GHI_spectrum .* PVlaminate(:,2))./GHI_tot;
            R_year = [R_year;PV_material];
            clear PVlaminate
            
            % Shapefile of PV integrated buildings (code 13)
            warning off
            polyPV = shaperead([INPUTdir,'\Shapefiles\',pvpolygon],'BoundingBox',bbox);     % Roof PV polygon generated from PV potential simulation   
            Mask_PV = single(zeros(size(H,1),size(H,2)));
            if ~isempty(polyPV)
                for i = 1 : length(polyPV)
                    % Get the current shape
                    shape = polyPV(i);
                    % Update the mask
                    mask_PV = maskgenerator(bbox,shape.X,shape.Y,1);
                    Mask_PV = Mask_PV | mask_PV;
                end
                binary_DSM(Mask_PV) = 13;
            end
            clear mask_PV Mask_PV polyPV
        end
       
        clear Materials R_spring R_summer R_fall R_winter PV_material X_1 Y_1 in on range_mask roof_mask combined_mask

        %% Create the hourly Reflectivity map based on the materials
        time_start = tic;

        sizeR = size(R_year,2); Reflectivity = zeros(sizeR,1); % Initiate the array

        % Specify the number of workers
        numWorkers = 48; 

        % Create a new parallel pool with the specified number of workers
        parpool('local', numWorkers);

        % Define the ParforProgMon object
        ppm = ParforProgMon('Hourly Reflectivity: ', sizeR);

        parfor i = 1:sizeR
            R_map = R_year(:,i);
            M_map = binary_DSM;
            for j = 2:size(R_map,1)+1
                M_map(M_map == j) = R_map(j-1,1);
            end
            Reflectivity(i) = sum(M_map .* vf,'all');
            ppm.increment();
        end

        delete(gcp('nocreate'));
        clear M_map
        TIME.ReflectivityMap = toc(time_start);
        
        %% Create bins for altitude and azimuth to reduce the number of simulations
        % For similar altitude and azimuth combinations, they are clustered
        % into the same bin for analysis

        alt_bins = [min(sun_alt):2:max(sun_alt), Inf];                  % The default is with bin size of 2
        azim_bins = [min(sun_azim):5:max(sun_azim), Inf];               % The defalut is with bin size of 10

        [alt_bin_idx, ~] = discretize(sun_alt, alt_bins);
        [azm_bin_idx, ~] = discretize(sun_azim, azim_bins);

        % Calculate mean of each bin and store it in new arrays
        mean_alt_bins = zeros(1, length(alt_bins) - 1);
        mean_azim_bins = zeros(1, length(azim_bins) - 1);

        for k = 1:length(alt_bins) - 1
            subset = sun_alt(alt_bin_idx == k);
            if isempty(subset)
                mean_alt_bins(k) = NaN;  
            else
                mean_alt_bins(k) = mean(subset);
            end
        end

        for k = 1:length(azim_bins) - 1
            subset = sun_azim(azm_bin_idx == k);
            if isempty(subset)
                mean_azim_bins(k) = NaN;  
            else
                mean_azim_bins(k) = mean(subset);
            end
        end

        % Replace the original bins with the mean bins
        alt_bins(1:end-1) = mean_alt_bins;
        azim_bins(1:end-1) = mean_azim_bins;

        % Create unique_ori with the mean bin values (altitude,azimuth)
        [unique_bins,~,ic] = unique([alt_bin_idx, azm_bin_idx], 'rows');
        unique_ori = [alt_bins(unique_bins(:,1))', azim_bins(unique_bins(:,2))'];

        % Remove the index of the class which does not exist
        Material_all = unique(binary_DSM);
        Material_indx = Material_all-1;
        R_year = R_year(Material_indx,:);

        clear alt_bins azim_bins alt_bin_idx azm_bin_idx unique_bins

        %% Initiations for the albedo calculation
        time_start = tic;

        % Initiatae the illumination profile cell
        albedo = zeros(length(sun_alt),1); arrayLength = length(unique_ori);
        F_A1 = zeros(length(sun_alt),1); F_A2 = zeros(length(sun_alt),1);
        
        % Select the chunkSize
        chunkSize = 20;
        binary_DSM = binary_DSM(:);
        numChunks = ceil(arrayLength / chunkSize);
        ori_chunks = splitIntoChunks2D(unique_ori, chunkSize, 1);
        ill_array = cell(arrayLength,1);
        for i = 1:arrayLength
            ill_array{i} = find(ismember(ic,i));
        end

        indexMatrix = zeros(numChunks, 2);
        for i = 1:numChunks
            indexMatrix(i, 1) = (i - 1) * chunkSize + 1;
            indexMatrix(i, 2) = i * chunkSize;
        end
        indexMatrix(end,2) = arrayLength;

        %% Calculate the hourly illumination profile and the albedo
        numWorkers = 20; % Specify the desired number of workers

        % Create a new parallel pool with the specified number of workers
        parpool('local', numWorkers);
        
        for k = 1:numChunks
            orientation = ori_chunks{k};
            if k == numChunks; chunkSize = size(orientation,1); end

            % Initiate the cells
            P_shd = cell(chunkSize,1); P_ill_vis = cell(chunkSize,1);
            illumination_profile = cell(chunkSize,1);

            parfor i = 1:chunkSize
                illumination_profile{i} = illumincationProfile(orientation(i,1),orientation(i,2),X,Y,H);
            end
            parfor i = 1:chunkSize
                [P_ill_vis{i},P_shd{i}] = Pcalculator(theta_v,orientation(i,1),orientation(i,2),proj_albedo_x,proj_albedo_y,roughness);
            end            

            ill_idx = ill_array(indexMatrix(k,1):indexMatrix(k,2));
            albedo_temp_total = cell(chunkSize, 1);
            F_A1_temp = cell(chunkSize, 1);
            F_A2_temp = cell(chunkSize, 1);

            for j = 1:chunkSize
                time_idx = ill_idx{j};
                albedo_temp = zeros(length(time_idx),1);
                F_ill_temp = zeros(length(time_idx),1);
                F_shd_temp = zeros(length(time_idx),1);

                % Calculate the view factor for illuminated and shaded
                % cells which are visible to the albedometer
                F_ill = illumination_profile{j} .* shadingProfile .* vf;            % illuminated and visible
                F_shd = (1 - illumination_profile{j}) .* shadingProfile .* vf;      % shaded and visible
              
                % Effevtive view factors (considering self-shadowing)
                F_shd = F_shd + F_ill .* P_shd{j};
                F_ill = F_ill .* P_ill_vis{j};                

                parfor p = 1:length(time_idx)
                    F_all = clearH(time_idx(p)) .* (F_ill + (1 / (1 + H_albedo(time_idx(p)))) .* F_shd)...
                        + overcastH(time_idx(p)) .* (1 / (1 + H_albedo(time_idx(p)))) .* F_overcast;
                    F_material = accumarray(binary_DSM, F_all(:));
                    F_material = F_material(Material_all);

                    % Calculate the albedo and effective view factors
                    albedo_temp(p) = sum(F_material .* R_year(:,time_idx(p)));
                    F_ill_temp(p) = clearH(time_idx(p)) * sum(F_ill,'all');
                    F_shd_temp(p) = clearH(time_idx(p)) * sum(F_shd,'all') + overcastH(time_idx(p)) * F_overcast_sum;
                end
                albedo_temp_total{j} = albedo_temp;
                F_A1_temp{j} = F_ill_temp;
                F_A2_temp{j} = F_shd_temp;
            end

            for j = 1:chunkSize
                time_idx = ill_idx{j};
                albedo(time_idx)= albedo_temp_total{j};
                F_A1(time_idx) = F_A1_temp{j};
                F_A2(time_idx) = F_A2_temp{j};
            end
            k
        end
        
        % Calculate the GHI weighted albedo
        albedo_weighted = sum(albedo .* GHI_tot') / sum(GHI_tot);

        delete(gcp('nocreate'));
        TIME.albedoCalculation = toc(time_start);
        
        filename = ['albedo',num2str(a),num2str(b),'_500m.mat'];
        save(filename,'albedo','albedo_weighted','F_A1','F_A2','H_albedo','Reflectivity','roughness','TIME','binary_DSM','shadingProfile');
        clearvars -except x_grid y_grid HOMEdir INPUTdir Currentdir Tile_name a b DSM_name boolean pvpolygon gridsteps
        pause(20);
    end
end


%% Functions
function [view_factors,Beta,proj_albedo_x,proj_albedo_y] = calculate_view_factors(grid_X, grid_Y, cell_size, observer_height, surface_height)

% Calculate the center of the grid
center_X = (grid_X + 1) / 2;
center_Y = (grid_Y + 1) / 2;

% Generate the grid of X and Y coordinates
[X, Y] = meshgrid(1:grid_X, 1:grid_Y);

% Calculate the distances from the observer to the centers of the pixels
X_distance = (X - center_X) * cell_size;
Y_distance = (Y - center_Y) * cell_size;
sqrdist = X_distance.^2 + Y_distance.^2;
clear X_distance Y_distance

proj_albedo_x = single(center_X - X);
proj_albedo_y = single(center_Y - Y);

% Adjust observer height with the surface height at each grid point
height_difference = observer_height - surface_height;

% Calculate the distance considering surface height
distance = sqrt(sqrdist + height_difference.^2);

% Calculate the angles Beta1 and Beta2 which are equal to each other
% when the surface normals are in parallel; here Beta is also the
% viewing angle
Beta = atan2d(sqrt(sqrdist), height_difference);

% Calculate the view factor for each pixel (F =
% (cosd(Beta1)*cosd(Beta2)/(pi*d^2)*integral(dA))
view_factors = cosd(Beta).^2./pi./distance.^2.*(cell_size^2);

end



function binary_map = illumincationProfile(Al,Az,X,Y,H_BP)
% This function is based on the method in the following paper:
% http://doi.org/10.1016/j.isprsjprs.2019.06.009

% Grid size set to 1 based on the cell size of the DSM
p = 1;

% Calculate the rotation matrix given the azimuth and altitude of the sun
TransF = [cosd(Az), sind(Az), cosd(Al), sind(Al)];

% Calculate UV and this rotates the DSM along the z axis
U = X.*TransF(1,1) - Y.*TransF(1,2);
V = X.*TransF(1,2) + Y.*TransF(1,1);

% Calculate PQ and put DSM into UQP coordinate system in which the DSM is
% perpendicular to the sun rays. The coordinates in the UQ plane are
% rounded to integers.
P = V.*TransF(1,3) + H_BP.*TransF(1,4);
Q = -V.*TransF(1,4) + H_BP.*TransF(1,3);
U = round(U/p); Q = round(Q/p);

% Adjust P by the minimum value
P = P - min(P(:));

% Reshape the matrix to arrays and sort
uqpts = [U(:), Q(:), P(:)];
[sorted_uqp,sorted_idx] = sortrows(uqpts, 3, 'descend');
[sorted_uqp, sorted_idx2] = sortrows(sorted_uqp, [1,2], 'ascend');
[~,ia,ic] = unique(sorted_uqp(:,1:2),'rows','stable');

% Check which points are on the top with a thickness of 1/sind(Al) meter
tolerance = 1 / sind(Al) ;
uqp_H = sorted_uqp(ia,3);
uqp_highest = [sorted_uqp(ic,1:2), uqp_H(ic)];
height_dif = uqp_highest(:,3) - sorted_uqp(:,3);
idx_H = height_dif<tolerance;

% Generate the shading profile for that hour
binary_map = single(zeros(size(uqpts,1),1));
binary_map(idx_H) = 1;
binary_map = binary_map(invperm(sorted_idx2));
binary_map = binary_map(invperm(sorted_idx));
binary_map = reshape(binary_map,size(H_BP));                        % illumination profile

    function inv = invperm(p)
        [~, inv] = sort(p);
    end
end


function [P_ill_vis,P_shd] = Pcalculator(theta_v,sun_alt,sun_azim,proj_albedo_x,proj_albedo_y,roughness)

theta_i = 90 - sun_alt;

% Create matrix for theta_max and theta_min
theta_max = theta_v; theta_max(theta_v<theta_i) = theta_i;
theta_min = theta_v; theta_min(theta_v>theta_i) = theta_i;

% Calculate k for lamda calculation
proj_solar_vector = [sind(sun_azim), cosd(sun_azim), 0];
Phi = calculate_viewing_angles(proj_solar_vector, proj_albedo_x, proj_albedo_y);
gama = 4.41 * Phi ./ (4.41 * Phi + 1);

% Calculate Lambdas
Lambda_max = LambdaCalculator(theta_max,roughness);
Lambda_min = LambdaCalculator(theta_min,roughness);
Lambda = LambdaCalculator(theta_v,roughness);

% Calculate Probabilities
P_ill_vis = 1 ./ (1 + Lambda_max + gama .* Lambda_min);
%P_vis = 1 ./ (1 + Lambda);
P_shd =  1 ./ (1 + Lambda) - P_ill_vis;



    function Lambda = LambdaCalculator(theta,r)

        cot_theta = cotd(abs(theta));

        part1 = r./(cot_theta.*sqrt(2*pi));
        part2 = exp(-(cot_theta).^2./(2*r^2));
        part3 = 1/2.*erfc(cot_theta./(r*sqrt(2)));

        Lambda = part1 .* part2 - part3;

    end

% Calculate Phi
    function phi = calculate_viewing_angles(proj_solar_vector, proj_to_albedometer_x, proj_to_albedometer_y)
        % Convert proj_solar_vector to a normalized vector
        proj_solar_vector = proj_solar_vector / norm(proj_solar_vector);

        % Calculate the determinant of the two vectors
        determinant = proj_solar_vector(1) * proj_to_albedometer_y - proj_solar_vector(2) * proj_to_albedometer_x;

        % Cosine and sine of the viewing angle
        cos_phi = (proj_solar_vector(1) * proj_to_albedometer_x + proj_solar_vector(2) * proj_to_albedometer_y) / norm(proj_solar_vector);
        sin_phi = determinant / norm(proj_solar_vector);

        % Viewing angle in radians [0, pi]
        phi = atan2(sin_phi, cos_phi);
        phi(phi < 0) = phi(phi < 0) + pi;

    end
end

function [Sdq,theta] = DSM_roughness(Surf)

    % grid density
    dx = 1;
    dy = 1;

    % Calculate differences between adjacent cells
    dSurf_dx = diff(Surf, 1, 2) / dx;  % Differences along columns (x direction)
    dSurf_dy = diff(Surf, 1, 1) / dy;  % Differences along rows (y direction)

    % Calculate mean slopes
    AvgSx = mean(dSurf_dx(:));
    AvgSy = mean(dSurf_dy(:));

    % Calculate variances
    Sx = var(dSurf_dx(:), 1);  % Variance along x
    Sy = var(dSurf_dy(:), 1);  % Variance along y

    % Calculate roughness
    Sdq = sqrt((Sx + Sy) / 2); % Mean of the variances, then square root
    theta = atan(Sdq);
end


function chunks = splitIntoChunks2D(matrix, chunkSize, dim)
numChunks = ceil(size(matrix, dim) / chunkSize);
chunks = cell(numChunks, 1);
for i = 1:numChunks
    startIdx = (i-1) * chunkSize + 1;
    endIdx = min(i * chunkSize, size(matrix, dim));
    if dim == 1
        chunks{i} = matrix(startIdx:endIdx, :);
    else
        chunks{i} = matrix(:, startIdx:endIdx);
    end
end
end

function normals = cellnormal(surface_height,bbox)
% Calculate the normal vector of the DSM cells
ptCloud = heightmaptopointcloud(surface_height, bbox, 1);
normals = pcnormals(ptCloud);

% Flip normals that are pointing downwards
[row, col] = find(normals(:,:,3) < 0);
for i = 1:length(row)
    normals(row(i), col(i), :) = -normals(row(i), col(i), :);
end
end
