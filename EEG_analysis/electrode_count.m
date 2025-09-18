

addpath('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\MultipleTestingToolbox\MultipleTestingToolbox')
addpath('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\Mass_Univariate_ERP_Toolbox-master')
mappeddatadrive='V:';

lengthofbaseline=[500];  %%baseline 500 or 1000 ms prestim
flip=[0]; %%bilater=0; flipX=3;

%%laterality: not flipped
laterality=0;
flipxmuliplier=1;
laterality_label=[''];

%% common baseline
baseline_name='common';
baselinetype='_common_r';

%% baseline start
baseline_start_real=59;
baseline_start=19;
baselinesuffixes='baseline5to0';

%% include baseline in stats
statsinclusion=['_baselinestats_'];
statsinclusionverbose=['including baseline statistics'];
analysis_start=41;

%% otherstuff
baseline_end_real=78;
baseline_end=38;

xtick_spacing_filename = ['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\filter_times_win100ms_shift25ms_full.mat'];
load (xtick_spacing_filename)
final_bin=length(T)-1; %first 500 ms
% runduration=num2str((T(end)-2)*1000+50);
runduration=[num2str(4000) 'ms'];

%% number of parcels
nparcels=[ 80 ];


for overlays=[1:3]
    if isequal(overlays, 1)
        overlay='CP';
        start_pt=44;
        
    elseif isequal(overlays, 2)
        overlay='CNP';
        start_pt=1;
        
    else
        overlay='Subtraction';
        start_pt=1;
        
    end
    
    
    for maps=start_pt:100
        
        
        inflationstep=5;
        createMontage = 0;
        views = 4;
        MapIndex=maps;
        plot_electrodes=[];
        parcel_flavor=['kmeans_mask_refitted']
        
        data_location=[mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies'];
        data_location2=[mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior'];
        save_location_1 = ['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\'];
        save_location_2_1= [num2str(lengthofbaseline) 'ms_' baseline_name '_baseline_' statsinclusion];
        save_location_2_2=['100msbins75msoverlap'];
        save_location_2_3=['filtered_2to2s'];
        save_location_2=[save_location_2_1 save_location_2_3];
        save_location=[save_location_1 save_location_2_2 '\' save_location_2  '\' num2str(nparcels) '_parcels_' parcel_flavor  '\' laterality_label overlay runduration];
        
        
        if ~exist(save_location, 'dir')
            mkdir(save_location)
        end
        
        
        significant_vertices_filename=[save_location '\allverticevalues_' num2str(MapIndex) '.mat'];
        disp(['Running stats for ' laterality_label ' ' overlay ' using map ' num2str(MapIndex) ' with ' num2str(nparcels) ' parcels, with a ' num2str(lengthofbaseline) 'ms ' baseline_name ' baseline '  statsinclusionverbose])
        
        %%takes care of patients
        load([mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_21_04_15.mat']);
        
        
        if strcmp(baselinetype, '_common_r')
            file_suffixes=eval(['mean_suffixes_77bins_100ms_filt_from2to2_' baselinesuffixes baselinetype]);
            %%%these are the suffixes we want, even
            %%%though they are for 77 bins instead of
            %%%157 bins. the 77 bins are a subsample of
            %%%the 157 bins
            
        else
            file_suffixes=eval(['meanpower_suffixes_77bins_100ms_filter_from2to2_' baselinesuffixes baselinetype]);
        end
        patients=patients(1:end-1);
        
        parcellation_folder = ['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\parcellations\parcellation_' num2str(nparcels) '_' parcel_flavor '\store_' num2str(MapIndex) '\'];

       
        
        %% group level
        
        % specify the location to save to
        
        colors = distinguishable_colors(length(patients), .5 * ones(1, 3));
        
        
        VTK_folder_name = 'VTK Leah';
        
        
        
        
        %% Determine if overlay argument was properly passed
        overlay_labels = {'Electrode Distribution', 'Arrival Times', 'CP', 'CNP', 'Electrode Density', 'Subtraction', 'Comparison', 'Parcellation'};
        %{

        %}
        
        % If overlay numeric code, replace with the label
        if isnumeric(overlay)
            overlay = overlay_labels{overlay};
        end
        
        
        switch overlay
            case 'Electrode Distribution'
                image_storage_folder = [save_location '/Electrode_Distribution'];
            case 'CP'
                image_storage_folder = [save_location '/Gamma_power_CP'];
            case 'CNP'
                image_storage_folder = [save_location '/Gamma_power_CNP'];
            case 'Electrode Density'
                image_storage_folder = [save_location '/Electrode_Density'];
            case 'Subtraction'
                image_storage_folder = [save_location '/Gamma_power_Subtraction'];
            case 'Comparison'
                image_storage_folder = [save_location '/Gamma_power_Comparison'];
            case 'Parcellation'
                image_storage_folder = [save_location '/Parcellation'];
        end
        
        
        %% Check more inputs
        valid_keywords = {'IsVisual', 'VertexVals1', 'VertexVals2', 'CompareMethod'};
        
        valid_compare_methods = {'XCorr', 'SignProdSum'};
        
        params.IsVisual = false(1, numel(patients));
        params.VertexVals1 = [];
        params.VertexVals2 = [];
        params.CompareMethod = '';
        
        
        
        
        
        if ~isempty(params.VertexVals1) && isempty(params.VertexVals2)
            throw(MException('displayElectrodesInflated:MissingVertexVals',...
                'Error: VertexVals1 supplied but VertexVals2 missing.'))
        end
        
        if isempty(params.VertexVals1) && ~isempty(params.VertexVals2)
            throw(MException('displayElectrodesInflated:MissingVertexVals',...
                'Error: VertexVals2 supplied but VertexVals1 missing.'))
        end
        
        if ~isempty(params.VertexVals1) && exist(params.VertexVals1, 'file')
            temp_data = load(params.VertexVals1);
            params.VertexVals1 = [];
            params.VertexVals1.vertex_valuesL = temp_data.vertex_valuesL;
            params.VertexVals1.vertex_valuesR = temp_data.vertex_valuesR;
        end
        
        if ~isempty(params.VertexVals2) && exist(params.VertexVals2, 'file')
            temp_data = load(params.VertexVals2);
            params.VertexVals2 = [];
            params.VertexVals2.vertex_valuesL = temp_data.vertex_valuesL;
            params.VertexVals2.vertex_valuesR = temp_data.vertex_valuesR;
        end
        
        if ~all(islogical(params.IsVisual))
            throw(MException('displayElectrodesInflated:IsVisualInvalidType',...
                'Error: ''IsVisual'' keyword argument must be a logical array.'))
        end
        
        if ~isempty(params.CompareMethod) && ~any(strcmp(params.CompareMethod, valid_compare_methods))
            throw(MException('displayElectrodesInflated:InvalidCompareMethod',...
                [sprintf('Error: ''%s''', params.CompareMethod), ' is an invalid comparison method. ',...
                'Valid comparison methods are:', sprintf(' ''%s''', valid_compare_methods{:})]))
        end
        
        
        %% Setup for brainmap creation
        
        
        % change to local version
        mni2fs_dir = 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\mni2fs-master';
        mni_nifti_path = fullfile(mni2fs_dir, 'MNI_T1_1mm_stripped.nii');
        load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
        mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
        Tmni = mnit1_1mm_stripped.transform;
        
        % Define patient electrode colors
        if ~strcmp(overlay, 'Comparison')
            colors = distinguishable_colors(length(patients), .5 * ones(1, 3));
        end
        
        % Define the face transparency value
        alpha = 0.7;
        
        % Side laterality - which side of the brain to draw
        start_side = 1;
        if laterality == 0 % bilateral
            sides_num = 2;
        elseif laterality == 1 % left only
            sides_num = 1;
        elseif laterality == 2 % right only
            start_side = 2;
            sides_num = 2;
        elseif laterality == 3 % do both left and right on the left brain
            sides_num = 1;
        end
        
        % Create figure
        figHandle = gobjects(1, sides_num);
        
                    total_num_frames = 1;

        
        %% Generating all frames for each side of the brain
        for side_index = start_side:sides_num    % Generate frames for L and R
            
            % Determine what side we are working on
            if side_index == 1
                side = 'L';
                hem = 'lh';
            else
                side = 'R';
                hem = 'rh';
            end
         
          
            
            % Initialize empty vertex_values variable
            if overlay ~= 8
                vertex_values=cell(length(patients),1);
                elecVerticesAll=cell(length(patients),1);
                
                % if you are plotting arrival times, then each frame is a different set
                % of electrodes
                if strcmp(overlay, 'Arrival Times')
                    electrode_mapping_times = 10;
                else
                    electrode_mapping_times = 1;
                end
                
           
                    for p = 1 : length(patients)
                        
                        % Display status
                        display(patients{p})
                        
                        % Get electrode info for the patient
                        % this part plots the electrodes for each patient and the overlay data
                        % for each patient
                        patient_folder = [data_location, '/', patients{p}];
                        patient_folder2= [data_location2, '/', patients{p} '/icEEG/HNCT Auditory ID Combined'];
                        % Make the Montage
                        % if the left/right electrode montages have not been created, do so
                        % now, indicated by createMontage = 0 or = 1
                        if createMontage == 1 && ~strcmp(patients{p}, '436BP')
                            create_electrode_montage(patient_folder, laterality == 3, params.IsVisual(p))
                        end
                        
                        % Now load some data files needed for electrode placement and surface shading
                        cd([data_location, '/', patients{p}]);
                        load([side '_MontageMap.mat']);
                        electrode = eval([side '_MontageMap' ]);
                        countelectrode=size(electrode,1);
                        count_tally2=[count_tally2; countelectrode];
                    end
            end
        end
    end
end
                  

%% Below are functions called upon within this function

%==================================================================
%
%Title:         savebrainimages
%Description:   This function takes a brain you have prepared and shines a
%               light on it in 3 views and saves them
%
%Inputs:        "side"          [str] left = 'L', right = 'R'
%               "frame_index"   [int] the "frame" of the "movie" you are
%                               on, which can be #clusters, time period, or
%                               bin number for power analysis
%==================================================================

function frames_set = savebrainimages(side_index,frame_index,views, figHandle)

figure(figHandle(side_index))

if views == 4
    % For each frame, take a snapshot for each of these views
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    % These arrays define the viewing and light perspective of the frame
    zoom_factors = [1.25, 1.25, 1.25, 1.25];
    if side_index == 1
        side = 'L';
        Light_Pos = [1200 -200 500; -700 -800 -100;-200, -100, -800; 200, 800, 0];
        View_Pos = [90 0;-90, 0;0 -90;180 0];
        %         View_Pos_xyz = [200 0 0;-200 0 0;0 0 -200;0 1000 0];
    elseif side_index == 2
        side = 'R';
        Light_Pos = [-1200 -300 700; 1000 -1100 -200; -200, -100, -800;-200, 800, 0];
        View_Pos = [-90 0;90, 0;0 -90;180 0];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif views == 2
    zoom_factors = [1.25, 1.25];
    if side_index == 1
        side = 'L';
        ViewAngles = {'Lateral', 'Medial'};
        Light_Pos = [1200 -200 500; -700 -800 -100];
        View_Pos = [90 0;-90, 0];
        %         Light_Pos = [-700 -800 -300];
        %         View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral', 'Medial'};
        Light_Pos = [-1200 -300 700; 1000 -1100 -200];
        View_Pos = [-90 0;90, 0];
        %         Light_Pos = [-1200 -300 -200];
        %         View_Pos = [-90 -15];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif strcmp(views, 'Lateral')
    zoom_factors = 1.25;
    if side_index == 1
        side = 'L';
        ViewAngles = {'Lateral'};
        Light_Pos = [1200 -200 500];
        View_Pos = [90 0];
        %         Light_Pos = [-700 -800 -300];
        %         View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral'};
        Light_Pos = [-1200 -300 700];
        View_Pos = [-90 0];
        %         Light_Pos = [-1200 -300 -200];
        %         View_Pos = [-90 -15];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif strcmp(views, 'Medial')
    zoom_factors = 1.25;
    if side_index == 1
        side = 'L';
        ViewAngles = {'Medial'};
        Light_Pos = [-700 -800 -100];
        View_Pos = [-90, 0];
        %         Light_Pos = [-700 -800 -300];
        %         View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Medial'};
        Light_Pos = [1000 -1100 -200];
        View_Pos = [90, 0];
        %         Light_Pos = [-1200 -300 -200];
        %         View_Pos = [-90 -15];
    else
        Light_Pos=0;
        View_Pos=0;
    end
end

for i = 1:length(ViewAngles)
    lightHandle = light('Position', Light_Pos(i,:));
    
    view(View_Pos(i,:));
    if i == 1
        % with this command, matlab will freeze the aspect ratio of
        % your image instead of changing the object size to fill the
        % window each time
        axis vis3d
    end
    camzoom(zoom_factors(i))
    drawnow
    fn = [side '_' num2str(frame_index) '_' ViewAngles{i}];
    frames_set(i) = getframe;
    %         frames_set(side_index,frame_index,i) = getframe;
    
    % before you save this picture, add a colorbar
    if frame_index == 1
        colorbar
        print(gcf, fn, '-dtiff');
        %             saveas(gcf, fn);
    end
    
    camzoom(1/zoom_factors(i))
    
    % you have to turn the light off before applying the new light for
    % the next position. The light will be the most recent axes child
    axes_children = get(gca,'children');
    delete(axes_children(1))
    
    % also take off the colorbar so that the next image will not
    % have a colorbar until right before you save it
    if frame_index == 1
        delete(colorbar)
    end
end
end

%==================================================================
%
%Title:         FormatMeshSurface
%Description:   This function loads a surface txt file output from Bioimage
%               Suite into Matlab, and arranges it into a readable format.
%               Repeat twice (once for cortex and once for smooth). The
%               output is a reorganized matrix of vertices.
%
%Inputs:        "vertices_num"  [int] The number of vertices in the file
%               "faces_num"     [int] The number of faces in the file
%               "data_location" [string] should point to the root folder
%                               for the seizure
%                               i.e. 'C:\IcEEGProject_ProcessedData\25\25_1'
%               "filename"      [string] Location of the vertice data
%==================================================================

function [mesh_vertices, mesh_faces]=FormatMeshSurface(vertices_num, faces_num, data_location, filename)
% cd(patient_folder);
% cd ..
cd(data_location)
data_location_previous=pwd;
vertices_range_str=strcat('A1..I',num2str(ceil(vertices_num/3)));
raw_mesh_vertices=dlmread(fullfile( data_location_previous, filename),'',vertices_range_str);
faces_range_str=strcat('B',num2str(ceil(vertices_num/3)+1),'..D',num2str(ceil(vertices_num/3)+faces_num));
raw_mesh_faces=dlmread(fullfile( data_location_previous, filename),'',faces_range_str );
mesh_faces=raw_mesh_faces+1;

vertices_range=1:ceil(vertices_num/3);
mesh_vertices=zeros(ceil(vertices_num/3)*3,3);
mesh_vertices(3*vertices_range-2,1:3)=raw_mesh_vertices(vertices_range, 1:3);
mesh_vertices(3*vertices_range-1,1:3)=raw_mesh_vertices(vertices_range, 4:6);
mesh_vertices(3*vertices_range,1:3)=raw_mesh_vertices(vertices_range, 7:9);
mesh_vertices(ismember(mesh_vertices,zeros(1,3),'rows'),:)=[];
end



%==================================================================
%
%Title:         ProjectElectrode2TransSurf
%Description:
%
%Inputs:
%==================================================================

% Project electrodes to the surface
%Input: number of vertices of the projected surface, matrix of vertices, all
%electrodes, electrodes that is going to be projected, Powerdata, frame index
%Output:electtrode vertices on the new surface, the color/power value of that electrode
function [elecVertices]=ProjectElectrode2TransSurf(vertices_num, vertices, all_electrodes, electrode)
elecVertices = zeros(size(electrode,1),1);
for elec_index = 1:size(electrode,1)
    xCoord = ones(vertices_num,1)*electrode(elec_index,2);
    yCoord = ones(vertices_num,1)*electrode(elec_index,3);
    zCoord = ones(vertices_num,1)*electrode(elec_index,4);
    [minVal, minInd] = min(sqrt((vertices(:,1) - xCoord).^2 + (vertices(:,2) - yCoord).^2 + (vertices(:,3) - zCoord).^2));
    elecVertices(elec_index, 1) = minInd;
end
% load('labels');
% overlay_frame_data = overlay_frame_data';
% electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);
end


%==================================================================
%
%Title:         create_electrode_montage
%Description:   This function takes a comma delimited sheet with all the
%               names of the electrodes and coordinates and reads it into
%               Matlab. It then saves the electrodes according to left and
%               right side with the (x,y,z) coordinates and the index of
%               that electrode in 'labels'. Essentially the electrode names
%               in 'labels' comes from the EEG .edf file and the coordinates
%               come from the .mgrid file so this reconciles the naming of
%               electrodes between the two
%
%Inputs:        "patient_folder"    [str] location of this patient's files
%               "flipX"             [int] 0 means it's the regular montage
%                                         1 means use the X-flipped montage
%==================================================================

function create_electrode_montage(patient_folder,flipX, IsVisual)
cd (patient_folder);

% if you are creating a montage for the x-flipped electrodes, you want to
% open and then save everything with a "_flipX" attached to it
if flipX == 0
    flipX_str = '';
elseif flipX == 1
    flipX_str = '_flipX';
end

% Read in Electrode locations (Map.xls or Map_flipX.xls) and Montage file (Montage.xls).
if IsVisual
    load('labels.mat');
else
    try
        load('labels_all_gray_matter.mat');
    catch
        load('labels_first_surgery_all_gray_matter.mat');
    end
end
if size(labels, 1) == 1
    labels = labels';
end
b1 = linspace(1, length(labels), length(labels))';
b1 = num2cell(b1);
b2 = labels(:,1);
montage = [b1 b2];
if IsVisual
    if ~flipX
        [a1, a2, map]=xlsread('Map.xls');
    else
        [a1, a2, map]=xlsread('Map_flipX.xls');
    end
else
    if ~flipX
        [a1, a2, map] = xlsread('Map_all_gray_matter.xlsx');
    else
        [a1, a2, map] = xlsread('Map_flipX_all_gray_matter.xlsx');
    end
end
% Convert the long channel names (i.e "A_L_Most_Mid_Frontal_Polar_12" into the abbreviated names found in the montage (i.e. "A12")
try
    for i=1:size(map,1)
        str             =   map{i,1};
        str(isspace(str)) = '';
        str_delim       =   regexp(str,'_','split');
        if length(str_delim)==1 || length(str_delim) == 2
            if startsWith(str_delim{1,1}, {'L', 'R'})
                map{i,5}=strtrim(str_delim{1,1}(1));
                map{i,1}=strtrim(str_delim{1,2});
                map{i,1}(map{i,1} == ';') = [];
            else
                map{i,5}=strtrim(str_delim{1,1}(end));
            end
        else
            firstLetter     =   str_delim{1};
            positionLetter  =   strtrim(str_delim{2});
            digit_split     =   split(str, ';');
            digitStr        =   digit_split{end};
            replaceStr      =   strcat(firstLetter, digitStr);
            map{i,1}        =   replaceStr; % i.e. "A12"
            map{i,5}        =   positionLetter; % Either L or R
        end
    end
catch ex
    disp(['Error reading Map file! There was a problem on line ' num2str(i) '. Verify that file is in the correct format (i.e each row is Letter_Side*Number. Unacceptable entries: J_Sup_Parietal8,PEG1)']);
    throw(ex);
end

Montage=montage;
%%
%     try
% Combine montage with 3-d electrode locations
for i=1:size(montage,1)
    MontageMap(i,1)=montage{i,1}; % Abbreviated Channel ("A12")
    for j=1:size(map,1)
        if strcmp(strtrim(montage{i,2}), map{j,1})
            if i == 113
                leah = [];
            end
            MontageMap(i,2)=map{j,2}; % X Coordinates
            MontageMap(i,3)=map{j,3}; % Y Coordinates
            MontageMap(i,4)=map{j,4}; % Z Coordinates
            %                     MontageMap(i,5) = newTimes{i,1}; % Arrival Time
            Position(i,1)=map{j,5}; % Side (R or L)
        end
    end
end
%         % catch
%         %     disp(['Error on MontageMap write! There was a problem on line ' num2str(i) '. Verify that file is in the correct format (i.e each row is Letter_Side*Number. Unacceptable entries: J_Sup_Parietal8)']);
%         %     error('MATLAB:myCode:dimensions', '');
%     end
j=1;k=1;

% check MontageMap to see if there are any errors that result in
% (0,0,0) coordinates for (x,y,z)
for i = 1:size(MontageMap,2)
    if MontageMap(i,1) == 0 && MontageMap(i,2) == 0 && MontageMap(i,3) == 0
        display(MontageMap{i,1});
        error('This electrode name in Map or Map_flipX is not matching, double check the original Map sheet for misnamed electrodes compared to labels');
    end
end


% Break the montage into L and R sides
for i=1: length(Position)
    if (Position(i)=='L')
        L_MontageMap(j,1:4)=MontageMap(i,1:4);
        j=j+1;
    else
        R_MontageMap(k,1:4)=MontageMap(i,1:4);
        k=k+1;
    end
end

if (j==1)
    L_MontageMap=[];
elseif (k==1)
    R_MontageMap=[];
end
% Save data into mat files
MontageMap=[L_MontageMap; R_MontageMap];
MontageMap=sortrows(MontageMap,1);
save(['Montage' flipX_str '.mat'], 'Montage');
save(['MontageMap' flipX_str '.mat'], 'MontageMap');
save(['L_MontageMap' flipX_str '.mat'], 'L_MontageMap');
save(['R_MontageMap' flipX_str '.mat'], 'R_MontageMap');

end

%==================================================================
%
%Title:         electrode_data_overlay
%Description:   This function loads a surface txt file output from Bioimage
%               Suite into Matlab, and arranges it into a readable format.
%               Repeat twice (once for cortex and once for smooth). The
%               output is a reorganized matrix of vertices.
%
%Inputs:        "vertices_num"  [int] The number of vertices in the file
%               "faces_num"     [int] The number of faces in the file
%               "data_location" [string] should point to the root folder
%                               for the seizure
%                               i.e. 'C:\IcEEGProject_ProcessedData\25\25_1'
%               "filename"      [string] Location of the vertice data
%==================================================================

function [electrode_vertex_values] = electrode_data_overlay(overlay,frame_index,all_electrodes, file_suffix, IsVisual, final_bin, location)
overlay_frame_data = [];
if any(strcmp(overlay,  {'Electrode Distribution', 'Arrival Times', 'Electrode Density'}))
    
    % Give every overlay value a 1
    % once aggregated across patients it gives the density map of
    % electrodes at each vertex
    if IsVisual
        load('labels.mat');
    else
        try
            load('labels_all_gray_matter.mat');
        catch
            load('labels_first_surgery_all_gray_matter.mat');
        end
    end
    overlay_frame_data = ones(length(labels),1);
    
elseif any(strcmp(overlay, {'CP', 'CNP', 'Subtraction'}))
    subjIDprep1=strfind(location, 'EEG_behavior/')+13;
    subjIDprep2=strfind(location, '/icEEG/')-1;
    subjID=location(subjIDprep1:subjIDprep2);
    location2=['V:\auditory\database\' subjID];
    % getting power zscores from mean power of each electrode
    %\\172.23.254.106\Data25\auditory\database\SUBJID\meanpower_traces_157bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_2to2_baseline5to0_commonr.mat
    file_suffix2=erase(file_suffix, 'from');
    if isfile([location '\meanpower_traces_157bins_zscore_recentered_rejoutliers', file_suffix2,  '.mat'])
        load([location '\meanpower_traces_157bins_zscore_recentered_rejoutliers', file_suffix2,  '.mat'])
    elseif isfile([location '\meanpower_traces_156bins_zscore_recentered_rejoutliers', file_suffix2,  '.mat'])
        load([location '\meanpower_traces_156bins_zscore_recentered_rejoutliers', file_suffix2,  '.mat'])
    elseif isfile([location2 '\meanpower_traces_157bins_zscore_recentered_rejoutliers', file_suffix2,  '.mat'])
        load([location2 '\meanpower_traces_157bins_zscore_recentered_rejoutliers', file_suffix2,  '.mat'])
    else
        load([location2 '\meanpower_traces_156bins_zscore_recentered_rejoutliers', file_suffix2,  '.mat'])
        
    end
    %     load(['meanpower_traces' file_suffix, '.mat'])
    
    if strcmp(overlay, 'CP')
        trialtype = 1; % confirmed perceived
    elseif strcmp(overlay, 'CNP')
        trialtype = 2; % confirmed not perceived
    end
    if ~strcmp(overlay, 'Subtraction')
        overlay_frame_data = squeeze(meanpower_traces(:, trialtype, 3, :));
    else
        overlay_frame_data = squeeze(meanpower_traces(:, 1, 3, :)) - squeeze(meanpower_traces(:, 2, 3, :));
    end
end

% once you have the overlay data, you need to match it to an electrode from
% labels
overlay_frame_data = overlay_frame_data';
electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);
end

%% Determine the transparency for every face/vertex

%==================================================================
%
%Title:         CalculateTransparency
%Description:   This function takes the input values from overlay frame
%               data for each electrode and gives it to all the vertices
%               within radius2 distance from each electrodes vertex
%
%Inputs:        "vertices_num"  The number of vertices in the file
%               "elecVertices"  vertex indices of each electrode to be
%                               plotted on this side of the brain
%               "electrode_vertex_values"
%                               the overlay value assigned to each of those
%                               vertices in elecVertices
%               "overlay"      	[int] this number is an input to the parent
%                               function. If overlay = 0 then it does not
%                               do any weighing based on other electrodes
%                               within a radius2 proximity. All vertices
%                               within radius 2 receive a 1
%
%Outputs:       "vertexCdata"   matrix the length of the number of vertices
%                               with the electrode_vertex_values assigned
%                               to all vertices within radius2 of each
%                               electrode
%==================================================================

% This function should be combined with OverlayColor to speed up processing
% time...

function [vertexCdata] = CalculateTransparency(vertices_num, vertices, elecVertices, electrode_vertex_values, overlay)

maxAlpha = 1;
radius1 = 1;
radius2 = 15;
electrodes_num = length(elecVertices);
electrode_vector = vertices(elecVertices, :);
vertexCdata = zeros(vertices_num, 1);

for vertIndex = 1 : vertices_num
    
    xCoord = ones(electrodes_num, 1) * vertices(vertIndex, 1);
    yCoord = ones(electrodes_num, 1) * vertices(vertIndex, 2);
    zCoord = ones(electrodes_num, 1) * vertices(vertIndex, 3);
    
    % Calculate "Transparency" in other words the "weight" based on
    % distance from the vertex
    [distanceToElectrode, electrodeIndice] = sort(sqrt((electrode_vector(:, 1) - xCoord) .^ 2 + (electrode_vector(:, 2) - yCoord) .^ 2 + (electrode_vector(:, 3) - zCoord) .^ 2));
    
    % Initialize transparency and electrode index vectors for this vertex
    
    % Wendys method
    sharedTransparency = [];
    sharedElectrode = [];
    
    %if any(strcmp(overlay, {'Electrode Distribution', 'Electrode Density'}))
    if any(strcmp(overlay, {'Electrode Distribution'}))
        if any(distanceToElectrode <= radius2)
            vertexCdata(vertIndex,1) = 1;
        end
    else
        % Get electrode indices corresponding to different vertex proximity conditions
        for n = 1 : length(electrodeIndice)
            if distanceToElectrode(n) <= radius1
                
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                
                sharedTransparency = [sharedTransparency; maxAlpha];
                
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
                
            elseif distanceToElectrode(n) <= radius2
                
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                sharedTransparency = [sharedTransparency; maxAlpha - (((distanceToElectrode(n) - radius1) / (radius2 - radius1)) * maxAlpha)];
                
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
            end
        end
        
        % Get electrode values which contribute to this vertex
        zIn = electrode_vertex_values(sharedElectrode)';
        
        % Calculate weighted values for each electrode dependent on the distance to the vertex
        % Wendys old and slow method
        %weightedZ = [];
        %for h = 1:length(zIn)
        %    weightedZ = [weightedZ zIn(h) * (sharedTransparency(h))];
        %end
        
        % Noahs new and fast method which is mathematically identical
        weightedZ = zIn .* sharedTransparency;
        
        
        % Aggregate electrode values at the vertex
        % Wendys electrode aggregation method: sum of weighted electrode values for the vertex
        weightedZ = nansum(weightedZ);
        
        % If there is no weighted zScore that means the vertex was greater than
        % 15mm from every electrode and has no power or zScore to display in
        % this frame
        if isempty(zIn)
            vertexCdata(vertIndex, 1) = 0;
        elseif isnan(zIn)
            vertexCdata(vertIndex, 1) = 0;
        else
            vertexCdata(vertIndex, 1) = weightedZ;
        end
    end
end
end
