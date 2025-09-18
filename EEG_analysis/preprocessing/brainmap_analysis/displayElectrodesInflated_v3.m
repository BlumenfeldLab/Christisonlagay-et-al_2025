function displayElectrodesInflated_v3(data_location, patients, overlay,...
    createMontage, laterality, views, inflationstep, file_suffixes, varargin)

%{
Chris Micek Edits

  Plots patient electrode locations and gamma power/electrode density
  overlays onto the inflated MNI brain. Requires adding the 'Functions for
  Inflated Brain Display' folder and its subfolders to the MATLAB path for
  correct operation. Also requires that each patient folder in 'data_location'
  contain the files 'labels[_all_gray_matter].mat' and
  'meanpower_traces_57bins_zcore_recentered_rejoutliers{file_suffix}.mat'.

  If running a patient for the first time, ensure that 'createMontage' is
  1.

  UPDATED 7/2019 from displayElectrodesInflated.m to increase robustness when
  faced with different combinations of user input, and adds some new
  functionality. Changes are summarized as follows:

  1. Instead of 'ElectrodeDensity' and 'Subtraction' being their own
     input arguments, they now have their own separate overlay codes.

  2. Different overlays now have both string and numeric codes, to make it
     easier to keep track of the different available options. Using the
     string codes is preferred, but numeric codes can still be used; they
     now start at 1 instead of 0 so they are consistent with MATLAB
     indexing.

  3. The code works with old visual task data without alteration; supplying
     the 'IsVisual' optional argument specifies which patients have visual
     data. See 'IsVisual' in the OPTIONAL INPUTS section below for more
     details.

  4. Adds a 'Comparison' overlay that can be used to compare two datasets;
     see OPTIONAL INPUTS below for more info.

  5. The 'image_storage_folder' optional input argument is still present; if
     used, it must be the first optional argument. All others are given as
     key-value pairs.

INPUTS:
  REQUIRED:
    - data_location: Location of patient data folders for brain surface
                     overlays. On Chris's computer, this is 'E:\Chris
                     HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP
                     movies'.

    - patients: Cell array of patient ID codes. There must be a corresponding
                folder in 'data_location' for each.

    - overlay: The type of overlay to use for the brain surface. This can
               be supplied as either a string or numeric code. String codes
               are preferred; numeric codes now begin at 1 to correspond
               with MATLAB indexing, and to avoid confusion with the
               previous syntax a prompt will appear asking if you are sure
               of your input before continuing whenever they are used.

                 'Electrode Distribution'/1 = nothing, just draw electrodes 
                                              on blank brain. Saves to 
                                              'data_location'\Group\Electrode_Distribution

                 'Arrival Times'/2 = arrival times

                 'CP'/3 = zScores confirmed perceived. Saves to
                          'data_location'\Group\Gamma_power_CP

                 'CNP'/4 = zscores confirmed not perceived. Saves to
                           'data_location'\Group\Gamma_power_CNP

                 'Electrode Density'/5 = overlay electrode coverage
                                         density along with electrode
                                         locations. You will be prompted 
                                         to select a 'vertex_values.mat'
                                         file containing power
                                         contributions from the electrodes
                                         you want to include. Saves images to
                                         'data_location'\Group\Electrode_Density.

                 'Subtraction'/6 = overlay subtraction of confirmed not perceived
                                   power from confirmed perceived power (CP - CNP).
                                   Saves images to
                                   'data_location'\Group\Gamma_power_Subtraction.

                 'Comparison'/7 = overlay of a measure of comparison
                                  between two datasets, specified in the
                                  'CompareMethod' optional argument. The
                                  vertex value file paths for the two datasets
                                  to be compared must be specified in
                                  the 'VertexVals1' and 'VertexVals2' optional
                                  arguments. Note that it is assumed that,
                                  if 'VertexVals1' contains data from M
                                  patients, and 'VertexVals2' contains data
                                  from N patients, the first M patients in
                                  'patients' have data from 'VertexVals1',
                                  and the remaining N patients have data
                                  from 'VertexVals2'. Saves images to
                                  'data_location'\Group\Gamma_power_Comparison.

   - createMontage: Should be 1 if a never-before analyzed patient is
                    included in 'patients':
                      1 = create each patient's montage and map for EEG
                          channels (i.e. read in and save the electrode names
                          and locations from their 'Map... .xlsx' file)
                      0 = nothing

    - laterality: Which hemisphere(s) to plot overlay data onto:
                    0 = bilateral
                    1 = left only
                    2 = right only
                    3 = project both left and right electrodes on the left
                        brain (flipX)

    - views: Which views of the inflated surface to include in the saved
             figures:
               4 = 4 views (lateral, medial, ventral, and posterior)
               2 = 2 views (Left & right lateral/medial)
               'Lateral' = Left & right lateral
               'Medial' = Left & right medial

    - inflationstep: How inflated the brain surface should be; options range
                     from 1 (least inflated) to 6 (most inflated)

    - file_suffixes: Cell array indicating the file suffixes to look for when
                     loading each patient's
                     'meanpower_traces_57bins_zcore_recentered_rejoutliers{file_suffix}.mat'
                     file. Each index of 'file_suffixes' must correspond with
                     the patient at the same index in 'patients'.

  OPTIONAL:
    If used, the 'image_storage_folder' argument below must be the first
    optional argument. All others take the form of key-value pairs.
    - image_storage_folder: Absolute path to a folder where images should
                            be saved. If specified, overrides default image
                            save location and uses the new location
                            instead.

    The options below are provided in the form of key-value pairs:
    - 'IsVisual': A logical vector the same size as 'patients', specifying
                  which patients in 'patients' have visual task data that
                  is being analyzed. If an element is true, the code runs
                  the same way Wendy's old code did for that particular
                  patient. Specifically,
                  1. the code looks for 'labels_depth_most.mat' in the
                     patient folder, and discards the depth electrodes
                     specified for analysis, and
                  2. if 'createMontage' is true, looks for 'labels.mat' in
                     the patient folder instead of
                     'labels_all_gray_matter.mat' when searching for which
                     electrode contacts to include.
                  
                  If 'IsVisual' is not provided, the default value is a
                  vector of false values for each patient.

    - 'CompareMethod': If 'overlay' is 'Comparison', you must specify one
                       of the following comparison methods:
                       - 'XCorr': Uses the MATLAB xcorr function to compute
                                  the normalized cross-correlation of the
                                  data in 'VertexVals1' with respect to the
                                  data in 'VertexVals2', using the 'coeff'
                                  option so all values are between -1 and
                                  +1. The output images are marked with
                                  time lags with respect to the start of
                                  the second dataset. NOTE that
                                  cross-correlation is not a commutative
                                  operation; switching the order of the
                                  operands will result in a reversed image
                                  sequence.

                       - 'SignProdSum': For each time point, takes the
                                        product of the signs of both
                                        datasets, and multiplies this with
                                        the aggregated sums of the absolute
                                        values of the z-scores for both
                                        datasets (abs(z-score) sums /
                                        sqrt(number of patients
                                        contributing to that vertex). This
                                        results in positive values where
                                        the datasets have the same
                                        polarity, and negative values where
                                        they are opposite. The magnitudes
                                        at each vertex follow the same
                                        weighting as the normal processing
                                        pipeline for group data, assuming
                                        values for the two datasets are
                                        part of one group.

    - 'VertexVals1': If 'overlay' is 'Comparison', a string with the
                     absolute/relative path to the 'vertex_values.mat' file
                     containing the vertex data of the first dataset to
                     compare. If the operand order of 'CompareMethod' is
                     important, 'VertexVals1' will be the first argument.

    - 'VertexVals2': If 'overlay' is 'Comparison', a string with the
                     absolute/relative path to the 'vertex_values.mat' file
                     containing the vertex data of the second dataset to
                     compare. If the operand order of 'CompareMethod' is
                     important, 'VertexVals2' will be the second argument.

SAVES:
  To the 'image_storage_folder' (default locations if optional parameter is
  unused):
    - Images for each brain view individually, each hemisphere with its
      views, and all views combined, with the overlay specified. If overlay
      is 2 or 3, contains pictures for each time bin of
      'meanpower_traces_57bins_zscore_recentered_rejoutliers{file_suffix}.mat'.

    - 'vertex_values.mat': Contains the variables 'vertex_valuesL' and
                           'vertex_valuesR', each vertices x patients x
                           time_bins matrices of overlay values for each
                           vertex of the respective hemisphere of the
                           inflated brain surface.

  If 'createMontage' is 1, additionally saves the following to each
  patient's brain overlay data folder (if 'laterality' is 3, appends
  the '_flipX' suffix to each file name). If 'createMontage' is 0, relies
  on the following files being present in each patient's overlay data folder:
    - 'Montage.mat': Contains the 'Montage' variable, a channels x 2 cell
                     array, where the first column is the indices of
                     included channel labels from
                     'labels_all_gray_matter.mat', and the second contains
                     the corresponding label.

    - 'MontageMap.mat': Contains 'MontageMap', a channels x 4 double array,
                        where the first column is the label index from
                        'Montage.mat', and columns 2 - 4 are that contact's x,
                        y, and z coordinates in voxel-space.

    - 'L_MontageMap.mat': Contains the data from 'MontageMap.mat' belonging
                          to the left hemisphere.

    - 'R_MontageMap.mat': Contains the data from 'MontageMap.mat' belonging
                          to the right hemisphere.
    
%}

% list of all patients
% patients = {'193AF','154NP','200JW','201MU','210VG','212DD','213SP','215JF','217CG'};
% patients ={'311KH', '320JG2', '383TS', '390DI', '404LL', '415DO', '436BP', '418JB', 'NS140', '444VS', '445JK', '446LO', '451NH', '468SZ'}
% displayElectrodesInflated(data_location, patients, overlay, createMontage, laterality, views, inflationstep, electrodeDensity, subtraction, file_suffixes)
%
% displayElectrodesInflated('E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies', {'311KH', '320JG2', '383TS', '390DI', '404LL', '415DO', '436BP', '418JB', 'NS140', '444VS', '445JK', '446LO', '451NH'}, 3, 0, 0, 4, 5, false, false, file_suffixes)
% image_storage_folder = [data_location '/Group/Wendy picture storage2'];

overlay_labels = {'Electrode Distribution', 'Arrival Times', 'CP', 'CNP', 'Electrode Density', 'Subtraction', 'Comparison'};

if isnumeric(overlay)
    disp('NOTE: Overlay numeric codes now begin at 1, not 0. Please ensure your overlay value is still correct!')
    inptchk = input(['You typed ' num2str(overlay) ' as your overlay code. Is this correct? Type another to change it,',...
        ' or leave your answer blank to keep it as-is:\n']);
    if ~isempty(inptchk)
        overlay = overlay_labels{inptchk};
    else
        overlay = overlay_labels(overlay);
    end
else
    if ~any(strcmp(overlay, overlay_labels))
        throw(MException('displayElectrodesInflated:InvalidOverlay',...
            [sprintf('Error: ''%s''', overlay), ' is an invalid overlay. Valid overlays are:',...
            sprintf(' ''%s''', overlay_labels{:})]))
    end
end

switch overlay
    case 'Electrode Distribution'
        image_storage_folder = [data_location '/Group/Electrode_Distribution'];
    case 'CP'
        image_storage_folder = [data_location '/Group/Gamma_power_CP'];
    case 'CNP'
        image_storage_folder = [data_location '/Group/Gamma_power_CNP'];
    case 'Electrode Density'
        image_storage_folder = [data_location '/Group/Electrode_Density'];
    case 'Subtraction'
        image_storage_folder = [data_location '/Group/Gamma_power_Subtraction'];
    case 'Comparison'
        image_storage_folder = [data_location '/Group/Gamma_power_Comparison'];
end

valid_keywords = {'IsVisual', 'VertexVals1', 'VertexVals2', 'CompareMethod'};

valid_compare_methods = {'XCorr', 'SignProdSum'};

params.IsVisual = false(1, numel(patients));
params.VertexVals1 = [];
params.VertexVals2 = [];
params.CompareMethod = '';

if exist(varargin{1}, 'dir')
    image_storage_folder = varargin{1};
    varargin(1) = [];
end

if mod(numel(varargin), 2)
    throw(MException('displayElectrodesInflated:MismatchedKeyword',...
        'Error: Number of keywords and values is not equal!'))
end

if numel(varargin) > 1
    keywords = varargin(1:2:end);
    vals = varargin(2:2:end);
    
    for k = 1:numel(keywords)
        if ismember(keywords{k}, valid_keywords)
            params.(keywords{k}) = vals{k};
        else
            throw(MException('displayElectrodesInflated:InvalidKeyword',...
                [sprintf('Error: ''%s''', keywords{k}), ' is an invalid keyword argument.']))
        end
    end
end

if strcmp(overlay, 'Comparison') && isempty(params.VertexVals1)
    throw(MException('displayElectrodesInflated:MissingVertexVals',...
        'Error: ''Comparison'' overlay selected but VertexVals1 missing.'))
end

if strcmp(overlay, 'Comparison') && isempty(params.VertexVals2)
    throw(MException('displayElectrodesInflated:MissingVertexVals',...
        'Error: ''Comparison'' overlay selected but VertexVals2 missing.'))
end

if strcmp(overlay, 'Comparison') && isempty(params.CompareMethod)
    throw(MException('displayElectrodesInflated:MissingVertexVals',...
        'Error: ''Comparison'' overlay selected but CompareMethod missing.'))
end

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

xtick_spacing_filename = 'spectrogram_times_57bins.mat';
mni2fs_dir = 'E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master';
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);

% file_suffix = '';
% file_suffix = '_sounds_restricted';

mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);

Tmni = mnit1_1mm_stripped.transform;

alpha = 0.7;
zScoreThresh = 2;

% list of patients with only right sided contacts

% VTK_folder_name = 'VTK Leah';

% the colors for the different patients' electrodes
% colors = {[0    0.4470    0.7410],...% blue - 311KH
%     [0.8500    0.3250    0.0980],... % red orange - 320JG2
%     [0.9290    0.6940    0.1250],... % gold - 383TS
%     [0.4940    0.1840    0.5560],... % purple - 390DI
%     [0.2100    0.3700    0.2300],... % olive green - 404LL
%     [0.3010    0.7450    0.9330],... % cyan - 415DO
%     [0.6350    0.0780    0.1840],... % maroon - 436BP
%     [0.8745    0.2353    0.3882],... % fuschia - 418JB
%     [0         1         0],...      % lime green - NS140
%     [1.00      0.00      1.00],...   % magenta - 444VS
%     [0.95      0.95      0.00],...   % bright yellow - 445JK
%     [0.2510    0.8784    0.8157],... % turquoise - 446LO
%     [1         1         1],...      % white - 451NH
%     [0         0         0],...      % black - 468SZ
%     [0.2510    0.1647    0.6784],... % royal purple - 477NS
%     [0.6980    0.7176    0.9765],... % periwinkle blue - 444VS2
%     [0.1804    0.3608    0.3490]};   % blue spruce - 485DG

if ~strcmp(overlay, 'Comparison')
    colors = distinguishable_colors(length(patients), .5 * ones(1, 3));
end

%% Side laterality - which side of the brain to draw
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

figHandle = gobjects(1, sides_num);

%% hard coded number of frames for each overlay setting, need to change later

if any(strcmp(overlay, {'Electrode Distribution', 'Arrival Times', 'Electrode Density'}))
    total_num_frames = 1;
elseif any(strcmp(overlay, {'CP', 'CNP', 'Subtraction', 'Comparison'}))
    load(xtick_spacing_filename)
    if strcmp(overlay, 'Comparison') && strcmp(params.CompareMethod, 'XCorr')
        total_num_frames = 2 * numel(T) - 1;
        t_interval = T(2) - T(1);
        T2 = -(numel(T) - 1) * t_interval:t_interval:(numel(T) - 1) * t_interval;
    else
        total_num_frames = numel(T);
    end
end

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
    
    %------------------------------------------------------------------
    %   Set-up some local variables
    %------------------------------------------------------------------
    
    surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);
    fs_surf(side_index) = export(gifti(surf_fn));
    
    v_num = size(fs_surf(side_index).vertices, 1);
    
    fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] *...
        Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
    fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
    
    surfrender_fn = fullfile(mni2fs_dir,['/surf/' hem '.inflated' num2str(inflationstep) '.surf.gii']);
    inflated_surf(side_index) = export(gifti(surfrender_fn));
    inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] *...
        Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
    inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);
    
    %     cortex_filename{side_index}  =    [side '_cortex.txt'];
    %     smooth_filename{side_index}  =    [side '_smooth.txt'];
    %     cd([data_location '\' VTK_folder_name]);
    %     cortex_v_num(side_index)     =    eval(['metadata.' side '_cortex_v_num']);
    %     cortex_f_num(side_index)     =    eval(['metadata.' side '_cortex_f_num']);
    %     smooth_v_num(side_index)     =    eval(['metadata.' side '_smooth_v_num']);
    %     smooth_f_num(side_index)    =    eval(['metadata.' side '_smooth_f_num']);
    %
    %     % Input and interpret the vertice data for the cortex and smooth
    %     % surfces
    %     [cortex(side_index).vertices, cortex(side_index).faces] = FormatMeshSurface(cortex_v_num(side_index), cortex_f_num(side_index), data_location, cortex_filename{side_index});
    %     [smooth(side_index).vertices, smooth(side_index).faces] = FormatMeshSurface(smooth_v_num(side_index), smooth_f_num(side_index), data_location, smooth_filename{side_index});
    
    
    %% create a figure of the left or right side of the brain
    % you will be working with this for different electrodes, overlays,
    % views, and light
    
    figHandle(side_index) = figure('Position',[70,70,700,700]); % Original = [50,50,600,600]Creates a figure graphic object
    %     g=trisurf(cortex.faces,cortex.vertices(:,1),cortex.vertices(:,2),cortex.vertices(:,3));
    
    temp_surf = [];
    temp_surf.hem = hem; % choose the hemesphere 'lh' or 'rh'
    temp_surf.inflationstep = inflationstep; % 1 no inflation, 6 fully inflated
    temp_surf.decimation = 0;
    temp_surf = mni2fs_brain(temp_surf);
    set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
    
    %     g = trisurf(inflated_surf(side_index).faces,inflated_surf(side_index).vertices(:,1),...
    %         inflated_surf(side_index).vertices(:,2),inflated_surf(side_index).vertices(:,3));
    %
    %     lighting flat;
    %     set(g, 'FaceColor', [1 1 1], 'LineStyle', 'none')%, 'EdgeLighting','gouraud')%, 'FaceLighting', 'gouraud');
    set(gcf, 'color', [1 1 1]);
    hold on
    %     xlabel('X');
    %     ylabel('Y');
    %     zlabel('Z');
    
    axis off;
    axis equal;
    vertex_values = zeros(v_num, length(patients));
    
    % if you are plotting arrival times, then each frame is a different set
    % of electrodes
    if strcmp(overlay, 'Arrival Times')
        electrode_mapping_times = 10;
    else
        electrode_mapping_times = 1;
    end
    
    for e = 1:electrode_mapping_times
        scatter_count = 0;
        for p = 1:length(patients)
            
            %% Get electrode info for the patient
            
            % this part plots the electrodes for each patient and the overlay data
            % for each patient
            patient_folder = [data_location, '/', patients{p}];
            
            %% Make the Montage
            % if the left/right electrode montages have not been created, do so
            % now, indicated by createMontage = 0 or = 1
            if createMontage == 1
                % in this case, just do the normal thing: create the normal
                % (non-X-flipped) montage for left and right
                flipX = 0;
                create_electrode_montage(patient_folder,flipX, params.IsVisual(p))
                if  laterality == 3 % for x-flipped in addition to regular
                    % in this case, create the x-flipped montage for the left and
                    % right electrodes but you will later only use the RIGHT
                    % electrodes with the coordinates flipped to the LEFT
                    flipX = 1;
                    create_electrode_montage(patient_folder,flipX, params.IsVisual(p))
                end
            end
            
            %         % Combine EEG Channel with Electrode Position
            %         cd(data_location)
            %
            %         [overlay_frame_data] = OrganizeElectrodePosition(data_location, patient_folder, p, overlay, createMontage, overlay_frame_data);
            
            display(patients{p})
            
            % Now load some data files needed for electrode placement and surface shading
            cd([data_location, '/', patients{p}]);
            load([side '_MontageMap.mat']);
            electrode = eval([side '_MontageMap' ]);
            
            %         electrode_num = size(electrode,1);
            % load the X-flipped montage and add to list of electrodes
            if laterality == 3
                load('R_MontageMap_flipX.mat');
                electrode = cat(1,electrode,R_MontageMap);
            end
            
            %% Discard depth electrodes or plot only some electrodes for arrival times
            if params.IsVisual(p) && ~isempty(electrode)
                % find out which depth electrodes should be discarded from the
                % analysis
                if ~strcmp(patients{p}, '193AF')
                    load('labels.mat')
                    load('labels_depth_most.mat')
                else
                    load('labels_first_surgery.mat')
                    labels_depth = {};
                end
                
                % find the indices of the white matter electrodes in labels
                labels_depth_idx = [];
                for l = 1:length(labels_depth)
                    labels_depth_idx = [labels_depth_idx find(strcmp(labels_depth{l},labels))];
                end
                
                % remove them from the electrode list that you will later access
                [side_labels_depth, depth_idx, electrode_idx] = intersect(labels_depth_idx,electrode(:,1));
                electrode(electrode_idx,:) = [];
                
                % if you are trying to plot arrival times, then you will need
                % to plot only the electrodes in that time period
                if strcmp(overlay, 'Arrival Times')
                    load('arrival_times_by_period_100ms.mat')
                    electrodes_this_period = find(arrival_times_by_period(:,e) == 1);
                    electrode = electrode(ismember(electrode(:,1),electrodes_this_period),:);
                end
            end
            
            %% Plot patient's electrodes
            
            %             electrode = [electrode(:, 1), [electrode(:, 2:end), ones(size(electrode, 1), 1)] * Tmni'];
            %             electrode = electrode(:, 1:4);
            
            % assign this patient's electrodes one of the colors
            %             colorElectrode = colors{p};
            
            if ~strcmp(overlay, 'Comparison')
                colorElectrode = colors(p, :);
            else
                if ismatrix(params.VertexVals1.vertex_valuesL)
                    num_pts_1 = 1;
                else
                    num_pts_1 = size(params.VertexVals1.vertex_valuesL, 2);
                end
                    
                if p <= num_pts_1 
                    colorElectrode = 'b';
                else
                    colorElectrode = 'r';
                end
            end
            
            if isempty(electrode) % Some files only have electrodes on the one side
                disp('-- No electrode data!');
                electrodes_present(side_index,p) = 0;
            else
                all_electrodes = electrode(:,1);
                electrodes_present(side_index,p) = 1;
                
                %                 [elecVertices2]= ProjectElectrode2TransSurf(smooth_v_num, smooth(side_index).vertices, all_electrodes, electrode, side_index);
                %                 electrode(:,2:4)                       = smooth(side_index).vertices(elecVertices2,1:3);
                [elecVertices] = ProjectElectrode2TransSurf(v_num, fs_surf(side_index).vertices, all_electrodes, electrode);
                
                %Plot the electrodes in a specified color and count them
                %in case you need to remove them later
                scatter3(inflated_surf(side_index).vertices(elecVertices,1),...
                    inflated_surf(side_index).vertices(elecVertices,2),...
                    inflated_surf(side_index).vertices(elecVertices,3), 40, colorElectrode,'filled');
                scatter_count = scatter_count + 1;
            end
            
            %% cycle through all the frames to aggregate zscores for each vertex
            if any(strcmp(overlay, {'CP', 'CNP', 'Subtraction', 'Comparison'}))
                for frame_index = 1:total_num_frames
                    
                    if isempty(electrode) % Some files only have electrodes on one hemisphere, hence the if statement
                        disp('-- No electrode data!');
                        FaceVertexAlphaData    = zeros(v_num, 1);
                        vertexCdata            = zeros(v_num, 1);
                        %             frames_set(side_index,1,:) = GenerateFrame(cortex,FaceVertexAlphaData,vertexCdata,[],color_map,colorrange,side,0);
                    else
                        %                 all_electrodes = electrode(:,1);
                        
                        %                 [elecVertices2,electrode_vertex_values]= ProjectElectrode2TransSurf(smooth_v_num, smooth.vertices, all_electrodes, electrode, frame_index, side_index, overlay_frame_data);
                        %                 electrode(:,2:4)                       = smooth.vertices(elecVertices2,1:3);
                        %                 [elecVertices,electrode_vertex_values] = ProjectElectrode2TransSurf(cortex_v_num, cortex.vertices, all_electrodes, electrode, frame_index, side_index, overlay_frame_data);
                        
                        display([patients{p} ' Aggregating Frame ' num2str(frame_index)])
                        
                        if ~strcmp(overlay, 'Comparison')
                            [electrode_vertex_values] = electrode_data_overlay(overlay,frame_index,all_electrodes, file_suffixes{p}, params.IsVisual(p));
                            [vertexCdata]        = CalculateTransparency(v_num, inflated_surf(side_index).vertices, elecVertices, electrode_vertex_values, overlay);
                            
                            %                 [FaceVertexAlphaData,ColorFace]        = CalculateTransparency(cortex_v_num,cortex.vertices, elecVertices);
                            %                 [vertexCdata, shared_vertices]         = OverlayColor(cortex_v_num,cortex.vertices,elecVertices,electrode_vertex_values, ColorFace,p);
                            
                            % compile the value and transparency assigned to each
                            % vertex across patients
                            %                 transparency(:,p,frame_index) = FaceVertexAlphaData;
                            vertex_values(:,p,frame_index) = vertexCdata;
                            
                            %                 %Plot the electrodes in a specified color
                            %                 if frame_index == 1
                            %                     scatter3(cortex.vertices(elecVertices,1),cortex.vertices(elecVertices,2),cortex.vertices(elecVertices,3),40, colors{p},'filled');
                            %                 end
                            
                            %             h = patch('Faces',cortex.faces, 'Vertices', cortex.vertices, 'FaceVertexCData',vertexCdata,'FaceColor','interp','FaceAlpha','interp',...
                            %                 'FaceVertexAlphaData', FaceVertexAlphaData, 'LineStyle', 'none', 'FaceLighting','gouraud', 'EdgeLighting', 'gouraud');
                            
                        end
                        
                    end
                end
            end
        end
        
        %% Frame Generation if saving without overlaying color onto the brain
        
        % if you are doing arrival times, you need to save the figures at
        % the end of every time period and redraw the electrodes
        if strcmp(overlay, 'Arrival Times')
            disp(e)
            % Grab the figure of the side of the brain you are working on
            figure(figHandle(side_index))
            
            cd(image_storage_folder)
            frames_set(side_index,e,:) = savebrainimages(side,e);
            
            % remove the electrodes to make the next one
            children = get(gca,'children');
            delete(children(1:scatter_count,1))
        end
        
    end
    
    if ~strcmp(overlay, 'Arrival Times')
        %------------------------------------------------------------------
        %   Frame Generation with overlaying color onto the brain surface
        %------------------------------------------------------------------
        
        % go to the correct folder
        cd(image_storage_folder)
        load('cmap6to8WRX.mat');
        
        if ~any(strcmp(overlay, {'Electrode Density', 'Comparison'}))
            vertex_values_sum = squeeze(nansum(vertex_values, 2));
            divisor = squeeze(sum(vertex_values ~= 0, 2));
            divisor(divisor == 0) = 1;
            vertex_values_sum = vertex_values_sum ./ sqrt(divisor);
            
            % load the colormap and set axis of colormap
            colormap(mycmap)
            caxis([-6 8]);
            
        elseif strcmp(overlay, 'Electrode Density')
            if side_index == 1
                [filename, pathname] = uigetfile({'*.mat'; 'MAT-files'}, 'Select vertex values .mat file');
                fullfilename = strcat(pathname,filename);
                load(fullfilename)
                vertex_valuesL = squeeze(vertex_valuesL(:, :, 1));
                vertex_valuesL(isnan(vertex_valuesL)) = 0;
                vertex_valuesL(vertex_valuesL ~= 0) = 1;
                vertex_values_sum_L = squeeze(nansum(vertex_valuesL, 2));
                vertex_values_sum = vertex_values_sum_L;
                
                if laterality ~= 3
                    vertex_valuesR = squeeze(vertex_valuesR(:, :, 1));
                    vertex_valuesR(isnan(vertex_valuesR)) = 0;
                    vertex_valuesR(vertex_valuesR ~= 0) = 1;
                    vertex_values_sum_R = squeeze(nansum(vertex_valuesR, 2));
                    
                    colormap jet
                    caxis([min([vertex_values_sum_L; vertex_values_sum_R]),...
                        max([vertex_values_sum_L; vertex_values_sum_R])])
                else
                    colormap jet
                    caxis([min(vertex_values_sum_L(:)), max(vertex_values_sum_L(:))])
                end
                
            else
                vertex_values_sum = vertex_values_sum_R;
                colormap jet
                caxis([min([vertex_values_sum_L; vertex_values_sum_R]),...
                    max([vertex_values_sum_L; vertex_values_sum_R])])
            end
            
            
        elseif strcmp(overlay, 'Comparison')
            if strcmp(params.CompareMethod, 'XCorr')
                
                colormap jet
                caxis([-1 1])
                
                if side_index == 1
                    vertex_values_sum1 = squeeze(nansum(params.VertexVals1.('vertex_valuesL'), 2));
                    divisor1 = squeeze(sum(params.VertexVals1.('vertex_valuesL') ~= 0, 2));
                    divisor1(divisor1 == 0) = 1;
                    vertex_values_sum1 = vertex_values_sum1 ./ sqrt(divisor1);

                    vertex_values_sum2 = squeeze(nansum(params.VertexVals2.('vertex_valuesL'), 2));
                    divisor2 = squeeze(sum(params.VertexVals2.('vertex_valuesL') ~= 0, 2));
                    divisor2(divisor2 == 0) = 1;
                    vertex_values_sum2 = vertex_values_sum2 ./ sqrt(divisor2);
                    
                    vertex_values_sum_L = NaN(size(vertex_values_sum1, 1), total_num_frames);
                    
                    for ii = 1:size(vertex_values_sum1, 1)
                        vertex_values_sum_L(ii, :) = xcorr(vertex_values_sum1(ii, :), vertex_values_sum2(ii, :), 'coef');
                    end
                    
                    vertex_values_sum = vertex_values_sum_L;
                    vertex_values = vertex_values_sum_L;
                    
                    if laterality ~= 3
                        vertex_values_sum1 = squeeze(nansum(params.VertexVals1.('vertex_valuesR'), 2));
                        divisor1 = squeeze(sum(params.VertexVals1.('vertex_valuesR') ~= 0, 2));
                        divisor1(divisor1 == 0) = 1;
                        vertex_values_sum1 = vertex_values_sum1 ./ sqrt(divisor1);

                        vertex_values_sum2 = squeeze(nansum(params.VertexVals2.('vertex_valuesR'), 2));
                        divisor2 = squeeze(sum(params.VertexVals2.('vertex_valuesR') ~= 0, 2));
                        divisor2(divisor2 == 0) = 1;
                        vertex_values_sum2 = vertex_values_sum2 ./ sqrt(divisor2);

                        vertex_values_sum_R = NaN(size(vertex_values_sum1, 1), total_num_frames);

                        for ii = 1:size(vertex_values_sum1, 1)
                            vertex_values_sum_R(ii, :) = xcorr(vertex_values_sum1(ii, :), vertex_values_sum2(ii, :), 'coef');
                        end
                    end
                    
                else
                    vertex_values_sum = vertex_values_sum_R;
                    vertex_values = vertex_values_sum_R;
                end
                
            elseif strcmp(params.CompareMethod, 'SignProdSum')
                
                vertex_values_sum1 = squeeze(nansum(params.VertexVals1.(['vertex_values' side]), 2));
                divisor1 = squeeze(sum(params.VertexVals1.(['vertex_values' side]) ~= 0, 2));
                divisor1(divisor1 == 0) = 1;
                vertex_values_sum1_signs = sign(vertex_values_sum1);
                
                vertex_values_sum2 = squeeze(nansum(params.VertexVals2.(['vertex_values' side]), 2));
                divisor2 = squeeze(sum(params.VertexVals2.(['vertex_values' side]) ~= 0, 2));
                divisor2(divisor2 == 0) = 1;
                vertex_values_sum2_signs = sign(vertex_values_sum2);
                
                vertex_values_sum = (vertex_values_sum1_signs .* vertex_values_sum2_signs) .* ...
                    abs(vertex_values_sum1 + vertex_values_sum2) ./ sqrt(divisor1 + divisor2);
                vertex_values = vertex_values_sum;
                colormap(mycmap)
                caxis([-6 8])
            end
            
        end
        
        for frame_index = 1:total_num_frames
            
            %% overlay the colors
            
            disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status
            
            % sum the count of vertices across patients
            if (any(strcmp(overlay, {'CP', 'CNP', 'Subtraction'})) && ndims(vertex_values) == 3) || ...
                    any(strcmp(overlay, {'Electrode Density', 'Comparison'}))% the surface is colored according to zscores of power
                
                
                if any(electrodes_present(side_index,:) == 1)
                    if ~strcmp(overlay, 'Electrode Density') && ~(strcmp(overlay, 'Comparison') && ...
                            strcmp(params.CompareMethod, 'XCorr'))
                        % for now just overlay the color as number of patients overlapping in
                        % this area without any transparency
                        h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                            'FaceVertexCData',vertex_values_sum(:, frame_index),'FaceColor','interp',...
                            'LineStyle', 'none', 'FaceAlpha', 'flat');
                        alpha_data = ones(numel(h.FaceVertexCData), 1) * alpha;
                        alpha_data(abs(h.FaceVertexCData) <= zScoreThresh) = 0;
                        lightset = [0.6 0.5 0.1];
                        material(lightset);
                        set(h, 'FaceVertexAlphaData', alpha_data)
                    else
                        h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                            'FaceVertexCData',vertex_values_sum(:, frame_index),'FaceColor','interp',...
                            'LineStyle', 'none', 'FaceAlpha', alpha);
                        lightset = [0.6 0.5 0.1];
                        material(lightset);
                    end
                end
                
                %             elseif overlay == 0 % the surface is colored according to overlap of space around the electrodes
                %                 vertex_values_sum = NaN(cortex_v_num(side_index),1);
                %                 for i = 1:size(vertex_values,1)
                %                     vertex_values_sum(i,1) = sum(vertex_values(i,:));
                %                 end
                %                 colormap jet
                %                 caxis([-1 7])
            else
                lightset = [0.6 0.5 0.1];
                material(lightset);
            end
            
            %% Save each side of the brain in 3 views
            
            frames_set(side_index,frame_index,:) = savebrainimages(side_index,frame_index,views, figHandle);
            
            if any(electrodes_present(side_index,:) == 1) && exist('h','var')
                % get rid of the color patch so you can lay on the color patch for the
                % next frame
                delete(h)
            end
            
            % at the end of each frame generation save frames_set in case you
            % get kicked off and want to continue from before
            %         save('frames_set.mat','frames_set','-v7.3')
        end
    end
    
    if side_index == 1
        vertex_valuesL = vertex_values;
    elseif side_index == 2
        vertex_valuesR = vertex_values;
    end
end



if laterality == 0
    save('vertex_values.mat','vertex_valuesL','vertex_valuesR','-v7.3')
elseif laterality == 3 || laterality == 1
    save('vertex_values.mat','vertex_valuesL','-v7.3')
elseif laterality == 2
    save('vertex_values.mat','vertex_valuesR','-v7.3')
end

%% Put the 3 views into one image

if strcmp(overlay, 'Arrival Times')
    all_frames_num = electrode_mapping_times;
else
    all_frames_num = total_num_frames;
end

for frame_index = 1:all_frames_num
    
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    for side_index = start_side:sides_num
        if side_index == 1
            side = 'L';
        else
            side = 'R';
        end
        for view_index = 1:numel(ViewAngles)
            eval([side '_' ViewAngles{view_index} ' = [];'])
        end
    end
    
    %% name the frames for each side
    for side_index = start_side:sides_num
        % Determine what side we are working on
        if(side_index == 1)
            side = 'L';
            ViewAngles = {'Medial'};
        else
            side = 'R';
            ViewAngles = {'Lateral'};
        end
        
        if views == 4
            ViewAngles={'Lateral','Medial','Ventral','Posterior'};
        elseif views == 2
            ViewAngles={'Lateral','Medial'};
        elseif strcmp(views, 'Lateral')
            ViewAngles={'Lateral'};
        elseif strcmp(views, 'Medial')
            ViewAngles={'Medial'};
        end
        
        for view_index=1:length(ViewAngles)
            Frame = frames_set(side_index,frame_index,view_index).cdata;
            % Do cropping here
            Frame(all(all(Frame == 255, 3), 2), :, :) = [];
            Frame(:, all(all(Frame == 255, 3), 1), :) = [];
            eval([side '_' ViewAngles{view_index} ' = Frame;'])
        end
    end
    
    %% create the frames for each side
    %     if views == 4
    for side_index = start_side:sides_num
        if side_index == 1
            
            %             figure
            L_Lateral_Padded = [L_Lateral, ones(size(L_Lateral, 1),...
                max([0, size(L_Medial, 2) - size(L_Lateral, 2)]), 3) * 255];
            %             imshow(L_Lateral_Crop)
            
            L_Medial_Padded = [L_Medial, ones(size(L_Medial, 1),...
                max([0, size(L_Lateral, 2) - size(L_Medial, 2)]), 3) * 255];
            %             figure
            %             imshow(L_Medial_Crop)
            
            if laterality ~= 3
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(R_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];
            else
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(L_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];
            end
            %             figure
            %             imshow(L_Ventral_Crop)
            
            %             L_Ventral_padded = [ones(15,1090,3)*255; L_Ventral_Crop,ones(572,839,3)*255];
            %             L_Lateral_Medial = [L_Medial_Crop, L_Lateral_Crop];
            %             combined_views_left{frame_index} = [L_Lateral_Medial; L_Ventral_padded];
            L_Lateral_Medial = [L_Medial_Padded; L_Lateral_Padded];
            
            %             figure
            %             imshow(L_Lateral_Medial)
            
            if laterality == 0
                R_Posterior_Padded = [R_Posterior, ones(size(R_Posterior, 1),...
                    max([0, size(L_Ventral, 2) - size(R_Posterior, 2)]), 3) * 255];
                %               figure
                %               imshow(R_Posterior_Crop)
                
                L_Ventral_R_Posterior = [L_Ventral_Padded; R_Posterior_Padded];
                
                L_Lateral_Medial_Padded = [L_Lateral_Medial; ones(max([0, size(L_Ventral_R_Posterior, 1) - ...
                    size(L_Lateral_Medial, 1)]), size(L_Lateral_Medial, 2), 3) * 255];
                
                L_Ventral_R_Posterior_Padded = [L_Ventral_R_Posterior; ones(max([0, size(L_Lateral_Medial, 1) - ...
                    size(L_Ventral_R_Posterior, 1)]), size(L_Ventral_R_Posterior, 2), 3) * 255];
                %                     L_Ventral_R_Posterior_Padded = [L_Ventral_R_Posterior_Padded; ones(155, 545, 3)*255];
                %             figure
                %             imshow(L_Ventral_R_Posterior)
                combined_views_left{frame_index} = [L_Ventral_R_Posterior_Padded, L_Lateral_Medial_Padded];
            elseif laterality == 3
                L_Posterior_Padded = L_Posterior;
                
                L_Ventral_L_Posterior = [L_Ventral_Padded;L_Posterior_Padded];
                
                L_Lateral_Medial_Padded = [L_Lateral_Medial; ones(max([0, size(L_Ventral_L_Posterior, 1) - ...
                    size(L_Lateral_Medial, 1)]), size(L_Lateral_Medial, 2), 3) * 255];
                
                L_Ventral_L_Posterior_Padded = [L_Ventral_L_Posterior; ones(max([0, size(L_Lateral_Medial_Padded, 1) - ...
                    size(L_Ventral_L_Posterior, 1)]), size(L_Ventral_L_Posterior, 2), 3) * 255];
                
                %                     figure
                %                     imshow(L_Ventral_L_Posterior)
                combined_views_left{frame_index} = [L_Ventral_L_Posterior_Padded, L_Lateral_Medial_Padded];
            end
            
            figure
            imshow(combined_views_left{frame_index});
            set(gca,'position',[0 0 1 1],'units','normalized')
            
            % adding time stamp text to the figure
            if laterality == 3 && any(strcmp(overlay,{'CP', 'CNP', 'Subtraction', 'Comparison'}))
                if strcmp(overlay, 'Comparison') && strcmp(params.CompareMethod, 'XCorr')
                    time_txt = ['Lag = ' num2str(round((T2(frame_index))*1000)) 'ms'];
                else
                    time_txt = ['Time = ' num2str(round((T(frame_index)-1)*1000)) 'ms'];
                end
                uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                    'FontSize', 48, 'String', time_txt);
                
                text(size(combined_views_left{frame_index}, 2) - 1.5 * uc.Extent(3),...
                    size(combined_views_left{frame_index}, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 48)
            end
            
            fn = ['combined_viewsL' '_' num2str(frame_index)];
            %             export_fig(fn, '-tiff');
            print(fn, '-dtiff');
            close
        elseif side_index == 2
            
            R_Lateral_Padded = [R_Lateral, ones(size(R_Lateral, 1),...
                max([0, size(R_Medial, 2) - size(R_Lateral, 2)]), 3) * 255];
            %             figure
            %             imshow(R_Lateral_Crop)
            
            R_Medial_Padded = [R_Medial, ones(size(R_Medial, 1),...
                max([0, size(R_Lateral, 2) - size(R_Medial, 2)]), 3) * 255];
            %             figure
            %             imshow(R_Medial_Crop)
            R_Lateral_Medial = [R_Medial_Padded; R_Lateral_Padded];
            
            %             figure
            %             imshow(L_Lateral_Medial_padded)
            
            R_Ventral_Padded = [R_Ventral, ones(size(R_Ventral, 1),...
                max([0, size(L_Posterior, 2) - size(R_Ventral, 2)]), 3) * 255];
            %             figure
            %             imshow(R_Ventral_Crop)
            
            %             R_Ventral_padded = [ones(572,849,3)*255, R_Ventral_Crop];
            %             R_Lateral_Medial = [R_Lateral_Crop, R_Medial_Crop];
            %             combined_views_right{frame_index} = [R_Lateral_Medial; R_Ventral_padded];
            
            if laterality == 0
                L_Posterior_Padded = [L_Posterior, ones(size(L_Posterior, 1),...
                    max([0, size(R_Ventral, 2) - size(L_Posterior, 2)]), 3) * 255];
                %               figure
                %               imshow(L_Posterior_Crop)
                
                R_Ventral_L_Posterior = [R_Ventral_Padded; L_Posterior_Padded];
                %               figure
                %               imshow(R_Ventral_R_Posterior)
                
                R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(max([0, size(R_Ventral_L_Posterior, 1) - ...
                    size(R_Lateral_Medial, 1)]), size(R_Lateral_Medial, 2), 3) * 255];
                
                R_Ventral_L_Posterior_Padded = [R_Ventral_L_Posterior; ones(max([0, size(R_Lateral_Medial, 1) - ...
                    size(R_Ventral_L_Posterior, 1)]), size(R_Ventral_L_Posterior, 2), 3) * 255];
                
                combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_L_Posterior_Padded];
            elseif laterality == 2
                R_Posterior_Padded = R_Posterior; %(80:495,155:395,:);
                %               figure
                %               imshow(R_Posterior_Crop)
                
                R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(116,545,3)*255];
                
                R_Ventral_R_Posterior = [R_Ventral_Padded;R_Posterior_Padded];
                %               figure
                %               imshow(R_Ventral_R_Posterior)
                
                combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_R_Posterior];
                combined_views_right{frame_index} = [combined_views_right{frame_index}; 255*ones(15, 776, 3)];
            end
            
            figure
            imshow(combined_views_right{frame_index});
            set(gca,'position',[0 0 1 1],'units','normalized')
            fn = ['combined_viewsR' '_' num2str(frame_index)];
            %             export_fig(fn, '-tiff');
            saveas(gcf,[fn '.tif'])
            %             print(fn, '-dtiff');
            close
        end
    end
    %     elseif views == 2
    %         L_Medial_Padded = L_Medial; %(75:485,:,:);
    %         R_Lateral_Padded = R_Lateral; %(95:505,:,:);
    %     end
    
    %------------------------------------------------------------------
    %   Movie Assembly Procedure
    %------------------------------------------------------------------
    % Add all 6 components of each movie frame together onto one plot after
    % cropping. Then after assembling the frame, add to an avi movie buffer
    % output
    
    if laterality == 0
        figure
        %         if views == 4
        combined_views = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left{frame_index}, 1) - ...
            size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
            [combined_views_left{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
            size(combined_views_left{frame_index}, 1)]), size(combined_views_left{frame_index}, 2), 3) * 255]];
        %         elseif views == 2
        %             combined_views = [L_Medial_Padded R_Lateral_Padded];
        %         end
        imshow(combined_views);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis tight
        frame(frame_index,:,:,:) = combined_views;
        if exist('j','var')
            fn = ['combined_views_full_' num2str(frame_index) '_' num2str(j)];
        else
            fn = ['combined_views_full_' num2str(frame_index)];
        end
        
        if any(strcmp(overlay, {'CP', 'CNP', 'Subtraction', 'Comparison'}))
            if strcmp(overlay, 'Comparison') && strcmp(params.CompareMethod, 'XCorr')
                time_txt = ['Lag = ' num2str(round((T2(frame_index))*1000)) 'ms'];
            else
                time_txt = ['Time = ' num2str(round((T(frame_index)-1)*1000)) 'ms'];
            end
            if views == 4
                uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                    'FontSize', 48, 'String', time_txt);
                
                text(size(combined_views, 2) - 1.5 * uc.Extent(3),...
                    size(combined_views, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 48)
            else
                %                 uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                %                     'FontSize', 20, 'String', time_txt);
                %                 text(size(combined_views, 2) - 1.5 * uc.Extent(3),...
                %                     size(combined_views, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 20)
            end
        elseif strcmp(overlay, 'Arrival Times')
            T = [100 200 300 400 500 600 700 800 900 1000];
            time_txt = ['Time = ' num2str(T(frame_index)) ' ms'];
            text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
        end
        
        saveas(gcf,[fn '.tiff'])
        %         export_fig(fn, '-tiff');
        %         print(fn, '-dtiff');
        % clear combined_views_with_colorbar; clear combined_views_left; clear combined_colorbars;
        % clear Beta_L_Lateral_Crop; clear Beta_L_Medial_Crop; clear Beta_L_Ventral_Crop;
        % clear Beta_R_Lateral_Crop; clear Beta_R_Medial_Crop; clear Beta__Ventral_Crop;
        % clear combined_views;
        
        close
    end
end

% if you want to create a montage of evenly space post-stimulus frames,
% then uncomment this code

if views == 2
    %indices of frames to capture in montage
    frame_idx = [33 37 41 45 49 53 57];
    
    large_image = [];
    for i = 1:length(frame_idx)
        large_image = cat(1,squeeze(frame(frame_idx(i),:,:,:)),large_image);
    end
    figure
    imshow(large_image)
    set(gca,'position',[0 0 1 1],'units','normalized')
    saveas(gcf,'frames_montage.tiff')
    saveas(gcf,'frames_montage.eps')
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

function [electrode_vertex_values] = electrode_data_overlay(overlay,frame_index,all_electrodes, file_suffix, IsVisual)
overlay_frame_data = [];
if any(strcmp(overlay,  {'Electrode Distribution', 'Arrival Times', 'Electrode Density'}))
    %% Give every overlay value a 1
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
    % elseif overlay == 1
    %     %% First arrival times - need to add code for 2nd and 3rd------------------
    %     load('arrival_time.mat')%, 'half_max_times');
    %     newTimes = [];
    %     for i = 1:size(half_max_times,2)
    %         arrivalOne = half_max_times{1,i};
    %         if ~isempty(half_max_times{1,i})
    %             firstTime = arrivalOne(1);
    %         else
    %             firstTime = 0;
    %         end
    %         newTimes = [newTimes, firstTime];
    %     end
    %     newTimes = newTimes';
    %     overlay_frame_data = newTimes;
elseif any(strcmp(overlay, {'CP', 'CNP', 'Subtraction'}))
    %% getting power zscores from mean power of each electrode
    load(['meanpower_traces_57bins_zscore_recentered_rejoutliers', file_suffix, '.mat'])
    if strcmp(overlay, 'CP')
        trialtype = 1; % confirmed perceived
    elseif strcmp(overlay, 'CNP')
        trialtype = 2; % confirmed not perceived
    end
    if ~strcmp(overlay, 'Subtraction')
        overlay_frame_data = squeeze(meanpower_traces(:,trialtype,3,:));
    else
        overlay_frame_data = squeeze(meanpower_traces(:,1,3,:)) - squeeze(meanpower_traces(:,2,3,:));
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
%               within radius2 distance from each electrode's vertex
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

%This function should be combined with OverlayColor to speed up processing
%time...

function [vertexCdata]= CalculateTransparency(vertices_num,vertices, elecVertices, electrode_vertex_values, overlay)
maxAlpha    = 1;
radius1     = 1;
radius2     = 15;
electrodes_num      = length(elecVertices);
electrode_vector    = vertices(elecVertices,:);
vertexCdata = zeros(vertices_num, 1);
j=1;
k=1;

for vertIndex = 1 : vertices_num
    sharedTransparency = [];
    sharedElectrode = [];
    
    xCoord = ones(electrodes_num,1)*vertices(vertIndex, 1);
    yCoord = ones(electrodes_num,1)*vertices(vertIndex, 2);
    zCoord = ones(electrodes_num,1)*vertices(vertIndex, 3);
    
    % Calculate "Transparency" in other words the "weight" based on
    % distance from the vertex
    [distanceToElectrode, electrodeIndice] = sort(sqrt((electrode_vector(:,1) - xCoord).^2 + (electrode_vector(:,2) - yCoord).^2 + (electrode_vector(:,3) - zCoord).^2));
    
    if any(strcmp(overlay, {'Electrode Distribution', 'Electrode Density'}))
        if any(distanceToElectrode <= radius2)
            vertexCdata(vertIndex,1) = 1;
        end
    else
        
        for n = 1:length(electrodeIndice)
            if distanceToElectrode(n) <= radius1
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                sharedTransparency = [sharedTransparency; maxAlpha];
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
            elseif distanceToElectrode(n) <= radius2
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                sharedTransparency = [sharedTransparency; (maxAlpha - ((distanceToElectrode(n) - radius1)/(radius2-radius1))*maxAlpha)];
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
            end
        end
        % Grabs the zScore values associated with the electrodes that
        % contribute power to this vertex
        zIn = electrode_vertex_values(sharedElectrode)';
        weightedZ = []; %transparencySum = [];
        % Calculates the weight of the z score according to the transparencies and zScores associated with this vertex
        for h = 1:length(zIn)
            %                 weightedZ = [weightedZ zIn(h)*(sharedTransparency(h).^2)]
            weightedZ = [weightedZ zIn(h)*(sharedTransparency(h))];
            %                 transparencySum = [transparencySum sharedTransparency(h)];
        end
        
        
        %             indDelete = [];
        %             for y = 1:length(zIn)
        %                 if isnan(weightedZ(y))
        %                     indDelete = [indDelete y];
        %                 end
        %             end
        %             weightedZ(indDelete) = [];
        %             transparencySum(indDelete) = [];
        
        
        weightedZ = nansum(weightedZ);
        %             transparencySum = sum(transparencySum);
        %             weightedZ = weightedZ/transparencySum;
        % If there is no weighted zScore that means the vertex was greater than
        % 15mm from eevery electrode and has no power or zScore to display in
        % this frame
        if isempty(zIn)
            vertexCdata(vertIndex,1) = 0;
        elseif isnan(zIn)
            vertexCdata(vertIndex,1) = 0;
        else
            vertexCdata(vertIndex,1) = weightedZ;
        end
    end
end
end