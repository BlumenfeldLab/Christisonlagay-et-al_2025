function displayElectrodesInflated_v2_horizontal_display(data_location, save_location, patients,...
    overlay, flipX_montage, createMontage, laterality, color_map, color_map_inner_bounds,...
    color_map_outer_bounds, views, inflationstep, plot_electrodes, file_suffixes,...
    xtick_spacing_filename, display_color_bar, stim_offset, varargin)
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

    - save_location: Location of brainmap data parent folder to save to.
                     Folders will be created in this parent folder to store
                     the brainmaps.

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
                 'Subject Density'/8 = creates a density plot of subjects
    
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

    - color_map: Full filename for colormap to use for display.

    - color_map_inner_bound: Threshold for visualization of vertex values
                           on the brainmap. Vertices which express absolute
                           values below this threshold will not be
                           visualized.

    - color_map_outer_bounds: 2x1 numeric vector specifying the outer
                             boundaries to which the colormap will be
                             snapped.

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
    
    - plot_electrodes: boolean indicating whether or not to display the
                        electrodes as color-filled spheres on the vertex mesh
    
    - xtick_spacing_filename: string indicating proper time label filename.
                        these files are present in \\server2.med.yale.internal\data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\hnct_iceeg_scripting\eeg_analysis
    
    - display_color_bar: boolean indicating whether to append a colorbar to
                        the right of the combined images
    
    - flipX_montage: 1xn boolean array where n = length(patients)
                        indicating for which patients to use the flipX instead of the standard
                        montage. this is for bringing all patients into the 'same relative
                        stimulus space' as for the tactile task some patients are stimulated on
                        the left hand and some are on the right
    
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
    
    % path to all subject
    if ispc
        root = '//server2.med.yale.internal/Data25/';
    else
        root = '/mnt/Data25/';
    end
    
    %% Determine if overlay argument was properly passed
    overlay_labels = {'Electrode Distribution', 'Arrival Times', 'CP', 'CNP', 'Electrode Density', 'Subtraction', 'Comparison', 'Subject Density'};

    % If overlay numeric code, replace with the label
    if isnumeric(overlay)
        overlay = overlay_labels{overlay};
    end

    %% Define directory path to which images will be saved
    if laterality == 3; save_location = [save_location, '/flipX']; end
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
        case 'Subject Density'
            image_storage_folder = [save_location '/Subject_Density'];
    end

    %% Check more inputs
    valid_keywords = {'IsVisual', 'VertexVals1', 'VertexVals2', 'CompareMethod'};
    valid_compare_methods = {'XCorr', 'SignProdSum'};
    params.IsVisual = false(1, numel(patients));
    params.VertexVals1 = [];
    params.VertexVals2 = [];
    params.CompareMethod = '';
    if numel(varargin) > 1
        if exist(varargin{1}, 'dir')
            image_storage_folder = varargin{1};
            varargin(1) = [];
        end
    end
    if ~isfolder(image_storage_folder)
        mkdir(image_storage_folder)
    else
        disp('WARNING: image_storage_folder already exists. Ensure that you are not attempting to overwrite valuable data!')
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

    %% Setup for brainmap creation

    % Load the mni template and transform
    mni2fs_dir = [root, 'HNCT/icEEG Analysis/Analysis/EEG_behavior/Functions for Inflated Brain Display/mni2fs-master'];
    mni_nifti_path = [root, 'HNCT/icEEG Analysis/Analysis/EEG_behavior/Functions for Inflated Brain Display/', 'MNI_T1_1mm_stripped.nii'];
    load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
    mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
    Tmni = mnit1_1mm_stripped.transform;
    
    % Define patient electrode colors
    if ~strcmp(overlay, 'Comparison')
        addpath([root, 'HNCT/icEEG Analysis/Analysis/EEG_behavior/hnct_iceeg_scripting/eeg_analysis/'])
        colors = distinguishable_colors(length(patients), .5 * ones(1, 3));
    end

    % Define the face transparency value
%     alpha = 0.7;
      alpha = 0.8;

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

    % Define total number of frames for the present overlay setting
    if any(strcmp(overlay, {'Electrode Distribution', 'Arrival Times', 'Electrode Density', 'Parcellation', 'Subject Density'}))
        total_num_frames = 1;
    elseif any(strcmp(overlay, {'CP', 'CNP', 'Subtraction', 'Comparison'}))
        xtick_dir = [root, 'HNCT/icEEG Analysis/Analysis/EEG_behavior/hnct_iceeg_scripting/eeg_analysis/',xtick_spacing_filename];
        load(xtick_dir)
        if strcmp(overlay, 'Comparison') && strcmp(params.CompareMethod, 'XCorr')
            total_num_frames = 2 * numel(T) - 1;
            t_interval = T(2) - T(1);
            T2 = -(numel(T) - 1) * t_interval:t_interval:(numel(T) - 1) * t_interval;
        else
            total_num_frames = numel(T);
        end
    end

    %% Generating all frames for each side of the brain

    % Generate frames for L and R
    side_labels = {'L', 'R'};
    for side_index = start_side : sides_num

        % Determine what side we are working on
        if side_index == 1
            side = 'L';
            hem = 'lh';
        else
            side = 'R';
            hem = 'rh';
        end

        % Setup surface variables
        surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);  %%need directory
        fs_surf(side_index) = export(gifti(surf_fn));  %%we need here through line 472 for parcellation
        v_num = size(fs_surf(side_index).vertices, 1);
        fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] *...
            Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
        fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);

        surfrender_fn = fullfile(mni2fs_dir,['/surf/' hem '.inflated' num2str(inflationstep) '.surf.gii']);
        inflated_surf(side_index) = export(gifti(surfrender_fn));
        inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] *...
            Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni'; %'
        inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);

        % create a figure of the left or right side of the brain
        figHandle(side_index) = figure('Position', [70, 70, 700, 700]);

        % Initialize a temporary surface variable
        temp_surf = [];
        temp_surf.hem = hem; % choose the hemesphere 'lh' or 'rh'
        temp_surf.inflationstep = inflationstep; % 1 no inflation, 6 fully inflated
        temp_surf.decimation = 0;
        temp_surf = mni2fs_brain(temp_surf);
        set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
        set(gcf, 'color', [1 1 1]);
        hold on
        axis off;
        axis equal;

        % Initialize empty vertex_values variable
        if overlay ~= 8
            vertex_values = zeros(v_num, length(patients));
        end

        % if you are plotting arrival times, then each frame is a different set
        % of electrodes
        if strcmp(overlay, 'Arrival Times')
            electrode_mapping_times = 10;
        else
            electrode_mapping_times = 1;
        end

        for e = 1 : electrode_mapping_times

            scatter_count = 0;
            for p = 1 : length(patients)
            
                % Display current patient
                display(patients{p})
                 disp(data_location)
                % Define patient folder
                patient_folder = [data_location, patients{p}];
                % Define labels folder
                % labels_folder = [root,'HNCT/icEEG Analysis/Analysis/EEG_behavior/Group analysis/Parcellation_2022/',patients{p}];
                labels_folder = [root,'HNCT/icEEG Analysis/Analysis/EEG_behavior/',patients{p},'/icEEG/HNCT Tactile ID Combined'];
              
                % Create montage
                % if the left/right electrode montages have not been created, do so
                % now, indicated by createMontage = 0 or = 1
                if createMontage == 1 && ~strcmp(patients{p}, '436BP')
                    create_electrode_montage(patient_folder, labels_folder, laterality == 3, params.IsVisual(p))
                end

                % Now load some data files needed for electrode placement and surface shading
                cd([data_location, patients{p}]);
                if flipX_montage(p) && laterality == 0
                    load([side_labels{3 - side_index} '_MontageMap_flipX.mat']);
                    electrode = eval([side_labels{3 - side_index} '_MontageMap' ]);
                else
                    load([side '_MontageMap.mat']);
                    electrode = eval([side '_MontageMap' ]);
                end

                % load the X-flipped montage and add to list of electrodes
                if laterality == 3
                    load('R_MontageMap_flipX.mat');
                    electrode = cat(1,electrode,R_MontageMap);
                end

                % Discard depth electrodes or plot only some electrodes for arrival times
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
                    for l = 1 : length(labels_depth)
                        labels_depth_idx = [labels_depth_idx find(strcmp(labels_depth{l}, labels))];
                    end

                    % remove them from the electrode list that you will later access
                    [side_labels_depth, depth_idx, electrode_idx] = intersect(labels_depth_idx, electrode(:, 1));
                    electrode(electrode_idx, :) = [];

                    % if you are trying to plot arrival times, then you will need
                    % to plot only the electrodes in that time period
                    if strcmp(overlay, 'Arrival Times')
                        load('arrival_times_by_period_100ms.mat')
                        electrodes_this_period = find(arrival_times_by_period(:, e) == 1);
                        electrode = electrode(ismember(electrode(:, 1), electrodes_this_period), :);
                    end
                end

                % Plot patients electrodes
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
                    [elecVertices] = ProjectElectrode2TransSurf(v_num, fs_surf(side_index).vertices, all_electrodes, electrode);

                    % Plot the electrodes in a specified color and count them
                    % in case you need to remove them later
                    if plot_electrodes
                        scatter3(inflated_surf(side_index).vertices(elecVertices,1),...
                            inflated_surf(side_index).vertices(elecVertices,2),...
                            inflated_surf(side_index).vertices(elecVertices,3), 40, colorElectrode,'filled');
                    end
                    scatter_count = scatter_count + 1;
                end

                % cycle through all the frames to aggregate zscores for each vertex
                if any(strcmp(overlay, {'CP', 'CNP', 'Subtraction', 'Comparison'}))
                    for frame_index = 1 : total_num_frames
                        if isempty(electrode) % Some files only have electrodes on one hemisphere, hence the if statement
                            disp('-- No electrode data!');
                            FaceVertexAlphaData = zeros(v_num, 1);
                            vertexCdata = zeros(v_num, 1);
                        else
                            display([patients{p} ' Aggregating Frame ' num2str(frame_index)])
                            if ~any(strcmp(overlay, {'Electrode Density', 'Comparison', 'Subject Density'}))
                                [electrode_vertex_values] = electrode_data_overlay(overlay, frame_index, all_electrodes, file_suffixes{p}, params.IsVisual(p), total_num_frames);
                                [vertexCdata] = CalculateTransparency(v_num, inflated_surf(side_index).vertices, elecVertices, electrode_vertex_values, overlay);
                                vertex_values(:,p,frame_index) = vertexCdata;
                            end
                        end
                    end
                elseif strcmp(overlay, 'Electrode Density')
                    if isempty(electrode) % Some files only have electrodes on one hemisphere, hence the if statement
                        disp('-- No electrode data!');
                        FaceVertexAlphaData = zeros(v_num, 1);
                        vertexCdata = zeros(v_num, 1);
                    else
                        [electrode_vertex_values] = electrode_data_overlay(overlay, 1, all_electrodes, file_suffixes{p}, params.IsVisual(p), total_num_frames);
                        electrode_vertex_values(isnan(electrode_vertex_values)) = 0;
                        electrode_vertex_values(electrode_vertex_values ~= 0) = 1;
                        [vertexCdata] = CalculateTransparency(v_num, inflated_surf(side_index).vertices, elecVertices, electrode_vertex_values, overlay);
                        vertex_values(:, p) = vertexCdata; % we will load parcellation map here
                    end
                elseif strcmp(overlay, 'Subject Density')
                    if isempty(electrode) % Some files only have electrodes on one hemisphere, hence the if statement
                        disp('-- No electrode data!');
                        FaceVertexAlphaData = zeros(v_num, 1);
                        vertexCdata = zeros(v_num, 1);
                    else
                        [electrode_vertex_values] = electrode_data_overlay(overlay, 1, all_electrodes, file_suffixes{p}, params.IsVisual(p), total_num_frames);
                        electrode_vertex_values(isnan(electrode_vertex_values)) = 0;
                        electrode_vertex_values(electrode_vertex_values ~= 0) = 1;
                        [vertexCdata] = CalculateTransparency(v_num, inflated_surf(side_index).vertices, elecVertices, electrode_vertex_values, overlay);
                        vertex_values(:, p) = vertexCdata; % we will load parcellation map here
                    end
                end
            end

            % Frame Generation if saving without overlaying color onto the brain
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

            % go to the correct folder
            cd(image_storage_folder)
            load(color_map);

            if ~any(strcmp(overlay, {'Electrode Density', 'Comparison', 'Subject Density'}))
                vertex_values_sum = squeeze(nansum(vertex_values, 2));
                divisor = squeeze(sum(vertex_values ~= 0, 2));
                divisor(divisor == 0) = 1;

                % Take weighted mean across patients
                vertex_values_sum = vertex_values_sum ./ sqrt(divisor);

                % load the colormap and set axis of colormap
                colormap(mycmap)
                caxis(color_map_outer_bounds);

            elseif any(strcmp(overlay, {'Electrode Density', 'Subject Density'}))

                vertex_values_sum = sum(vertex_values, 2);
                colormap jet
                %caxis([0, max(vertex_values_sum)])
                caxis(color_map_outer_bounds)
                
                %old
                %{
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
                %}
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
                    caxis(color_map_outer_bounds)
                    %caxis([0, max(vertex_values_sum)])
                end

            end

            for frame_index = 1 : total_num_frames

                %% overlay the colors
                disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status

                % sum the count of vertices across patients
                if (any(strcmp(overlay, {'CP', 'CNP', 'Subtraction'})) && ndims(vertex_values) == 3) || ...
                        any(strcmp(overlay, {'Electrode Density', 'Comparison', 'Subject Density'}))% the surface is colored according to zscores of power

                    if any(electrodes_present(side_index, :) == 1)
                        if ~(strcmp(overlay, 'Comparison') && strcmp(params.CompareMethod, 'XCorr'))
                            % for now just overlay the color as number of patients overlapping in
                            % this area without any transparency
                            h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                                'FaceVertexCData',vertex_values_sum(:, frame_index),'FaceColor','interp',...
                                'LineStyle', 'none', 'FaceAlpha', 'flat');
                            if ~any(strcmp(overlay, {'Electrode Density', 'Subject Density'}))
                                alpha_data = ones(numel(h.FaceVertexCData), 1) * alpha;
                                alpha_data(h.FaceVertexCData >= color_map_inner_bounds(1) & h.FaceVertexCData <= color_map_inner_bounds(2)) = 0;
                            else
                                alpha_data = ones(numel(h.FaceVertexCData), 1);
                            end
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
                else
                    lightset = [0.6 0.5 0.1];
                    material(lightset);
                end

                %% Save each side of the brain in 3 views
                frames_set(side_index,frame_index, :) = savebrainimages(side_index, frame_index, views, figHandle);
                if any(electrodes_present(side_index,:) == 1) && exist('h','var')
                    % get rid of the color patch so you can lay on the color patch for the
                    % next frame
                    delete(h)
                end
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
    
    % begin bulk replacement
    combined_views_left = cell(all_frames_num, 1);
    combined_views_right = cell(all_frames_num, 1);
%     for frame_index = 32 : 34
    for frame_index = 1 : all_frames_num
        ViewAngles={'Lateral','Medial','Ventral','Posterior'};
        for side_index = start_side:sides_num
            if side_index == 1
                side = 'L';
            else
                side = 'R';
            end
            for view_index = 1 : numel(ViewAngles)
                eval([side '_' ViewAngles{view_index} ' = [];'])
            end
        end

        %% name the frames for each side
        for side_index = start_side : sides_num

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
                disp(['Pre-crop Frame size for ', side, '_', ViewAngles{view_index}, ': ', num2str(size(Frame))]);
                % Do cropping here
                Frame(all(all(Frame == 255, 3), 2), :, :) = [];
                disp(['Mid-crop Frame size for ', side, '_', ViewAngles{view_index}, ': ', num2str(size(Frame))]);
                Frame(:, all(all(Frame == 255, 3), 1), :) = [];
                disp(['Post-crop Frame size for ', side, '_', ViewAngles{view_index}, ': ', num2str(size(Frame))]);
                eval([side '_' ViewAngles{view_index} ' = Frame;'])
            end
        end
        
        %% create the frames for each side
        
        % horz_edits: begin bulk replacement to create the frames for each side
        for side_index = start_side : sides_num
            if side_index == 1
                
                L_Ventral = cat(1, L_Ventral, ones(5, size(L_Ventral, 2), 3) * 255);
                L_Ventral = cat(1, ones(5, size(L_Ventral, 2), 3) * 255, L_Ventral);
                
                L_Lateral = cat(1, L_Lateral, ones(floor((size(L_Ventral, 1) - size(L_Lateral, 1)) / 2), size(L_Lateral, 2), 3) * 255);
                L_Lateral = cat(1, ones(size(L_Ventral, 1) - size(L_Lateral, 1), size(L_Lateral, 2), 3) * 255, L_Lateral);
                
                L_Medial = cat(1, L_Medial, ones(floor((size(L_Ventral, 1) - size(L_Medial, 1))/ 2), size(L_Medial, 2), 3) * 255);
                L_Medial = cat(1, ones(size(L_Ventral, 1) - size(L_Medial, 1), size(L_Medial, 2), 3) * 255, L_Medial);

                L_Posterior = cat(1, L_Posterior, ones(floor((size(L_Ventral, 1) - size(L_Posterior, 1))/ 2), size(L_Posterior, 2), 3) * 255);
                L_Posterior = cat(1, ones(size(L_Ventral, 1) - size(L_Posterior, 1), size(L_Posterior, 2), 3) * 255, L_Posterior);

                %combined_views_left{frame_index} = cat(2, L_Lateral, L_Medial, L_Ventral, L_Posterior);
                combined_views_left{frame_index} = cat(2, L_Ventral, L_Lateral, L_Medial, L_Posterior);
                disp(['pre-colorbar left_views size at frame ', num2str(frame_index), ': ', num2str(size(combined_views_left{frame_index}))]);
                
                % add a colorbar
                if display_color_bar
                    color_bar_image = imread('L_1_Ventral.tif');
                    color_bar_image = color_bar_image(:, 850 : 1024, :);
                    vert_dim = size(combined_views_left{frame_index}, 1);
                    horz_dim = floor((vert_dim / size(color_bar_image, 1)) * size(color_bar_image, 2));
                    for r = 1 : 3
                        color_bar_image_resize(:, :, r) = imresize(color_bar_image(:, :, r), [vert_dim, horz_dim]);
                    end
                    combined_views_left{frame_index} = [combined_views_left{frame_index}, color_bar_image_resize];
                    disp(['post-colorbar left_views size at frame ', num2str(frame_index), ': ', num2str(size(combined_views_left{frame_index}))]);
                end
                
                figure
                imshow(combined_views_left{frame_index});
                set(gca, 'position', [0 0 1 1], 'units', 'normalized')

                % adding time stamp text to the figure
                if laterality == 3 && any(strcmp(overlay,{'CP', 'CNP', 'Subtraction', 'Comparison'}))
                    if strcmp(overlay, 'Comparison') && strcmp(params.CompareMethod, 'XCorr')
                        time_txt = ['Lag = ' num2str(round((T2(frame_index)) * 1000)) 'ms'];
                    else
                        t_num = round((T(frame_index) - stim_offset) * 1000);
                        if strcmp(xtick_spacing_filename, 'spectrogram_times_57bins.mat')
                            time_txt = [num2str(t_num - 125) ' to ' num2str(t_num + 125) ' ms'];
                        elseif strcmp(xtick_spacing_filename, 'spectrogram_times_40bins.mat')
                            time_txt = [num2str(t_num - 25) ' to ' num2str(t_num + 25) ' ms'];
                        elseif strcmp(xtick_spacing_filename, 'filter_times_win100ms_shift50ms.mat') ||...
                                strcmp(xtick_spacing_filename, 'filter_times_win100ms_shift25ms.mat')
                            time_txt = [num2str(t_num - 50) ' to ' num2str(t_num + 50) ' ms'];
                        end
                    end
                    uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                        'FontSize', 32, 'String', time_txt);
                    text(size(combined_views_left{frame_index}, 2) - uc.Extent(3) - 150,...
                        size(combined_views_left{frame_index}, 1) - uc.Extent(4) - 10, time_txt, 'Fontsize', 12)
                end
                
                % save figure
                saveas(gcf,['combined_viewsL' '_' num2str(frame_index), '.fig'])
                saveas(gcf,['combined_viewsL' '_' num2str(frame_index), '.tif'])
                close
                
            elseif side_index == 2
                
                R_Ventral = cat(1, R_Ventral, ones(5, size(R_Ventral, 2), 3) * 255);
                R_Ventral = cat(1, ones(5, size(R_Ventral, 2), 3) * 255, R_Ventral);

                R_Lateral = cat(1, R_Lateral, ones(floor((size(R_Ventral, 1) - size(R_Lateral, 1)) / 2), size(R_Lateral, 2), 3) * 255);
                R_Lateral = cat(1, ones(size(R_Ventral, 1) - size(R_Lateral, 1), size(R_Lateral, 2), 3) * 255, R_Lateral);

                R_Medial = cat(1, R_Medial, ones(floor((size(R_Ventral, 1) - size(R_Medial, 1))/ 2), size(R_Medial, 2), 3) * 255);
                R_Medial = cat(1, ones(size(R_Ventral, 1) - size(R_Medial, 1), size(R_Medial, 2), 3) * 255, R_Medial);

                R_Posterior = cat(1, R_Posterior, ones(floor((size(R_Ventral, 1) - size(R_Posterior, 1))/ 2), size(R_Posterior, 2), 3) * 255);
                R_Posterior = cat(1, ones(size(R_Ventral, 1) - size(R_Posterior, 1), size(R_Posterior, 2), 3) * 255, R_Posterior);
                
                %combined_views_right{frame_index} = cat(2, R_Posterior, R_Ventral, R_Medial, R_Lateral);
                combined_views_right{frame_index} = cat(2, R_Posterior, R_Medial, R_Lateral, R_Ventral);
                disp(['right_views size at frame ', num2str(frame_index), ': ', num2str(size(combined_views_right{frame_index}))]);
                
                figure
                imshow(combined_views_right{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                saveas(gcf,['combined_viewsR' '_' num2str(frame_index), '.fig'])
                saveas(gcf,['combined_viewsR' '_' num2str(frame_index), '.tif'])
            end
        end
        % horz_edits: end bulk replacement
        
        %------------------------------------------------------------------
        %   Movie Assembly Procedure
        %------------------------------------------------------------------
        % Add all 6 components of each movie frame together onto one plot after
        % cropping. Then after assembling the frame, add to an avi movie buffer
        % output
        fprintf('size(combined_views_right) is %s\n', mat2str(size(combined_views_right{frame_index})))

        
        fprintf('size(combined_views_left) is %s\n', mat2str(size(combined_views_left{frame_index})))
        if laterality == 0
            figure
            % horz_edits: replace combined_views assigment with the following
            combined_views = cat(2, combined_views_right{frame_index}, combined_views_left{frame_index});

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
                    t_num = round((T(frame_index) - stim_offset) * 1000);
                    if strcmp(xtick_spacing_filename, 'spectrogram_times_57bins.mat')
                        time_txt = [num2str(t_num - 125) ' to ' num2str(t_num + 125) ' ms'];
                    elseif strcmp(xtick_spacing_filename, 'spectrogram_times_40bins.mat')
                        time_txt = [num2str(t_num - 25) ' to ' num2str(t_num + 25) ' ms'];
                    elseif strcmp(xtick_spacing_filename, 'filter_times_win100ms_shift50ms.mat') ||...
                            strcmp(xtick_spacing_filename, 'filter_times_win100ms_shift25ms.mat') ||...
                            strcmp(xtick_spacing_filename, 'filter_times_win100ms_shift25ms_157bins.mat')
                        time_txt = [num2str(t_num - 50) ' to ' num2str(t_num + 50) ' ms'];
                    end
                end
                if views == 4
                    uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                        'FontSize', 32, 'String', time_txt);

                    text(size(combined_views, 2) - uc.Extent(3) - 150,...
                        size(combined_views, 1) - uc.Extent(4) - 10, time_txt, 'Fontsize', 12)
                end
            elseif strcmp(overlay, 'Arrival Times')
                T = [100 200 300 400 500 600 700 800 900 1000];
                time_txt = ['Time = ' num2str(T(frame_index)) ' ms'];
                text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
            end
            saveas(gcf, [fn '.fig'])
            saveas(gcf, [fn '.tiff'])            
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
            large_image = cat(1,squeeze(frame(frame_idx(i), :, :, :)),large_image);
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
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral', 'Medial'};
        Light_Pos = [-1200 -300 700; 1000 -1100 -200];
        View_Pos = [-90 0;90, 0];
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
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral'};
        Light_Pos = [-1200 -300 700];
        View_Pos = [-90 0];
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
    else
        Light_Pos=0;
        View_Pos=0;
    end
end

for i = 1:length(ViewAngles)
    lightHandle = light('Position', Light_Pos(i,:));
    view(View_Pos(i, :));
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

    % before you save this picture, add a colorbar
    if frame_index == 1
%         colorbar 
        colorbar('FontSize',22);
        print(gcf, fn, '-dtiff');
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
for elec_index = 1 : size(electrode, 1)
    xCoord = ones(vertices_num,1)*electrode(elec_index,2);
    yCoord = ones(vertices_num,1)*electrode(elec_index,3);
    zCoord = ones(vertices_num,1)*electrode(elec_index,4);
    [~, minInd] = min(sqrt((vertices(:,1) - xCoord).^2 + (vertices(:,2) - yCoord).^2 + (vertices(:,3) - zCoord).^2));
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

function create_electrode_montage(patient_folder, labels_folder, flipX, IsVisual)
% move to the patients folder %Julia Ding 11/11/2020
cd(labels_folder)

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

clear b1 b2 montage a1 a2 map

b1 = linspace(1, length(labels), length(labels))';
b1 = num2cell(b1);
b2 = labels(:,1);
montage = [b1 b2];

cd(patient_folder) 

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
% Original code
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

% Modified code section, updated on 10/4/2023 by Taruna Yadav. It was
% needed as Montage labels were not created correctly by the original code
% section.
% try
%     for i=1:size(map,1)
%         str             =   map{i,1};
%         str(isspace(str)) = '';
%         str_delim       =   regexp(str,'_','split');
%         firstLetter     =   str_delim{1};
%         positionLetter  =   strtrim(str_delim{1,1}(1));
%         digit_split     =   split(str, ';');
%         digitStr        =   digit_split{end};
%         replaceStr      =   strcat(firstLetter, digitStr);
%         map{i,1}        =   replaceStr; % i.e. "A12"
%         map{i,5}        =   positionLetter; % Either L or R
%     end
% catch ex
%     disp(['Error reading Map file! There was a problem on line ' num2str(i) '. Verify that file is in the correct format (i.e each row is Letter_Side*Number. Unacceptable entries: J_Sup_Parietal8,PEG1)']);
%     throw(ex);
% end
% Combine montage with 3-d electrode locations

Montage = montage;
for i = 1 : size(montage, 1)
    MontageMap(i, 1) = montage{i, 1}; % Abbreviated Channel ("A12")
    for j = 1 : size(map, 1)
        if strcmpi(strtrim(montage{i, 2}), map{j, 1})
            % ....Commented out on 10/3/2023 by taruna
%             if i == 113
%                 leah = [];
%             end
            %.....
            MontageMap(i, 2) = map{j, 2}; % X Coordinates
            MontageMap(i, 3) = map{j, 3}; % Y Coordinates
            MontageMap(i, 4) = map{j, 4}; % Z Coordinates
            Position(i, 1) = map{j, 5}; % Side (R or L)
        end
    end
end

% check MontageMap to see if there are any errors that result in
% (0,0,0) coordinates for (x,y,z)
for i = 1:size(MontageMap,2)
    if MontageMap(i, 1) == 0 && MontageMap(i, 2) == 0 && MontageMap(i, 3) == 0
        display(MontageMap{i,1});
        error('This electrode name in Map or Map_flipX is not matching, double check the original Map sheet for misnamed electrodes compared to labels');
    end
end

% Break the montage into L and R sides
j = 1;
k = 1;

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

function [electrode_vertex_values] = electrode_data_overlay(overlay,frame_index,all_electrodes, file_suffix, IsVisual, total_num_frames)
overlay_frame_data = [];

if any(strcmp(overlay,  {'Electrode Distribution', 'Arrival Times', 'Electrode Density', 'Subject Density'}))

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

    % getting power zscores from mean power of each electrode
    
    load(['meanpower_traces_', num2str(total_num_frames), 'bins_zscore_recentered_rejoutliers', file_suffix, '.mat'])
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

    if any(strcmp(overlay, {'Electrode Distribution', 'Electrode Density', 'Subject Density'}))
        if any(distanceToElectrode <= radius2)
            vertexCdata(vertIndex, 1) = sum(distanceToElectrode <= radius2);
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
if strcmp(overlay, 'Subject Density')
    vertexCdata(vertexCdata ~= 0) = maxAlpha;
end
end
