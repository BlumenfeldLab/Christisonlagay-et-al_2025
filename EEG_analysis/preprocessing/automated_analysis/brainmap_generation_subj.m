
% specify the data location
data_location = 'V:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\';

% specify the pipeline info file to use
load('V:\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_20_04_27.mat')

% specify subject indices to use
% subj_indices = [3 12 13 30 19];
subj_indices = [3];

for i = subj_indices
    
    % save location
    save_location = [data_location, patients{i}, '\Auditory ID 77 Bins Sounds Restricted by CP Accuracy 100ms'];
    
    % bounds determined first to include 95% of all data then adjusted s.t. 0
    % lies at a null |heat| value
    color_map_outer_bounds = [-6 8];
    color_map_inner_bounds = [-2, 2];
    color_map = 'V:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\cmap6to8WRX.mat';
    plot_electrodes = 1;
    xtick_spacing_filename = 'filter_times_win100ms_shift25ms.mat';
    display_color_bar = 1;
    
    % generic options
    views = 4;
    inflationstep = 5;
    
    % CP bilateral
    laterality = 0;
    createMontage = 1;
    overlay = 8;
    flipX_montage=boolean([0]);
    flipX_montage
    displayElectrodesInflated_v2(data_location, save_location, patients(i), overlay,...
        createMontage, flipX_montage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
        views, inflationstep, plot_electrodes, meanpower_traces_suffixes_77bins_100ms_filter_from2to2(i), xtick_spacing_filename, 1)
    
    
 
    
    % CNP bilateral
    createMontage = 0;
    overlay = 4;
    displayElectrodesInflated_v2(data_location, save_location, patients(i), overlay,...
        createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
        views, inflationstep, plot_electrodes, meanpower_traces_suffixes_77bins_100ms_filter_from2to2(i), xtick_spacing_filename, 1)
    
    %{
    % Subtraction bilateral
    overlay = 6;
    displayElectrodesInflated_v2(data_location, save_location, patients(i), overlay,...
        createMontage, laterality, color_map, color_map_inner_bound, color_map_outer_bounds,...
        views, inflationstep, plot_electrodes,  meanpower_traces_suffixes(i))

    % CP flipX
    createMontage = 1;
    laterality = 3;
    overlay = 3;
    displayElectrodesInflated_v2(data_location, save_location, patients(i), overlay,...
        createMontage, laterality, color_map, color_map_inner_bound, color_map_outer_bounds,...
        views, inflationstep, plot_electrodes, meanpower_traces_suffixes(i))
    
    % CNP flipX
    createMontage = 0;
    overlay = 4;
    displayElectrodesInflated_v2(data_location, save_location, patients(i), overlay,...
        createMontage, laterality, color_map, color_map_inner_bound, color_map_outer_bounds,...
        views, inflationstep, plot_electrodes, meanpower_traces_suffixes(i))
    
    % Subtraction flipX
    overlay = 6;
    displayElectrodesInflated_v2(data_location, save_location, patients(i), overlay,...
        createMontage, laterality, color_map, color_map_inner_bound, color_map_outer_bounds,...
        views, inflationstep, plot_electrodes, meanpower_traces_suffixes(i))
    %}
end
