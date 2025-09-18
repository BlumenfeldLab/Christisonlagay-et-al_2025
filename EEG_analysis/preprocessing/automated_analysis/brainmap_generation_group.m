%% setup
datadrive='V:';
% specify data location
data_location = [datadrive ']\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\'];

% specify the pipeline info to use
load([ datadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_21_04_22.mat'])
patients = patients(1 : 31);
meanpower_traces_suffixes = mean_suffixes_157bins_100ms_filter_baseline5to0_common_b(1 : 31);

% specify the location to save to
save_location='D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\subjectdensity';
% save_location = [data_location, 'Group\Inflated Surface\31 Patients Auditory ID 157 Bins Sounds Restricted by CP Accuracy 100ms\filter_baseline5to0_common_b\Subject Density'];

% bounds determined first to include 95% of all data then adjusted s.t. 0
% lies at a null |heat| value
%color_map_outer_bounds = [-3.435 6.325];
%color_map_inner_bounds = [-2.215, 2.36];
color_map_outer_bounds = [-15 15];
color_map_inner_bounds = [-2, 2];
%color_map = '\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\cmap6to8WRX.mat';
color_map = '\\Mwmj04vg66\d\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\extendedcolorbar.mat';


 

% generic settings
createMontage = 0;
views = 4;
inflationstep = 5;
plot_electrodes = 1;
xtick_spacing_filename = 'filter_times_win100ms_shift25ms_157bins.mat';
display_color_bar = 1;

% initialize the flipX_montage (flips patients electrodes for whom stimulation was on the left hand) as 0s (no flips)
flipX_montage = zeros(length(patients), 1);

%% standard brainmaps

% bilateral settings
laterality = 0;

%{
% electrode density
overlay = 5;
displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
    createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
    views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1)

% electrode map only
overlay = 1;
displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
    createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
    views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1)

% subject density
overlay = 8;
displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
    createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
    views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1)

% CP bilateral
overlay = 3;
displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
    createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
    views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1, 2)
%}

overlay = 8;

if overlay==8
    color_map = '\\Mwmj04vg66\d\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\modifiedjet2.mat';
end
displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
    createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
    views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1, 2)
% 
% 
% % CNP bilateral
% overlay = 4;
% displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
%     createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
%     views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1, 2)
% 
% % Subtraction bilateral
% overlay = 6;
% displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
%     createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
%     views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1, 2)
%     
% % flipX settings
% laterality = 3;
% 
% % CP flipX
% overlay = 3;
% displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
%     createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
%     views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1, 2)
% 
% % CNP flipX
% overlay = 4;
% displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
%     createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
%     views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1, 2)
% 
% % Subtraction flipX
% overlay = 6;
% displayElectrodesInflated_v2(data_location, save_location, patients, overlay, flipX_montage,...
%     createMontage, laterality, color_map, color_map_inner_bounds, color_map_outer_bounds,...
%     views, inflationstep, plot_electrodes, meanpower_traces_suffixes, xtick_spacing_filename, 1, 2)

%% brainmaps from vertex values
%{
% load vertex values if using
laterality = 0;
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 77 Bins Sounds Restricted by CP Accuracy 100ms\filter_from2to2_baseline1to0_bar3to7\Gamma_power_CP\vertex_values.mat')
overlay = 3;
displayElectrodesInflated_v2_from_vertex_values(data_location, save_location, patients,...
    overlay, laterality, color_map, color_map_inner_bounds,...
    color_map_outer_bounds, views, inflationstep, plot_electrodes,...
    xtick_spacing_filename, vertex_valuesL, vertex_valuesR, display_color_bar)

% load vertex values if using
load('\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 77 Bins Sounds Restricted by CP Accuracy 100ms\filter_from2to2_baseline5to0_common_b\Gamma_power_CNP\vertex_values.mat')
overlay = 4;
displayElectrodesInflated_v2_from_vertex_values(data_location, save_location, patients,...
    overlay, laterality, color_map, color_map_inner_bounds,...
    color_map_outer_bounds, views, inflationstep, plot_electrodes,...
    xtick_spacing_filename, vertex_valuesL, vertex_valuesR, display_color_bar)
%}

%% parcellation brainmaps from vertex values
%{
% colormap
nparcels = 100;
ROIcolors = distinguishable_colors(nparcels);
color_map = colormap(ROIcolors);

% load vertex values if using
laterality = 0;
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\parcellation_maps\mask_fs_surf_visualize_inflation5\verticesIndex_L_whitematter.mat')
vertexvaluesL = ones(175409, 1);
vertexvaluesL(whiteMatter_vertices) = 2;
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\parcellation_maps\mask_fs_surf_visualize_inflation5\verticesIndex_R_whitematter.mat')
vertexvaluesR = ones(171090, 1);
vertexvaluesR(whiteMatter_vertices) = 2;

% save location
save_location = 'E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\parcellation_maps\mask_fs_surf_visualize_inflation5\mask';

overlay = 8;
displayElectrodesInflated_v2_parcellation(data_location, save_location, {'311KH'}, 8, 1, 0, color_map,...
        0, 0, 4, inflationstep, 0, {'_sounds_restricted_CPonly'}, vertexvaluesL, vertexvaluesR, nparcels)
%}