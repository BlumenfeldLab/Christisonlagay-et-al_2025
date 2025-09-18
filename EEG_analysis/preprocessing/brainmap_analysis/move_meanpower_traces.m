

% define source and destination locations for file movement
source = '\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\';
destination = '\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\';

% 'load pipeline info
load('\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_21_04_15.mat')

for i = 1 : 31
    fprintf('\n%d / %d', i, 31)
    source_file = [source, patients{i}, '\icEEG\HNCT Auditory ID Combined\meanpower_traces_77bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_from2to2_baseline5to0_common_b.mat'];

    % remove old file
    %delete([destination, patients{i}, '\meanpower_traces_79bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_from2to2_common.mat'])

    % rename file
    %destination_file = [source, patients{i}, '\icEEG\HNCT Auditory ID Combined\meanpower_traces_40bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_50ms_norej_filter.mat'];
    %movefile(source_file, destination_file)

    try
        % move file
        destination_folder = [destination, patients{i}];
        copyfile(source_file, destination_folder)
    catch
        fprintf(sprintf('No file for %s\n', patients{i}))
    end
end
