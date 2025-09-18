
% define source and destination locations for file movement
if ispc
    source = '//172.23.254.106/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/';
    destination = '//172.23.254.106/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/HNCT ERP movies/';
    load('//172.23.254.106/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/normal_pipeline_file_info/normal_pipeline_file_info_21_04_22.mat')
else
    source = '/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/';
    destination = '/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/HNCT ERP movies/';
    load('/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/normal_pipeline_file_info_tactile/normal_pipeline_file_info_21_04_22.mat')
end

for i = 1 : 6
    fprintf('\n%d / 31', i)
    try
        try
            load([source, patients{i}, '/icEEG/HNCT Auditory ID Combined/meanpower_traces_157bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_2to2_baseline5to0_b.mat']);
        catch
            load([source, patients{i}, '/icEEG/HNCT Auditory ID Combined/meanpower_traces_156bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_2to2_baseline5to0_b.mat']);
        end

        % extract subset

        % 100/50 ms
        %meanpower_traces = meanpower_traces(:, :, :, 21 : 59);

        % 100/25 ms 77 bins
        meanpower_traces = meanpower_traces(:, :, :, 41 : 117);

        % 250/31.25 ms 57 bins
        %meanpower_traces = meanpower_traces(:, :, :, 33 : 89);

        % save
        save([source, patients{i}, '/icEEG/HNCT Auditory ID Combined/meanpower_traces_77bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_from2to2_baseline5to0_b.mat'], 'meanpower_traces')
    catch
        fprintf('\t no data for %s', patients{i})
    end
end
