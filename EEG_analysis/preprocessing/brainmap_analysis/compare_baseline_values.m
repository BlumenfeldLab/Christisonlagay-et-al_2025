root = '\\172.23.254.106\Data25\';

% define the pipeline info file
load([root, 'HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info_tactile\normal_pipeline_file_info_20_11_11.mat'])
data_folder = [root, 'HNCT\icEEG Analysis\Analysis\EEG_behavior\'];

%% z-score baseline -1 to 0 s
for p = 3
    
    % define save directory
    save_dir = [data_folder, patients{p}, '\icEEG\HNCT Tactile ID Combined\Electrode Group Filter 100ms Window 50ms Shift From 2to2 Baseline 5to0\'];
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    
    % load labels and baseline timeseries
    load([root, 'HNCT\icEEG Analysis\Analysis\EEG_behavior\', patients{p}, '\icEEG\HNCT Tactile ID Combined\labels_all_gray_matter.mat'])
    load([root, 'HNCT\icEEG Analysis\Analysis\EEG_behavior\', patients{p}, '\icEEG\HNCT Tactile ID Combined\baselines157bins_sounds_restricted_CPonly_100ms_filter_2to2_baseline5to0.mat'])
    load([root, 'HNCT\icEEG Analysis\Analysis\EEG_behavior\', patients{p}, '\icEEG\HNCT Tactile ID Combined\meanpower_traces_77bins_zscore_recentered_rejoutliers_sounds_restricted_CPonly_100ms_filter_from2to2_baseline5to0.mat'])
    
    % isolate gamma
    baselines = squeeze(baselines(:, :, 3, :));
    meanpower_traces = squeeze(meanpower_traces(:, :, 3, :));
    
    % collect same electrode traces
    group_traces = [];
    group_baselines = [];
    baseline_ranges = [];
    group_labels{1} = labels{1}(isstrprop(labels{1}, 'alpha'));
    i = 1;
    group_indices{i} = [];
    j = 1;
    for k = 1 : length(labels)
        electrode = labels{k}(isstrprop(labels{k}, 'alpha'));
        if ~strcmp(electrode, group_labels{i})
            group_labels{i + 1} = electrode;
            i = i + 1;
            j = 1;
            group_indices{i} = [];
        end
        group_indices{i}(end + 1) = k;
        group_traces{i}(j, :, :) = meanpower_traces(k, :, :);
        group_baselines{i}(j, :, :) = cat(3, baselines(k, :, :), NaN(1, 2, 57));
        j = j + 1;
        baseline_ranges(k, 1, :) = [min(baselines(k, 1, :), [], 'all'), max(baselines(k, 1, :), [], 'all')];
        baseline_ranges(k, 2, :) = [min(baselines(k, 2, :), [], 'all'), max(baselines(k, 2, :), [], 'all')];
    end
    
    % display
    t = linspace(-.95, .95, 77);
    baseline_ranges_group = zeros(length(group_labels), 2, 2);
    for i = 1 : length(group_labels)
        
        % compute baseline range
        baseline_ranges_group(i, 1, :) = [min(group_baselines{i}(:, 1, :),[], 'all'), max(group_baselines{i}(:, 1, :),[], 'all')];
        baseline_ranges_group(i, 2, :) = [min(group_baselines{i}(:, 2, :),[], 'all'), max(group_baselines{i}(:, 2, :),[], 'all')];
        
        % display timeseries
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2, 2, 1)
        plot(t, squeeze(group_baselines{i}(:, 1, :))', 'linewidth', 1.5)
        title(sprintf('range = [%.2f, %.2f]', baseline_ranges_group(i, 1, 1), baseline_ranges_group(i, 2, 2)), 'fontsize', 14)
        xlim([-1, 1])
        xlabel('time (ms)')
        ylabel('gamma power (uV^2)')
        subplot(2, 2, 2)
        plot(t, squeeze(group_baselines{i}(:, 2, :))', 'linewidth', 1.5)
        title(sprintf('range = [%.2f, %.2f]', baseline_ranges_group(i, 2, 1), baseline_ranges_group(i, 2, 2)), 'fontsize', 14)
        xlim([-1, 1])
        xlabel('time (ms)')
        ylabel('gamma power (uV^2)')
        subplot(2, 2, 3)
        plot(t, squeeze(group_traces{i}(:, 1, :))', 'linewidth', 1.5); hold on
        plot(t, zeros(size(group_traces{i}, 3), 1), 'k--', 'linewidth', 1.5); hold off
        xlim([-1, 1])
        xlabel('time (ms)')
        ylabel('gamma power (uV^2)')
        subplot(2, 2, 4)
        plot(t, squeeze(group_traces{i}(:, 2, :))', 'linewidth', 1.5); hold on
        plot(t, zeros(size(group_traces{i}, 3), 1), 'k--', 'linewidth', 1.5); hold off
        xlim([-1, 1])
        xlabel('time (ms)')
        ylabel('gamma power (uV^2)')
        suptitle(group_labels{i})
        
        % save figure
        saveas(gcf, [save_dir, group_labels{i}, '.tiff'])
        close all
        
        % compare baseline ranges across all electrodes
        figure('units','normalized','outerposition',[0 0 1 1])
        subplot(2, 1, 1)
        plot(1 : length(group_indices{i}), baseline_ranges(group_indices{i}, 1, 1), 'linewidth', 1.5); hold on
        plot(1 : length(group_indices{i}), baseline_ranges(group_indices{i}, 1, 2), 'linewidth', 1.5); hold off
        xlim([0 length(group_indices{i}) + 1])
        set(gca, 'xtick', 1 : length(group_indices{i}), 'xticklabels', labels(group_indices{i}), 'fontsize', 12)
        title('CP baseline range across electrodes')
        ylabel('gamma power (uV^2)')
        subplot(2, 1, 2)
        plot(1 : length(group_indices{i}), baseline_ranges(group_indices{i}, 2, 1), 'linewidth', 1.5); hold on
        plot(1 : length(group_indices{i}), baseline_ranges(group_indices{i}, 2, 2), 'linewidth', 1.5); hold off
        xlim([0 length(group_indices{i}) + 1])
        set(gca, 'xtick', 1 : length(group_indices{i}), 'xticklabels', labels(group_indices{i}), 'fontsize', 12)
        title('CNP baseline range across electrodes')
        ylabel('gamma power (uV^2)')
        
        % save figure
        saveas(gcf, [save_dir, 'electrode_range_', group_labels{i}, '.tiff'])
        close all
    end
    
    % save data
    save('electrode_groups.mat', 'group_labels', 'group_indices')
    save('baselines157bins_sounds_restricted_CPonly_100ms_filter_2to2_baseline5to0_ranges.mat', 'baseline_ranges')
end


