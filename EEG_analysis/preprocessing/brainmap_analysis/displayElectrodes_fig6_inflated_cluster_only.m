
% save location
save_location = 'E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 77 Bins Sounds Restricted by CP Accuracy 100ms\filter_from2to2_baseline1to0r\kmeans_clust\CP_w250ms_test';

% vertex values
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 77 Bins Sounds Restricted by CP Accuracy 100ms\filter_from2to2_baseline1to0r\Gamma_power_CP\vertex_values.mat')

% previous clustering
cluster_file = '/mnt/Data9/HNCT/auditory/results/31patients_77bins_100ms_filter_from2to2_baseline1to0/CP_clusters_kmedoids_w250ms/clustering_tmp_workspace.mat';

% preferences
cluster_method = 'kmeans'; % options: 'kmeans', 'kmedoids'
min_clusters = 2;
max_clusters = 6;
laterality = 0;
xtick_spacing_filename = 'filter_times_win100ms_shift25ms.mat';
stimulus_offset = 1;
use_prev_cluster = 0;

if ~isfolder(save_location)
  mkdir(save_location)
end

% load prior clusters
if use_prev_cluster
    load(cluster_file, 'changing_vertices_idx', 'constant_vertices_idx', 'vertices_by_cluster_reduced', 'vertices_by_cluster_reduced_0_250ms', 'vertices_by_cluster_reduced_250_500ms',...
        'vertices_by_cluster_reduced_500_750ms', 'vertices_by_cluster_reduced_750_1000ms')
    changing_vertices_idx_prev = changing_vertices_idx;
    constant_vertices_idx_prev = constant_vertices_idx;
    vertices_by_cluster_reduced_prev = vertices_by_cluster_reduced;
    vertices_by_cluster_reduced_0_250ms_prev = vertices_by_cluster_reduced_0_250ms;
    vertices_by_cluster_reduced_250_500ms_prev = vertices_by_cluster_reduced_250_500ms;
    vertices_by_cluster_reduced_500_750ms_prev = vertices_by_cluster_reduced_500_750ms;
    vertices_by_cluster_reduced_750_1000ms_prev = vertices_by_cluster_reduced_750_1000ms;
end

% load time and get stimulus bin
load(['E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\hnct_iceeg_scripting\eeg_analysis\', xtick_spacing_filename])
total_num_bins = length(T);

% Load relevant patient processing metadata
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_20_06_09.mat')
patients = patients(1 : end - 1); % remove patient 32

num = squeeze(sum(vertex_valuesL ~= 0, 2));
num(num == 0) = 1;

% Determine if we want to divide by sqrt(num patients) (Wendys method) per
% vertex or to divide by num patients
wendys_method = 1;
if wendys_method
    if laterality == 0
        vertex_values_sum_L = squeeze(nansum(vertex_valuesL, 2)) ./ sqrt(num);
        num = squeeze(sum(vertex_valuesR ~=0, 2));
        num(num == 0) = 1;
        vertex_values_sum_R = squeeze(nansum(vertex_valuesR, 2)) ./ sqrt(num);
    elseif laterality == 3
        vertex_values_sum = squeeze(nansum(vertex_valuesL, 2)) ./ sqrt(num);
    end
    y_limits = [-4 10];
else
    if laterality == 0
        vertex_values_sum_L = squeeze(nansum(vertex_valuesL, 2)) ./ num;
        num = squeeze(sum(vertex_valuesR ~=0, 2));
        num(num == 0) = 1;
        vertex_values_sum_R = squeeze(nansum(vertex_valuesR, 2)) ./ num;
    elseif laterality == 3
        vertex_values_sum = squeeze(nansum(vertex_valuesL, 2)) ./ num;
    end
    y_limits = [-.5 1];
end

% final colors for k-means using confirmed not perceived trials:
% dark blue DMN
% red visual assocation
% green visual cortex
% darker cyan fronto-parietal
ColorSetMaps = [
    1.00 1.00 1.00 % no cluster - white
    0.00 1.00 0.00 % Data 1 - green
    0.00 0.00 1.00 % Data 2 - blue
    1.00 0.00 0.00 % Data 3 - red
    0.00 1.00 1.00 % Data 4 - cyan
    1.00 0.00 1.00 % Data 5 - magenta
    0.75 0.75 0.00 % Data 6 - RGB
    0.25 0.25 0.25 % Data 7
    0.75 0.25 0.25 % Data 8
    0.95 0.95 0.00 % Data 9
    0.25 0.25 0.75 % Data 10
    ];
ColorSet = [
    0.00 1.00 0.00 % Data 1 - green
    0.00 0.00 1.00 % Data 2 - blue
    1.00 0.00 0.00 % Data 3 - red
    0.00 1.00 1.00 % Data 4 - cyan
    1.00 0.00 1.00 % Data 5 - magenta
    0.75 0.75 0.00 % Data 6 - RGB
    0.25 0.25 0.25 % Data 7
    0.75 0.25 0.25 % Data 8
    0.95 0.95 0.00 % Data 9
    0.25 0.25 0.75 % Data 10
    ];

%% previously displayElectrodes_phase3_fig6_inflated

cd(save_location)

% get stimulus bin
[~, stimulus_bin] = min(abs(T - stimulus_offset));

% combine left and right side vertex time series
if laterality == 0
    vertex_values_sum_bilateral = cat(1,vertex_values_sum_L,vertex_values_sum_R);
    poststim = vertex_values_sum_bilateral(:, stimulus_bin : end);
elseif laterality == 3
    poststim = vertex_values_sum(:, stimulus_bin : end);
end

% exclude vertices with 0 std dev
poststim_std = std(poststim, 0, 2);
if ~use_prev_cluster
    constant_vertices_idx = poststim_std == 0;
    changing_vertices_idx = poststim_std > 0;
else
    constant_vertices_idx = constant_vertices_idx_prev;
    changing_vertices_idx = changing_vertices_idx_prev;
end
vertices_by_cluster = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_0_250ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_250_500ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_500_750ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_750_1000ms = zeros(numel(changing_vertices_idx), max_clusters);

% take only the vertices that are not constant
poststim_reduced = poststim(changing_vertices_idx, :);
if laterality == 0
    vertex_values_sum_bilateral_reduced = vertex_values_sum_bilateral(changing_vertices_idx, :);
    vertex_values_sum_reduced = vertex_values_sum_bilateral_reduced;
else
    vertex_values_sum_reduced = vertex_values_sum(changing_vertices_idx, :);
end

% redo the above for the first 200 ms (4 frames) post-stimulus
time_res = mean(diff(T));
poststim_reduced_0_250ms = poststim_reduced(:, 1 : floor(.25 / time_res));
poststim_reduced_250_500ms = poststim_reduced(:, floor(.25 / time_res) : floor(.5 / time_res));
poststim_reduced_500_750ms = poststim_reduced(:, floor(.5 / time_res) : floor(.75 / time_res));
poststim_reduced_750_1000ms = poststim_reduced(:, floor(.75 / time_res) : end);

% Intialize vertices x clusters NaN matrices
vertices_by_cluster_reduced = NaN(sum(changing_vertices_idx) ,max_clusters);
vertices_by_cluster_reduced_0_250ms = NaN(sum(changing_vertices_idx), max_clusters);
vertices_by_cluster_reduced_250_500ms = NaN(sum(changing_vertices_idx), max_clusters);
vertices_by_cluster_reduced_500_750ms = NaN(sum(changing_vertices_idx), max_clusters);
vertices_by_cluster_reduced_750_1000ms = NaN(sum(changing_vertices_idx), max_clusters);

% Initialize vars for sillouette values, average sillouette values and p values per cluster
s = NaN(sum(changing_vertices_idx), max_clusters);
s_mean = NaN(max_clusters);
p_val = NaN(max_clusters);
s_0_250ms = NaN(sum(changing_vertices_idx), max_clusters);
s_0_250ms_mean = NaN(max_clusters);
p_val_0_250ms = NaN(max_clusters);
s_250_500ms = NaN(sum(changing_vertices_idx),max_clusters);
s_250_500ms_mean = NaN(max_clusters);
p_val_250_500ms = NaN(max_clusters);
s_500_750ms = NaN(sum(changing_vertices_idx),max_clusters);
s_500_750ms_mean = NaN(max_clusters);
p_val_500_750ms = NaN(max_clusters);
s_750_1000ms = NaN(sum(changing_vertices_idx),max_clusters);
s_750_1000ms_mean = NaN(max_clusters);
p_val_750_1000ms = NaN(max_clusters);

% Initialize vars for something...
meantrace = NaN(max_clusters, max_clusters, total_num_bins);
baseline_trace = NaN(max_clusters, stimulus_bin - 1);
poststim_trace = NaN(max_clusters, total_num_bins - stimulus_bin + 1);
meantrace_0_250ms = NaN(max_clusters, max_clusters, total_num_bins);
baseline_trace_0_250ms = NaN(max_clusters, stimulus_bin - 1);
poststim_trace_0_250ms = NaN(max_clusters, total_num_bins - stimulus_bin + 1);
meantrace_250_500ms = NaN(max_clusters, max_clusters, total_num_bins);
baseline_trace_250_500ms = NaN(max_clusters, stimulus_bin - 1);
poststim_trace_250_500ms = NaN(max_clusters, total_num_bins - stimulus_bin + 1);
meantrace_500_750ms = NaN(max_clusters, max_clusters, total_num_bins);
baseline_trace_500_750ms = NaN(max_clusters, stimulus_bin - 1);
poststim_trace_500_750ms = NaN(max_clusters, total_num_bins - stimulus_bin + 1);
meantrace_750_1000ms = NaN(max_clusters, max_clusters, total_num_bins);
baseline_trace_750_1000ms = NaN(max_clusters, stimulus_bin - 1);
poststim_trace_750_1000ms = NaN(max_clusters, total_num_bins - stimulus_bin + 1);

load(xtick_spacing_filename)
for i = min_clusters : max_clusters

    fprintf('%d clusters\n', i)

    % postimulus kmeans clustering
    if use_prev_cluster
        vertices_by_cluster_reduced(:, i) = vertices_by_cluster_reduced_prev(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            vertices_by_cluster_reduced(:,i) = kmeans(poststim_reduced, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            vertices_by_cluster_reduced(:,i) = kmedoids(poststim_reduced, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % silhouette values
    fprintf('\tCalculating silhouette values for %d clusters poststimulus\n', i)
    figure;
    s(:, i) = silhouette(poststim_reduced, vertices_by_cluster_reduced(:, i), 'correlation');
    for j = 1 : i
        cluster_vertex_idx = find(vertices_by_cluster_reduced(:, i) == j);
        s_mean(i, j) = mean(s(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_' num2str(i) '.tif']);

    % cluster time series
    figure
    legend_str = cell(1, i);
    for j = 1 : i
        % find the indices of vertices both left and right that are part of that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced(:,i) == j);
        vals = vertex_values_sum_reduced(cluster_vertex_idx, :);
        plot(T - stimulus_offset, nanmean(vals, 1), 'Color', ColorSet(j, :), 'Linewidth', 2); hold on;
        meantrace(i, j, :) = nanmean(vals, 1);
        baseline_trace(j, :) = abs(meantrace(i, j, 1 : stimulus_bin - 1));
        poststim_trace(j, :) = abs(meantrace(i, j, stimulus_bin : total_num_bins));
        p_val(i, j) = ranksum(baseline_trace(j, :), poststim_trace(j, :));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val(i, j)) ', s = ' num2str(s_mean(i, j))];
    end
    title(['Gamma Power Zscore, # Clusters = ' num2str(i)], 'Fontsize', 14);
    xlabel('Time (s)', 'FontSize', 14)
    ylabel('Power z-score', 'FontSize', 14)
    xlim([T(1) - stimulus_offset T(end) - stimulus_offset]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL, 'Color', 'k');
    line(xL, [0 0], 'Color', 'k');
    legend(legend_str,'Location', 'Northwest')
    legend('boxoff');
    set(gca, 'Fontsize', 14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf, ['Gamma_clusters_' num2str(i) '.tif']);
    clearvars legend_str

    % 0 - 250 ms post stimulus kmeans
    if use_prev_cluster
        vertices_by_cluster_reduced_0_250ms(:, i) = vertices_by_cluster_reduced_0_250ms_prev(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            [vertices_by_cluster_reduced_0_250ms(:,i), ~] = kmeans(poststim_reduced_0_250ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            [vertices_by_cluster_reduced_0_250ms(:,i), ~] = kmedoids(poststim_reduced_0_250ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % Graph silhouette values
    fprintf('\tCalculating silhouette values for %d clusters 0-250ms\n', i)
    figure
    [s_0_250ms(:, i), ~] = silhouette(poststim_reduced_0_250ms, vertices_by_cluster_reduced_0_250ms(:, i), 'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_0_250ms(:, i) == j);
        s_0_250ms_mean(i, j) = mean(s_0_250ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_0_250ms_' num2str(i) '.tif']);

    % Gamma part earlier in the code because of graphing reasons
    for j = 1 : i
        cluster_vertex_idx = find(vertices_by_cluster_reduced_0_250ms(:, i) == j);
        % plot that index
        plot(T - stimulus_offset, nanmean(vertex_values_sum_reduced(cluster_vertex_idx, :), 1), 'Color', ColorSet(j, :)); hold on

        meantrace_0_250ms(i, j, :) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx, :), 1);
        baseline_trace_0_250ms(j, :) = abs(meantrace_0_250ms(i, j, 1 : stimulus_bin - 1));
        poststim_trace_0_250ms(j, :) = abs(meantrace_0_250ms(i, j, stimulus_bin : total_num_bins));
        [p_val_0_250ms(i, j)] = ranksum(baseline_trace_0_250ms(j, :), poststim_trace_0_250ms(j, :));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_0_250ms(i,j)) ', s = ' num2str(s_0_250ms_mean(i,j))];
    end
    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest','AutoUpdate', 'off')
    legend('boxoff');
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T - stimulus_offset, squeeze(meantrace(i, j, :)), 'Color', [.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_0_250ms_' num2str(i) '.tif']);
    clearvars legend_str

    % 250 - 500 ms kmeans
    if use_prev_cluster
        vertices_by_cluster_reduced_250_500ms(:, i) = vertices_by_cluster_reduced_250_500ms_prev(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            [vertices_by_cluster_reduced_250_500ms(:,i), ~] = kmeans(poststim_reduced_250_500ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            [vertices_by_cluster_reduced_250_500ms(:,i), ~] = kmedoids(poststim_reduced_250_500ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % Graph silhouette values
    figure
    fprintf('\tCalculating silhouette values for %d clusters 250-500ms\n', i)
    [s_250_500ms(:, i), ~] = silhouette(poststim_reduced, vertices_by_cluster_reduced_250_500ms(:, i), 'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_250_500ms(:, i) == j);
        s_250_500ms_mean(i, j) = mean(s_250_500ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_400ms_' num2str(i) '.tif']);

    % Gamma part earlier in the code because of graphing reasons
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_250_500ms(:, i) == j);
        % plot that index
        plot(T-1,nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1), 'Color', ColorSet(j,:)); hold on;
        % do stats on the prestim and poststim gamma power zscore changes
        meantrace_250_500ms(i,j,:) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1);
        baseline_trace_250_500ms(j,:) = abs(meantrace_250_500ms(i,j,1:stimulus_bin - 1));
        poststim_trace_250_500ms(j,:) = abs(meantrace_250_500ms(i,j,stimulus_bin:total_num_bins));
        [p_val_250_500ms(i,j)] = ranksum(baseline_trace_250_500ms(j,:),poststim_trace_250_500ms(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_250_500ms(i,j)) ', s = ' num2str(s_250_500ms_mean(i,j))];
    end

    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest', 'AutoUpdate', 'off')
    legend('boxoff');
    % set the size of the legend and axes fonts
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T-1,squeeze(meantrace(i,j,:)), 'Color',[.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_250_500ms_' num2str(i) '.tif']);
    clearvars legend_str

    % 500 - 750 ms kmeans
    if use_prev_cluster
        vertices_by_cluster_reduced_500_750ms(:, i) = vertices_by_cluster_reduced_500_750ms_prev(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            [vertices_by_cluster_reduced_500_750ms(:,i), ~] = kmeans(poststim_reduced_500_750ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            [vertices_by_cluster_reduced_500_750ms(:,i), ~] = kmedoids(poststim_reduced_500_750ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % Graph silhouette values
    figure
    fprintf('\tCalculating silhouette values for %d clusters 500-750ms\n', i)
    [s_500_750ms(:, i), ~] = silhouette(poststim_reduced,vertices_by_cluster_reduced_500_750ms(:, i), 'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_500_750ms(:, i) == j);
        s_500_750ms_mean(i, j) = mean(s_500_750ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_400ms_' num2str(i) '.tif']);

    % Gamma part earlier in the code because of graphing reasons
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_500_750ms(:, i) == j);
        % plot that index
        plot(T-1,nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1), 'Color', ColorSet(j,:)); hold on;
        % do stats on the prestim and poststim gamma power zscore changes
        meantrace_500_750ms(i,j,:) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1);
        baseline_trace_500_750ms(j,:) = abs(meantrace_500_750ms(i,j,1:stimulus_bin - 1));
        poststim_trace_500_750ms(j,:) = abs(meantrace_500_750ms(i,j,stimulus_bin:total_num_bins));
        [p_val_500_750ms(i,j)] = ranksum(baseline_trace_500_750ms(j,:),poststim_trace_500_750ms(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_500_750ms(i,j)) ', s = ' num2str(s_500_750ms_mean(i,j))];
    end

    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest', 'AutoUpdate', 'off')
    legend('boxoff');
    % set the size of the legend and axes fonts
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T-1,squeeze(meantrace(i,j,:)), 'Color',[.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_500_750ms_' num2str(i) '.tif']);
    clearvars legend_str

    % 750 - 1000 ms kmeans
    if use_prev_cluster
        vertices_by_cluster_reduced_750_1000ms(:, i) = vertices_by_cluster_reduced_750_1000ms_prev(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            [vertices_by_cluster_reduced_750_1000ms(:,i), ~] = kmeans(poststim_reduced_750_1000ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            [vertices_by_cluster_reduced_750_1000ms(:,i), ~] = kmedoids(poststim_reduced_750_1000ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % Graph silhouette values
    figure
    fprintf('\tCalculating silhouette values for %d clusters 750-1000ms\n', i)
    [s_750_1000ms(:, i), ~] = silhouette(poststim_reduced,vertices_by_cluster_reduced_750_1000ms(:, i), 'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_750_1000ms(:, i) == j);
        s_750_1000ms_mean(i, j) = mean(s_750_1000ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_400ms_' num2str(i) '.tif']);

    % Gamma part earlier in the code because of graphing reasons
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_750_1000ms(:, i) == j);
        % plot that index
        plot(T-1,nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1), 'Color', ColorSet(j,:)); hold on;
        % do stats on the prestim and poststim gamma power zscore changes
        meantrace_750_1000ms(i,j,:) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1);
        baseline_trace_750_1000ms(j,:) = abs(meantrace_750_1000ms(i,j,1:stimulus_bin - 1));
        poststim_trace_750_1000ms(j,:) = abs(meantrace_750_1000ms(i,j,stimulus_bin:total_num_bins));
        [p_val_750_1000ms(i,j)] = ranksum(baseline_trace_750_1000ms(j,:),poststim_trace_750_1000ms(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_750_1000ms(i,j)) ', s = ' num2str(s_750_1000ms_mean(i,j))];
    end

    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest', 'AutoUpdate', 'off')
    legend('boxoff');
    % set the size of the legend and axes fonts
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T-1,squeeze(meantrace(i,j,:)), 'Color',[.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_750_1000ms_' num2str(i) '.tif']);
    clearvars legend_str
end

% put back the constant vertices so they are in the same order again
vertices_by_cluster(changing_vertices_idx, :) = vertices_by_cluster_reduced;
vertices_by_cluster_0_250ms(changing_vertices_idx, :) = vertices_by_cluster_reduced_0_250ms;
vertices_by_cluster_250_500ms(changing_vertices_idx, :) = vertices_by_cluster_reduced_250_500ms;
vertices_by_cluster_500_750ms(changing_vertices_idx, :) = vertices_by_cluster_reduced_500_750ms;
vertices_by_cluster_750_1000ms(changing_vertices_idx, :) = vertices_by_cluster_reduced_750_1000ms;

% separate vertices into left and right hemisphere
if laterality == 0
    vertices_by_cluster_L = vertices_by_cluster(1 : size(vertex_values_sum_L, 1), :);
    vertices_by_cluster_R = vertices_by_cluster(size(vertex_values_sum_L, 1) + 1 : end, :);

    vertices_by_cluster_L_0_250ms = vertices_by_cluster_0_250ms(1 : size(vertex_values_sum_L, 1), :);
    vertices_by_cluster_R_0_250ms = vertices_by_cluster_0_250ms(size(vertex_values_sum_L, 1) + 1 : end, :);

    vertices_by_cluster_L_250_500ms = vertices_by_cluster_250_500ms(1 : size(vertex_values_sum_L, 1), :);
    vertices_by_cluster_R_250_500ms = vertices_by_cluster_250_500ms(size(vertex_values_sum_L, 1) + 1 : end, :);

    vertices_by_cluster_L_500_750ms = vertices_by_cluster_500_750ms(1 : size(vertex_values_sum_L, 1), :);
    vertices_by_cluster_R_500_750ms = vertices_by_cluster_500_750ms(size(vertex_values_sum_L, 1) + 1 : end, :);

    vertices_by_cluster_L_750_1000ms = vertices_by_cluster_750_1000ms(1 : size(vertex_values_sum_L, 1), :);
    vertices_by_cluster_R_750_1000ms = vertices_by_cluster_750_1000ms(size(vertex_values_sum_L, 1) + 1 : end, :);
end

%save('clustering_tmp_workspace.mat', '-v7.3')
