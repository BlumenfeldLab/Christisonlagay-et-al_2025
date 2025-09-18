
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
cd(save_location)

% load prior clusters
if use_prev_cluster
    load(cluster_file, 'changing_vertices_idx', 'vertices_by_cluster_reduced', 'vertices_by_cluster_reduced_0_250ms', 'vertices_by_cluster_reduced_250_500ms',...
        'vertices_by_cluster_reduced_500_750ms', 'vertices_by_cluster_reduced_750_1000ms')
    changing_vertices_idx_prev = changing_vertices_idx;
    vertices_by_cluster_reduced_prev = vertices_by_cluster_reduced;
    vertices_by_cluster_reduced_0_250ms_prev = vertices_by_cluster_reduced_0_250ms;
    vertices_by_cluster_reduced_250_500ms_prev = vertices_by_cluster_reduced_250_500ms;
    vertices_by_cluster_reduced_500_750ms_prev = vertices_by_cluster_reduced_500_750ms;
    vertices_by_cluster_reduced_750_1000ms_prev = vertices_by_cluster_reduced_750_1000ms;
end

% load time and get stimulus bin
load(['E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\hnct_iceeg_scripting\eeg_analysis\', xtick_spacing_filename])
total_num_bins = length(T);
[~, stimulus_bin] = min(abs(T - stimulus_offset));

% compute vertex vaues and extract post stimulus interval
num = squeeze(sum(vertex_valuesL ~= 0, 2));
num(num == 0) = 1;
if laterality == 0
    vertex_values_sum_L = squeeze(nansum(vertex_valuesL, 2)) ./ sqrt(num);
    num = squeeze(sum(vertex_valuesR ~=0, 2));
    num(num == 0) = 1;
    vertex_values_sum_R = squeeze(nansum(vertex_valuesR, 2)) ./ sqrt(num);
    vertex_values_sum = cat(1,vertex_values_sum_L,vertex_values_sum_R);
    poststim = vertex_values_sum(:, stimulus_bin : end);
elseif laterality == 3
    vertex_values_sum = squeeze(nansum(vertex_valuesL, 2)) ./ sqrt(num);
    poststim = vertex_values_sum(:, stimulus_bin : end);
end

% exclude vertices with 0 std dev
poststim_std = std(poststim, 0, 2);
if ~use_prev_cluster
    changing_vertices_idx = poststim_std > 0;
else
    changing_vertices_idx = changing_vertices_idx_prev;
end
vertex_values_sum_reduced = vertex_values_sum(changing_vertices_idx, :);
poststim_reduced = poststim(changing_vertices_idx, :);
time_res = mean(diff(T));
poststim_reduced_0_250ms = poststim_reduced(:, 1 : floor(.25 / time_res));
poststim_reduced_250_500ms = poststim_reduced(:, floor(.25 / time_res) : floor(.5 / time_res));
poststim_reduced_500_750ms = poststim_reduced(:, floor(.5 / time_res) : floor(.75 / time_res));
poststim_reduced_750_1000ms = poststim_reduced(:, floor(.75 / time_res) : end);

% intialize vertices x clusters NaN matrices
vertices_by_cluster_reduced = NaN(sum(changing_vertices_idx) ,max_clusters);
vertices_by_cluster_reduced_0_250ms = NaN(sum(changing_vertices_idx), max_clusters);
vertices_by_cluster_reduced_250_500ms = NaN(sum(changing_vertices_idx), max_clusters);
vertices_by_cluster_reduced_500_750ms = NaN(sum(changing_vertices_idx), max_clusters);
vertices_by_cluster_reduced_750_1000ms = NaN(sum(changing_vertices_idx), max_clusters);

% inialize sillouette values, average sillouette values and cluster p-value
s = NaN(sum(changing_vertices_idx), max_clusters);
s_mean = NaN(max_clusters);
p_val = NaN(max_clusters);
s_0_250ms = NaN(sum(changing_vertices_idx), max_clusters);
s_mean_0_250ms = NaN(max_clusters);
p_val_0_250ms = NaN(max_clusters);
s_250_500ms = NaN(sum(changing_vertices_idx),max_clusters);
s_250_500ms_mean = NaN(max_clusters);
p_val_250_500ms = NaN(max_clusters);
s_500_750ms = NaN(sum(changing_vertices_idx),max_clusters);
s_mean_500_750ms = NaN(max_clusters);
p_val_500_750ms = NaN(max_clusters);
s_750_1000ms = NaN(sum(changing_vertices_idx),max_clusters);
s_mean_750_1000ms = NaN(max_clusters);
p_val_750_1000ms = NaN(max_clusters);

% initialize cluster mean timeseries
meantrace = NaN(max_clusters, max_clusters, total_num_bins);
meantrace_0_250ms = NaN(max_clusters, max_clusters, total_num_bins);
meantrace_250_500ms = NaN(max_clusters, max_clusters, total_num_bins);
meantrace_500_750ms = NaN(max_clusters, max_clusters, total_num_bins);
meantrace_750_1000ms = NaN(max_clusters, max_clusters, total_num_bins);

ColorSet = [0.00, 1.00, 0.00;
            0.00, 0.00, 1.00;
            1.00, 0.00, 0.00;
            0.00, 1.00, 1.00;
            1.00, 0.00, 1.00;
            0.75, 0.75, 0.00;
            0.25, 0.25, 0.25;
            0.75, 0.25, 0.25;
            0.95, 0.95, 0.00;
            0.25, 0.25, 0.75];
ColorSetMaps = ColorSet;
y_limits = [-6 10];
for i = min_clusters : max_clusters

    fprintf('%d clusters\n', i)

    % postimulus clustering
    if use_prev_cluster
        vertices_by_cluster_reduced(:, i) = vertices_by_cluster_reduced_prev(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            vertices_by_cluster_reduced(:, i) = kmeans(poststim_reduced, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            vertices_by_cluster_reduced(:, i) = kmedoids(poststim_reduced, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % silhouette values
    fprintf('\tCalculating silhouette values for %d clusters poststimulus\n', i)
    figure;
    s(:, i) = silhouette(poststim_reduced, squeeze(vertices_by_cluster_reduced(:, i)), 'correlation');
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['silhouettes_' num2str(i) '.tif']);
    for j = 1 : i
        s_mean(i, j) = mean(s(vertices_by_cluster_reduced(:, i) == j, i));
    end

    % cluster time series
    figure
    hold on
    legend_str = cell(1, i);
    h = NaN(i, total_num_bins);
    for j = 1 : i
        vals = vertex_values_sum_reduced(vertices_by_cluster_reduced(:, i) == j, :);
        plot(T - stimulus_offset, nanmean(vals, 1), 'Color', ColorSet(j, :), 'Linewidth', 2)
        meantrace(i, j, :) = nanmean(vals, 1);
        p_val(i, j) = ranksum(abs(squeeze(meantrace(i, j, 1 : stimulus_bin - 1))), abs(squeeze(meantrace(i, j, stimulus_bin : total_num_bins))));
        b_vals = vals(:, 1 : stimulus_bin - 1);
        b_mean = mean(b_vals(:));
        b_std = std(b_vals(:));
        for k = stimulus_bin : total_num_bins
            h(j, k) = ztest(vals(:, k), b_mean, b_std);
        end
        legend_str{j} = ['cluster ' num2str(j) ', p = ' num2str(p_val(i, j)) ', s = ' num2str(s_mean(i, j))];
    end
    h(h == 0) = NaN;
    for j = 1 : i
        plot(T - stimulus_offset, (-3 * h(j, :)) - .5 * (j - 1), '*', 'Color', ColorSet(j, :))
    end
    title(['# clusters = ' num2str(i)], 'Fontsize', 14);
    xlabel('time (s)', 'FontSize', 14)
    ylabel('z-score gamma power', 'FontSize', 14)
    xlim([T(1) - stimulus_offset, T(end) - stimulus_offset]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL, 'Color', 'k');
    line(xL, [0 0], 'Color', 'k');
    legend(legend_str,'Location', 'Northwest')
    legend('boxoff');
    set(gca, 'Fontsize', 14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf, ['cluster_timeseries_' num2str(i) '.tif']);

    % postimulus 0-250ms clustering
    if use_prev_cluster
        vertices_by_cluster_reduced_0_250ms(:, i) = vertices_by_cluster_reduced_prev_0_250ms(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            vertices_by_cluster_reduced_0_250ms(:,i) = kmeans(poststim_reduced_0_250ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            vertices_by_cluster_reduced_0_250ms(:,i) = kmedoids(poststim_reduced_0_250ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % silhouette values
    fprintf('\tCalculating silhouette values for %d clusters poststimulus 0-250ms\n', i)
    figure;
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['silhouettes_0_250_' num2str(i) '.tif']);
    s_0_250ms(:, i) = silhouette(poststim_reduced_0_250ms, squeeze(vertices_by_cluster_reduced_0_250ms(:, i)), 'correlation');
    for j = 1 : i
        s_mean_0_250ms(i, j) = mean(s_0_250ms(vertices_by_cluster_reduced_0_250ms(:, i) == j, i));
    end

    % cluster time series
    figure
    hold on
    legend_str = cell(1, i);
    h = NaN(i, total_num_bins);
    for j = 1 : i
        vals = vertex_values_sum_reduced(vertices_by_cluster_reduced_0_250ms(:, i) == j, :);
        plot(T - stimulus_offset, nanmean(vals, 1), 'Color', ColorSet(j, :), 'Linewidth', 2)
        meantrace_0_250ms(i, j, :) = nanmean(vals, 1);
        p_val(i, j) = ranksum(abs(squeeze(meantrace_0_250ms(i, j, 1 : stimulus_bin - 1))), abs(squeeze(meantrace_0_250ms(i, j, stimulus_bin : total_num_bins))));
        b_vals = vals(:, 1 : stimulus_bin - 1);
        b_mean = mean(b_vals(:));
        b_std = std(b_vals(:));
        for k = stimulus_bin : total_num_bins
            h(j, k) = ztest(vals(:, k), b_mean, b_std);
        end
        legend_str{j} = ['cluster ' num2str(j) ', p = ' num2str(p_val(i, j)) ', s = ' num2str(s_mean(i, j))];
    end
    h(h == 0) = NaN;
    for j = 1 : i
        plot(T - stimulus_offset, (-3 * h(j, :)) - .5 * (j - 1), '*', 'Color', ColorSet(j, :))
    end
    title(['# clusters = ' num2str(i)], 'Fontsize', 14);
    xlabel('time (s)', 'FontSize', 14)
    ylabel('z-score gamma power', 'FontSize', 14)
    xlim([T(1) - stimulus_offset, T(end) - stimulus_offset]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL, 'Color', 'k');
    line(xL, [0 0], 'Color', 'k');
    legend(legend_str,'Location', 'Northwest')
    legend('boxoff');
    set(gca, 'Fontsize', 14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf, ['cluster_timeseries_0_250_' num2str(i) '.tif']);

    % postimulus 250-500ms clustering
    if use_prev_cluster
        vertices_by_cluster_reduced_250_500ms(:, i) = vertices_by_cluster_reduced_prev_250_500ms(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            vertices_by_cluster_reduced_250_500ms(:,i) = kmeans(poststim_reduced_250_500ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            vertices_by_cluster_reduced_250_500ms(:,i) = kmedoids(poststim_reduced_250_500ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % silhouette values
    fprintf('\tCalculating silhouette values for %d clusters poststimulus 250-500ms\n', i)
    figure;
    s_250_500ms(:, i) = silhouette(poststim_reduced_250_500ms, squeeze(vertices_by_cluster_reduced_250_500ms(:, i)), 'correlation');
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['silhouettes_250_500_' num2str(i) '.tif']);
    for j = 1 : i
        s_mean_250_500ms(i, j) = mean(s_250_500ms(vertices_by_cluster_reduced_250_500ms(:, i) == j, i));
    end

    % cluster time series
    figure
    hold on
    legend_str = cell(1, i);
    h = NaN(i, total_num_bins);
    for j = 1 : i
        vals = vertex_values_sum_reduced(vertices_by_cluster_reduced_250_500ms(:, i) == j, :);
        plot(T - stimulus_offset, nanmean(vals, 1), 'Color', ColorSet(j, :), 'Linewidth', 2)
        meantrace_250_500ms(i, j, :) = nanmean(vals, 1);
        p_val(i, j) = ranksum(abs(squeeze(meantrace_250_500ms(i, j, 1 : stimulus_bin - 1))), abs(squeeze(meantrace_250_500ms(i, j, stimulus_bin : total_num_bins))));
        b_vals = vals(:, 1 : stimulus_bin - 1);
        b_mean = mean(b_vals(:));
        b_std = std(b_vals(:));
        for k = stimulus_bin : total_num_bins
            h(j, k) = ztest(vals(:, k), b_mean, b_std);
        end
        legend_str{j} = ['cluster ' num2str(j) ', p = ' num2str(p_val(i, j)) ', s = ' num2str(s_mean(i, j))];
    end
    h(h == 0) = NaN;
    for j = 1 : i
        plot(T - stimulus_offset, (-3 * h(j, :)) - .5 * (j - 1), '*', 'Color', ColorSet(j, :))
    end
    title(['# clusters = ' num2str(i)], 'Fontsize', 14);
    xlabel('time (s)', 'FontSize', 14)
    ylabel('z-score gamma power', 'FontSize', 14)
    xlim([T(1) - stimulus_offset, T(end) - stimulus_offset]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL, 'Color', 'k');
    line(xL, [0 0], 'Color', 'k');
    legend(legend_str,'Location', 'Northwest')
    legend('boxoff');
    set(gca, 'Fontsize', 14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf, ['cluster_timeseries_250_500_' num2str(i) '.tif']);

    % postimulus 500-750ms clustering
    if use_prev_cluster
        vertices_by_cluster_reduced_500_750ms(:, i) = vertices_by_cluster_reduced_prev_500_750ms(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            vertices_by_cluster_reduced_500_750ms(:,i) = kmeans(poststim_reduced_500_750ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            vertices_by_cluster_reduced_500_750ms(:,i) = kmedoids(poststim_reduced_500_750ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % silhouette values
    fprintf('\tCalculating silhouette values for %d clusters poststimulus 500-750ms\n', i)
    figure;
    s_500_750ms(:, i) = silhouette(poststim_reduced_500_750ms, squeeze(vertices_by_cluster_reduced_500_750ms(:, i)), 'correlation');
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['silhouettes_500_750_' num2str(i) '.tif']);
    for j = 1 : i
        s_mean_500_750ms(i, j) = mean(s_500_750ms(vertices_by_cluster_reduced_500_750ms(:, i) == j, i));
    end

    % cluster time series
    figure
    hold on
    legend_str = cell(1, i);
    h = NaN(i, total_num_bins);
    for j = 1 : i
        vals = vertex_values_sum_reduced(vertices_by_cluster_reduced_500_750ms(:, i) == j, :);
        plot(T - stimulus_offset, nanmean(vals, 1), 'Color', ColorSet(j, :), 'Linewidth', 2)
        meantrace_500_750ms(i, j, :) = nanmean(vals, 1);
        p_val(i, j) = ranksum(abs(squeeze(meantrace_500_750ms(i, j, 1 : stimulus_bin - 1))), abs(squeeze(meantrace_500_750ms(i, j, stimulus_bin : total_num_bins))));
        b_vals = vals(:, 1 : stimulus_bin - 1);
        b_mean = mean(b_vals(:));
        b_std = std(b_vals(:));
        for k = stimulus_bin : total_num_bins
            h(j, k) = ztest(vals(:, k), b_mean, b_std);
        end
        legend_str{j} = ['cluster ' num2str(j) ', p = ' num2str(p_val(i, j)) ', s = ' num2str(s_mean(i, j))];
    end
    h(h == 0) = NaN;
    for j = 1 : i
        plot(T - stimulus_offset, (-3 * h(j, :)) - .5 * (j - 1), '*', 'Color', ColorSet(j, :))
    end
    title(['# clusters = ' num2str(i)], 'Fontsize', 14);
    xlabel('time (s)', 'FontSize', 14)
    ylabel('z-score gamma power', 'FontSize', 14)
    xlim([T(1) - stimulus_offset, T(end) - stimulus_offset]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL, 'Color', 'k');
    line(xL, [0 0], 'Color', 'k');
    legend(legend_str,'Location', 'Northwest')
    legend('boxoff');
    set(gca, 'Fontsize', 14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf, ['cluster_timeseries_500_750_' num2str(i) '.tif']);

    % postimulus 750-1000ms clustering
    if use_prev_cluster
        vertices_by_cluster_reduced_750_1000ms(:, i) = vertices_by_cluster_reduced_prev_750_1000ms(:, i);
    else
        if strcmp(cluster_method, 'kmeans')
            vertices_by_cluster_reduced_750_1000ms(:,i) = kmeans(poststim_reduced_750_1000ms, i, 'dist', 'correlation', 'rep', 20);
        elseif strcmp(cluster_method, 'kmedoids')
            vertices_by_cluster_reduced_750_1000ms(:,i) = kmedoids(poststim_reduced_750_1000ms, i, 'Distance', 'correlation', 'replicates', 20);
        end
    end

    % silhouette values
    fprintf('\tCalculating silhouette values for %d clusters poststimulus 750-1000ms\n', i)
    figure;
    s_750_1000ms(:, i) = silhouette(poststim_reduced_750_1000ms, squeeze(vertices_by_cluster_reduced_750_1000ms(:, i)), 'correlation');
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['silhouettes_750_1000_' num2str(i) '.tif']);
    for j = 1 : i
        s_mean_750_1000ms(i, j) = mean(s_750_1000ms(vertices_by_cluster_reduced_750_1000ms(:, i) == j, i));
    end

    % cluster time series
    figure
    hold on
    legend_str = cell(1, i);
    h = NaN(i, total_num_bins);
    for j = 1 : i
        vals = vertex_values_sum_reduced(vertices_by_cluster_reduced_750_1000ms(:, i) == j, :);
        plot(T - stimulus_offset, nanmean(vals, 1), 'Color', ColorSet(j, :), 'Linewidth', 2)
        meantrace_750_1000ms(i, j, :) = nanmean(vals, 1);
        p_val(i, j) = ranksum(abs(squeeze(meantrace_750_1000ms(i, j, 1 : stimulus_bin - 1))), abs(squeeze(meantrace_750_1000ms(i, j, stimulus_bin : total_num_bins))));
        b_vals = vals(:, 1 : stimulus_bin - 1);
        b_mean = mean(b_vals(:));
        b_std = std(b_vals(:));
        for k = stimulus_bin : total_num_bins
            h(j, k) = ztest(vals(:, k), b_mean, b_std);
        end
        legend_str{j} = ['cluster ' num2str(j) ', p = ' num2str(p_val(i, j)) ', s = ' num2str(s_mean(i, j))];
    end
    h(h == 0) = NaN;
    for j = 1 : i
        plot(T - stimulus_offset, (-3 * h(j, :)) - .5 * (j - 1), '*', 'Color', ColorSet(j, :))
    end
    title(['# clusters = ' num2str(i)], 'Fontsize', 14);
    xlabel('time (s)', 'FontSize', 14)
    ylabel('z-score gamma power', 'FontSize', 14)
    xlim([T(1) - stimulus_offset, T(end) - stimulus_offset]);
    ylim(y_limits)
    xL = xlim;
    yL = ylim;
    line([0 0], yL, 'Color', 'k');
    line(xL, [0 0], 'Color', 'k');
    legend(legend_str,'Location', 'Northwest')
    legend('boxoff');
    set(gca, 'Fontsize', 14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf, ['cluster_timeseries_750_1000_' num2str(i) '.tif']);
end

% initialize cluster tracking variables
vertices_by_cluster = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_0_250ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_250_500ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_500_750ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_750_1000ms = zeros(numel(changing_vertices_idx), max_clusters);

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

    save('clustering.mat', 'vertices_by_cluster_L', 'vertices_by_cluster_R',...
        'vertices_by_cluster_L_0_250ms', 'vertices_by_cluster_R_0_250ms',...
        'vertices_by_cluster_L_250_500ms', 'vertices_by_cluster_R_250_500ms',...
        'vertices_by_cluster_L_500_750ms', 'vertices_by_cluster_R_500_750ms',...
        'vertices_by_cluster_L_750_1000ms', 'vertices_by_cluster_R_750_1000ms',...
        'changing_vertices_idx',...
        's', 's_mean', 'p_val',...
        's_0_250ms', 's_mean_0_250ms', 'p_val_0_250ms',...
        's_250_500ms', 's_250_500ms_mean', 'p_val_250_500ms',...
        's_500_750ms', 's_mean_500_750ms', 'p_val_500_750ms',...
        's_750_1000ms', 's_mean_750_1000ms', 'p_val_750_1000ms', '-v7.3')
else
    save('clustering.mat', 'vertices_by_cluster',...
        'vertices_by_cluster_0_250ms',...
        'vertices_by_cluster_250_500ms',...
        'vertices_by_cluster_500_750ms',...
        'vertices_by_cluster_750_1000ms',...
        'changing_vertices_idx',...
        's', 's_mean', 'p_val',...
        's_0_250ms', 's_mean_0_250ms', 'p_val_0_250ms',...
        's_250_500ms', 's_250_500ms_mean', 'p_val_250_500ms',...
        's_500_750ms', 's_mean_500_750ms', 'p_val_500_750ms',...
        's_750_1000ms', 's_mean_750_1000ms', 'p_val_750_1000ms', '-v7.3')
end
