
% save location
save_location = 'E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 77 Bins Sounds Restricted by CP Accuracy 100ms\filter_from2to2_baseline1to0r\CP_w250ms\';

% vertex values
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 77 Bins Sounds Restricted by CP Accuracy 100ms\filter_from2to2_baseline1to0r\Gamma_power_CP\vertex_values.mat')

% load cluster data
load('\\172.28.28.95\Data9\HNCT\auditory\results\31patients_77bins_100ms_filter_from2to2_baseline1to0\CP_clusters_kmeans_w250ms_nosig_differentsave\clustering_tmp_workspace.mat',...
    'vertices_by_cluster', 'vertices_by_cluster_0_250ms', 'vertices_by_cluster_250_500ms', 'vertices_by_cluster_500_750ms',...
    'vertices_by_cluster_750_1000ms')

cluster_method = 'kmeans'; % options: 'kmeans', 'kmedoids'
min_clusters = 2;
max_clusters = 6;
laterality = 0;
xtick_spacing_filename = 'filter_times_win100ms_shift25ms.mat';
stimulus_offset = 1;
use_prev_cluster = 0;

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
    
    % silhouette values
    fprintf('\tCalculating silhouette values for %d clusters poststimulus\n', i)
    figure;
    s(:, i) = silhouette(poststim(changing_vertices_idx, :), squeeze(vertices_by_cluster(changing_vertices_idx, i)), 'correlation');
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['silhouettes_' num2str(i) '.tif']);
    for j = 1 : i
        s_mean(i, j) = mean(s(vertices_by_cluster_reduced(:, i) == j, i));
    end
    
    % timeseries
    figure
    hold on
    legend_str = cell(1, i);
    h = NaN(i, total_num_bins);
    for j = 1 : i
        vals = vertex_values_sum(vertices_by_cluster(:, i) == j, :);
        plot(T - stimulus_offset, nanmean(vals, 1), 'Color', ColorSet(j, :), 'Linewidth', 2)
        meantrace = nanmean(vals, 1);
        p_val(i, j) = ranksum(abs(squeeze(meantrace(1 : stimulus_bin - 1))), abs(squeeze(meantrace(stimulus_bin : total_num_bins))));
        b_vals = vals(:, 1 : stimulus_bin - 1);
        b_mean = mean(b_vals(:));
        b_std = std(b_vals(:));
        for k = stimulus_bin : total_num_bins
            h(j, k) = ztest(vals(:, k), b_mean, b_std);
        end
        legend_str{j} = ['cluster ' num2str(j) ', p = ' num2str(p_val(i, j))];
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
end







