%%% Revised 7/31/2019 by Chris Micek to make the code more generalized,
%%% removing references to Wendy/Chris clusters from the visual/calibration
%%% data.


% wendy_data = load(['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\',...
% '4th surface\Clustering 9 patients 57bins recentered rejoutliers\CP\',...
% 'vertex_values_bilat_9pts57bins_recentered_rejoutliers_CP.mat'], 'vertices_by_cluster');

% wendy_data = load(['F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\',...
% 'Inflated Surface\9 Patients Visual 57 Bins Clustering\CP\',...
% 'vertex_values_bilat_9pts57bins_recentered_rejoutliers.mat'], 'vertices_by_cluster');

% clusters_1_data was previously wendy_data
clusters_1_data = load(['F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\',...
'Inflated Surface\Clustering 9 Patients Visual + 19 Patients Auditory ID Sounds Restricted by CP Accuracy 57 Bins'...
'\CP\vertex_vals_bilat_9vis_19aud_pts57bins_recentered_rejoutliers.mat'], 'vertices_by_cluster');

clusters_1_label = 'Visual + Auditory';

% chris_data = load(['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\4th surface\',...
%     'Clustering 9 patients 57bins calibration recentered rejoutliers\P\',...
%     'vertex_values_bilat_9pts57bins_calibration_recentered_rejoutliers_P.mat'], 'vertices_by_cluster');

% clusters_2_data was previously chris_data
clusters_2_data = load(['F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\',...
'Inflated Surface\Clustering 9 Patients Visual 57 Bins\CP\',...
'vertex_values_bilat_9pts57bins_recentered_rejoutliers.mat'], 'vertices_by_cluster');

clusters_2_label = 'Visual';

savedir = ['F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\',...
    '9 Patients Visual + 19 Patients Auditory ID Sounds Restricted by CP Accuracy 57 Bins Cluster Comparisons'];

max_clusters_1 = 10;
min_clusters_1 = 2;

max_clusters_2 = 4;
min_clusters_2 = 4;

coeff_mats = cell(min([max_clusters_1, max_clusters_2]), 1);

figs = gobjects(1, min([max_clusters_1, max_clusters_2]) - 1);
ax = gobjects(1, min([max_clusters_1, max_clusters_2]) - 1);
dice_min_max = NaN(2, min([max_clusters_1, max_clusters_2]) - 1);

for num_clusters = max([min_clusters_1, min_clusters_2]):min([max_clusters_1, max_clusters_2])

    % Dimensions: max_clusters_1 rows x max_cluster columns
%     dice_coeffs = NaN(4, num_clusters);
    dice_coeffs = NaN(num_clusters, num_clusters);

%     wendy_cluster_verts = wendy_data.vertices_by_cluster(:, 4);
    cluster_1_verts = clusters_1_data.vertices_by_cluster(:, num_clusters);

    cluster_2_verts = clusters_2_data.vertices_by_cluster(:, num_clusters);

%     for w = 1:4
    for c1 = 1:num_clusters
        for c2 = 1:num_clusters

            cluster_1_inds = cluster_1_verts == c1;
            cluster_2_inds = cluster_2_verts == c2;

            dice_coeffs(c1, c2) = 2 * sum(cluster_1_inds & cluster_2_inds) /...
                (sum(cluster_1_inds) + sum(cluster_2_inds));

        end
    end
    
    coeff_mats{num_clusters} = dice_coeffs;
    
    figs(num_clusters - 1) = figure;
    ax(num_clusters - 1) = gca;
    imagesc(dice_coeffs)
    colormap bone
    colorbar
    set(gca,'YDir','normal')
    xticks(1:num_clusters)
    box off
    title('Cluster Overlap S\orensen-Dice Index', 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')
%     xlabel('Calibration-Based Clusters', 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')
%     ylabel('Run-Based Clusters', 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')

%     xlabel('Auditory Clusters', 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')
    xlabel([clusters_1_label, ' Clusters'], 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')
    ylabel([clusters_2_label, ' Clusters'], 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')

%     yticks(1:4)
    yticks(1:num_clusters)
    dice_min_max(:, num_clusters - 1) = [min(dice_coeffs(:)); max(dice_coeffs(:))];
    
end

dice_min = min(dice_min_max(1, :));
dice_max = max(dice_min_max(2, :));

for a = (max([min_clusters_1, min_clusters_2]):min([max_clusters_1, max_clusters_2])) - 1
    caxis(ax(a), [dice_min dice_max])
    saveas(figs(a), fullfile(savedir, ['Dice Coeffs ', num2str(a + 1), ' ', clusters_1_label, ' ', num2str(a + 1),...
        ' ', clusters_2_label, ' Clusters.tiff']))

end