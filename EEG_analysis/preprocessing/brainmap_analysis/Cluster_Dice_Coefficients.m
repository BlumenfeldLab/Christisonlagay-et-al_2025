% wendy_data = load(['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\',...
% '4th surface\Clustering 9 patients 57bins recentered rejoutliers\CP\',...
% 'vertex_values_bilat_9pts57bins_recentered_rejoutliers_CP.mat'], 'vertices_by_cluster');

wendy_data = load(['F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\',...
'Inflated Surface\9 Patients Visual 57 Bins Clustering\CP\',...
'vertex_values_bilat_9pts57bins_recentered_rejoutliers.mat'], 'vertices_by_cluster');

% chris_data = load(['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\4th surface\',...
%     'Clustering 9 patients 57bins calibration recentered rejoutliers\P\',...
%     'vertex_values_bilat_9pts57bins_calibration_recentered_rejoutliers_P.mat'], 'vertices_by_cluster');

chris_data = load(['F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\',...
'Inflated Surface\17 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy Clustering\CP\',...
'vertex_values_bilat_17pts57bins_recentered_rejoutliers_sounds_restricted_CPonly.mat'], 'vertices_by_cluster');

savedir = ['F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\',...
    '9 Patients Visual + 17 Patients Auditory ID Sounds Restricted by CP Accuracy 57 Bins Cluster Comparisons'];

coeff_mats = cell(10, 1);

figs = gobjects(1, 9);
ax = gobjects(1, 9);
dice_min_max = NaN(2, 9);

for num_clusters = 2:10

    % Dimensions: Wendy_cluster rows x Chris_cluster columns
    dice_coeffs = NaN(4, num_clusters);

    wendy_cluster_verts = wendy_data.vertices_by_cluster(:, 4);

    chris_cluster_verts = chris_data.vertices_by_cluster(:, num_clusters);

    for w = 1:4
        for c = 1:num_clusters

            wendy_cluster_inds = wendy_cluster_verts == w;
            chris_cluster_inds = chris_cluster_verts == c;

            dice_coeffs(w, c) = 2 * sum(wendy_cluster_inds & chris_cluster_inds) /...
                (sum(wendy_cluster_inds) + sum(chris_cluster_inds));

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

    xlabel('Auditory Clusters', 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')
    ylabel('Visual Clusters', 'FontSize', 17.6, 'Color', 'black', 'FontWeight', 'bold')

    yticks(1:4)
    dice_min_max(:, num_clusters - 1) = [min(dice_coeffs(:)); max(dice_coeffs(:))];
    
end

dice_min = min(dice_min_max(1, :));
dice_max = max(dice_min_max(2, :));

for a = 1:num_clusters - 1
    caxis(ax(a), [dice_min dice_max])
    saveas(figs(a), fullfile(savedir, ['Dice Coeffs 4 Visual ' num2str(a + 1) ' Auditory Clusters.tiff']))
end
