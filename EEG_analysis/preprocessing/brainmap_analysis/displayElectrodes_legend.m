
load('E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info_19.10.03.mat')
figure
colors = distinguishable_colors(length(patients), .5 * ones(1, 3));
for i = 1 : length(patients)
    h(i) = plot(NaN, NaN, 'color', colors(i, :), 'LineWidth', 5); hold on
end
lgnd = legend(h, patients, 'Interpreter', 'none', 'color', 0.5 * ones(1, 3),...
    'FontSize', 14);
xticks([])
yticks([])