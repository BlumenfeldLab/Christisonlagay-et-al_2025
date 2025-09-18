
% specify pipeline file to use
pipeline_file = 'normal_pipeline_file_info_21_04_14.mat';

% specify subjects to process
subjects_to_process = [27 : 31];
if ispc
    source = '\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\';
    load(['\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\', pipeline_file])
else
    source = '/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/';
    load(['/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/normal_pipeline_file_info/', pipeline_file])
end
for i = subjects_to_process

    fprintf('\n%d / %d', i, 31)
    source_file = [source, patients{i}, '/icEEG/HNCT Auditory ID Combined/sorted_trials', sorted_trials_extended_r_suffixes{i}, '.mat'];
    load(source_file)

    % extract -2 to 2 s
    BC_trials = BC_trials(:, 1 + floor(1 * fs(i)) : floor(5 * fs(i)), :);
    WG_trials = WG_trials(:, 1 + floor(1 * fs(i)) : floor(5 * fs(i)), :);
    %blank_TN_trials = blank_TN_trials(:, 1 + floor(1 * fs(i)) : floor(3 * fs(i)), :);
    save([source, patients{i}, '/icEEG/HNCT Auditory ID Combined/sorted_trials_sounds_restricted_CPonly_2to2_r.mat'], 'BC_trials', 'WG_trials', '-v7.3')
end
