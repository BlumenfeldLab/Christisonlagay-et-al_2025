
% define the pipeline info file
pipeline_file = 'normal_pipeline_file_info_21_04_14.mat';

if ispc
    load(['//172.23.254.106/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/normal_pipeline_file_info/', pipeline_file])
    data_folder = '//172.23.254.106/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/';
else
    load(['/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/normal_pipeline_file_info/', pipeline_file])
    data_folder = '/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/';
end

gamma_script = 'common';  % can be 'normal', 'common', 'or 'common_equal'

subjects_to_process = 27 : 31;
for i = subjects_to_process
    
    % 415DO cannot be processed using version 2019 and must be processed
    % using 2017
    if strcmp(patients{i}, '415DO') && contains(version, '2019')
        continue;
    elseif (strcmp(patients{i}, '415DO') && contains(version, '2017'))...
            || ~strcmp(patients{i}, '415DO')
        
        % Display current subject
        fprintf(sprintf('\nSubject %s (%d/%d)', patients{i}, i, length(patients)));
        
        % Move to subject directory
        cd([data_folder, patients{i},'/icEEG/HNCT Auditory ID Combined'])
        folderpath = pwd;
        
        % Load subject labels
        try
            load('labels_all_gray_matter.mat')
        catch
            load('labels_first_surgery_all_gray_matter.mat')
        end
        
        % FREQUENCY ANALYSIS PARAMETERS
        full = false;
        raw = false;
        baseline_opts.method = 'z-score';
        baseline_opts.start_time = -.5;
        stimulus_offset = 2;
        power_opts.method = 'filter';
        power_opts.filter_window = .1;
        power_opts.filter_slide = .025;
        power_opts.window_size = floor(fs(i) *  power_opts.filter_window);
        power_opts.sliding_window = floor(fs(i) * power_opts.filter_slide);
        
        % SORTED TRIALS SUFFIXES
        load(['sorted_trials', sorted_trials_2to2_r_suffixes{i}, '.mat'])
        
        % COMPOSITE REJECTIONS SUFFIXES
        load(['composite_rejection_recentered', composite_rejection_recentered_r_suffixes{i}, '.mat']);
        
        % FILE AND FOLDER SUFFIXES
        file_suffix = '_sounds_restricted_CPonly_100ms_filter_2to2_baseline5to0_common_r';
        folder_suffix = ' Sounds Restricted by CP Accuracy 100ms Filter 2to2 Baseline5to0 Common R';
        
        switch gamma_script
            case 'normal'
                % gamma extraction
                gamma_power_overlap_NCF(labels, BC_trials, WG_trials, noisy_BC_trials_idx, noisy_WG_trials_idx,...
                    fs(i), raw, full, baseline_opts, power_opts, stimulus_offset, file_suffix, folder_suffix, folderpath);
            case 'common'
                % gamma extraction common baseline
                gamma_power_overlap_NCF_common_baseline(labels, BC_trials, WG_trials, noisy_BC_trials_idx, noisy_WG_trials_idx,...
                    fs(i), raw, full, baseline_opts, power_opts, stimulus_offset, file_suffix, folder_suffix, folderpath);
            case 'common_equal'
                % gamma extraction common baseline equal weight
                gamma_power_overlap_NCF_common_baseline_equal_weight(labels, BC_trials, WG_trials, noisy_BC_trials_idx, noisy_WG_trials_idx,...
                    fs(i), raw, full, baseline_opts, power_opts, stimulus_offset, file_suffix, folder_suffix, folderpath);
        end
    end
end
