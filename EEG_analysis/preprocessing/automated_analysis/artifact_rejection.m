
% specify the pipeline info file to use
pipeline_file = 'normal_pipeline_file_info_21_04_14.mat';

% path to all subject
if ispc
    data_folder = '//172.23.254.106/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/';
    load(['//172.23.254.106/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/normal_pipeline_file_info/', pipeline_file])
else
    data_folder = '/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/';
    load(['/mnt/Data25/HNCT/icEEG Analysis/Analysis/EEG_behavior/normal_pipeline_file_info/', pipeline_file])
end

% specify file and folder suffixes
file_suffix = '_sounds_restricted_CPonly_r';
folder_suffix = ' Sounds Restricted by CP Accuracy R';

% specify patient to use
patients_to_process = 23;  % index of patient to process in the 'patients' cell array of subject ID strings

for n = patients_to_process
    fprintf(sprintf('\nSubject (%d/%d)', n, length(patients)));

    % move to the patients data folder
    cd([data_folder, patients{n}, '/icEEG/HNCT Auditory ID Combined'])

    % load the electrode labels
    load('labels_all_gray_matter.mat')

    % load the sorted trials
    load(['sorted_trials', sorted_trials_r_suffixes{n}])

    % define some additional variables
    rootfolder = pwd;
    fs_subj = fs(n);

    % perform artifact rejection
    artifact_rejection_CJM(BC_trials, WG_trials, blank_TN_trials, labels, rootfolder,...
        file_suffix, folder_suffix, fs_subj);
end
