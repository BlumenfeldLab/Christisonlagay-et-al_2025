
% Load pipeline metadata 
load('\\172.23.254.106\data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info_tactile\normal_pipeline_file_info_20_12_15.mat')

% Remove patient 511WC as gray matter selection is not yet done

% Count gray matter electrodes
elec_count = 0;
for i = 1 : 6
    cd(['\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\', patients{i},...
        '\icEEG\HNCT Tactile ID Combined'])
    labels = [];
    try
        load('labels_all_gray_matter.mat');
    catch
        load('labels_first_surgery_all_gray_matter.mat');
    end
    elec_count = elec_count + numel(labels);
end