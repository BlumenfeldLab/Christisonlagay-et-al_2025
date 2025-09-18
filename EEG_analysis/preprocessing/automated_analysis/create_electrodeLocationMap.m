
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_20_08_31.mat');
data_location = 'E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies';
patients = patients(1 : 31);
createMontage = 0;
laterality = 0;
inflationstep = 5;
pre_proj_plot = 0;

electrodeLocationMap(data_location, patients, createMontage, laterality, inflationstep, pre_proj_plot)