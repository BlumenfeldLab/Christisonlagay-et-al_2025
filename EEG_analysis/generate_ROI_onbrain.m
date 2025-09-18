% close all
% Generate ROIs on the inflated brain for MNI marsbar areas
% Written by Kate Christison-Lagay, Sept 2020

toolboxpath = fileparts(which('mni2fs'));
mappeddatadrive='V';


% roi_labels = [
%  {'MNI_Amygdala_L_roi'}
% {'MNI_Angular_L_roi'}
% {'MNI_Calcarine_L_roi'}
% % {'MNI_Caudate_L_roi'}
% {'MNI_Cingulum_Ant_L_roi'}
% {'MNI_Cingulum_Post_L_roi'}
% {'MNI_Cingulum_Mid_L_roi'}
% {'MNI_Cuneus_L_roi'}
% {'MNI_Frontal_Inf_Oper_L_roi'}
% {'MNI_Frontal_Inf_Orb_L_roi'}
% {'MNI_Frontal_Sup_Medial_L_roi'}
% {'MNI_Frontal_Mid_L_roi'}
% {'MNI_Frontal_Mid_Orb_L_roi'}
% {'MNI_Lingual_L_roi'}
% {'MNI_Frontal_Sup_Orb_L_roi'}
% {'MNI_Fusiform_L_roi'}
% {'MNI_Hippocampus_L_roi'}
% {'MNI_Insula_L_roi'}
% {'MNI_Occipital_Inf_L_roi'}
% {'MNI_Occipital_Mid_L_roi'}
% {'MNI_Parietal_Sup_L_roi'}
% {'MNI_Occipital_Sup_L_roi'}
% % {'MNI_Olfactory_L_roi'}
% % {'MNI_Pallidum_L_roi'}
% {'MNI_Paracentral_Lobule_L_roi'}
% {'MNI_ParaHippocampal_L_roi'}
% {'MNI_Parietal_Inf_L_roi'}
% {'MNI_Postcentral_L_roi'}
% {'MNI_Precentral_L_roi'}
% {'MNI_Precuneus_L_roi'}
% % {'MNI_Putamen_L_roi'}
% {'MNI_Supp_Motor_Area_L_roi'}
% {'MNI_SupraMarginal_L_roi'}
% {'MNI_Temporal_Pole_Mid_L_roi'}
% {'MNI_Temporal_Pole_Sup_L_roi'}
% {'MNI_Temporal_Sup_L_roi'}
% {'MNI_Thalamus_L_roi'}
% {'MNI_Temporal_Mid_L_roi'}
% % {'MNI_Olfactory_L_roi'}
% {'MNI_Rectus_L_roi'}
% {'MNI_Heschl_L_roi'}
% {'FEF_L'}
% ];


roi_labels = [
 {'MNI_Amygdala_R_roi'}
{'MNI_Angular_R_roi'}
{'MNI_Calcarine_R_roi'}
% {'MNI_Caudate_R_roi'}
{'MNI_Cingulum_Ant_R_roi'}
{'MNI_Cingulum_Post_R_roi'}
{'MNI_Cingulum_Mid_R_roi'}
{'MNI_Cuneus_R_roi'}
{'MNI_Frontal_Inf_Oper_R_roi'}
{'MNI_Frontal_Inf_Orb_R_roi'}
{'MNI_Frontal_Sup_Medial_R_roi'}
{'MNI_Frontal_Mid_R_roi'}
{'MNI_Frontal_Mid_Orb_R_roi'}
{'MNI_Lingual_R_roi'}
{'MNI_Frontal_Sup_Orb_R_roi'}
{'MNI_Fusiform_R_roi'}
{'MNI_Hippocampus_R_roi'}
{'MNI_Insula_R_roi'}
{'MNI_Occipital_Inf_R_roi'}
{'MNI_Occipital_Mid_R_roi'}
{'MNI_Parietal_Sup_R_roi'}
{'MNI_Occipital_Sup_R_roi'}
% {'MNI_Olfactory_R_roi'}
% {'MNI_Pallidum_R_roi'}
{'MNI_Paracentral_Lobule_R_roi'}
{'MNI_ParaHippocampal_R_roi'}
{'MNI_Parietal_Inf_R_roi'}
{'MNI_Postcentral_R_roi'}
{'MNI_Precentral_R_roi'}
{'MNI_Precuneus_R_roi'}
% {'MNI_Putamen_R_roi'}
{'MNI_Supp_Motor_Area_R_roi'}
{'MNI_SupraMarginal_R_roi'}
{'MNI_Temporal_Pole_Mid_R_roi'}
{'MNI_Temporal_Pole_Sup_R_roi'}
{'MNI_Temporal_Sup_R_roi'}
{'MNI_Thalamus_R_roi'}
{'MNI_Temporal_Mid_R_roi'}
% {'MNI_Olfactory_R_roi'}
{'MNI_Rectus_R_roi'}
{'MNI_Heschl_R_roi'}
];


figure('Color','k','position',[20 72 800 600])

% Load and Render the FreeSurfer surface
S = [];
mni2fs_dir = [mappeddatadrive ':\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master'];
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
S.hem = 'rh'; % choose the hemesphere 'lh' or 'rh'
S.inflationstep = 5; % 1 no inflation, 6 fully inflated
S.plotsurf = 'inflated';
S.lookupsurf = 'mid';
S.decimation = false; % Decimate the surface for speed. (Use FALSE for publishable quality figures).
S = mni2fs_brain(S);
view([-90 0])
colorschoice=distinguishable_colors(length(roi_labels));
thresholds=100;

table_of_faces=table;
for rois=1:length(roi_labels)
    current_area=char(roi_labels(rois))
    if ~strcmp(current_area, 'FEF_L')
        S.mnivol = fullfile(['\\172.23.254.106\Data5\Sharif Data F\HNCT fMRI Study\fMRI Analysis\Analysis Codes\Formatted ROIs\Marsbar ROIs\' current_area '\' current_area '.nii']);
    else
        S.mnivol = fullfile('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\NIFTI\FEFsphere_10--39_-7_56 left_roi.nii');
    end
    S.roicolorspec =colorschoice(rois, :); % color. Can also be a three-element vector
    S.roialpha = .5; % transparency 0-1
    S.roismoothdata=0;
    D.thresh=100;
    D.reps=4;
    [S output_array] = mni2fs_roi_kcl2(S, rois, D);
    singletab=table(output_array, 'VariableNames', roi_labels(rois));
    table_of_faces=[table_of_faces singletab];
    
end

% 
% coordinatesfilename=[num2str(length(roi_labels)) '_marsbarROIs_min' num2str(thresholds) '_smoothing' num2str(D.reps) '_'  S.hem '.mat'];
% save(coordinatesfilename, 'table_of_faces', 'roi_labels');
