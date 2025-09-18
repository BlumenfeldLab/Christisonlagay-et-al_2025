clear all
mappeddatadrive='g:\mnt\Data25';

mni2fs_dir = [mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master'];
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
% load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
% Tmni = mnit1_1mm_stripped.transform;
inflationstep = 5;

brain.inflationstep=inflationstep;
brain.plotsurf='inflated';
brain.hem='lh';
brain.decimation=0;

brain=mni2fs_brain(brain);

brain.mnivol=mni2fs_load_nii('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\NIFTI\MNI_Cingulum_Ant_L_roi.nii');
brain.roicolorspec='r';

brain=mni2fs_roi(brain)

