figure('Color','k','position',[20 72 800 600])
toolboxpath = fileparts(which('mni2fs'));

%%%%go here! https://scalablebrainatlas.incf.org/human/HOA06

roi_labels=[
    {'background'}
    {'FP'}
    {'Ins'}
    {'SFG'}
    {'MFG'}
    {'IFGpt'}
    {'IFGpo'}
    {'PrG'}
    {'TP'}
    {'STGa'}
    {'STGp'}
    {'MTGa'}
    {'MTGp'}
    {'MTGtp'}
    {'ITGa'}
    {'ITGp'}
    {'ITGtp'}
    {'PoG'}
    {'SPL'}
    {'SmGa'}
    {'SmGp'}
    {'AG'}
    {'LOCs'}
    {'LOCi'}
    {'IcC'}
    {'FMC'}
    {'SMC'}
    {'ScC'}
    {'PcG'}
    {'CGa'}
    {'CGp'}
    {'PcC'}
    {'CC'}
    {'FOC'}
    {'PhGa'}
    {'PaGp'}
    {'LG'}
    {'TFCa'}
    {'TFCp'}
    {'TOF'}
    {'OFG'}
    {'FOpC'}
    {'COpC'}
    {'POpC'}
    {'PP'}
    {'H1/H2'}
    {'PT'}
    {'SccC'}
    {'OcP'}
    {'49'}
    {'50'}
    {'51'}
    {'52'}
    {'53'}
    {'54'}
    {'55'}
    {'56'}
    {'57'}
    {'58'}
    {'59'}
    ];
close all
% Replace the following path with the path to the mni2fs toolbox folder

toolboxpath = fileparts(which('mni2fs'));




figure('Color','k','position',[20 72 800 600])

% Load and Render the FreeSurfer surface
S = [];
mni2fs_dir = 'N:\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master';
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
S.hem = 'rh'; % choose the hemesphere 'lh' or 'rh'
S.inflationstep = 5; % 1 no inflation, 6 fully inflated
S.plotsurf = 'inflated';
S.lookupsurf = 'mid';
S.decimation = false; % Decimate the surface for speed. (Use FALSE for publishable quality figures).
S = mni2fs_brain(S);
view([-90 0])

% table_of_faces=table;


% S.roicolorspec =colorschoice(rois, :); % color. Can also be a three-element vector
S.roialpha = .5; % transparency 0-1
S.roismoothdata=1;
S.mnivol = fullfile('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\NIFTI\HarvardOxford-cort-maxprob-thr0-1mm.nii');
D.reps=6;
D.thresh=200;
[S table_of_facesHOA] = mni2fsHOA_roi_kcl(S, roi_labels, D);
% 
coordinatesfilename=[num2str(size(table_of_facesHOA,2)) '_HOAROIs_min' num2str(D.thresh) '_smoothing' num2str(D.reps) '_'  S.hem '.mat'];
save(coordinatesfilename, 'table_of_facesHOA');




