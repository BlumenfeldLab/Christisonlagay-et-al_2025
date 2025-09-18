% load roi
load('\\Mwmj04vg66\d\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\1_marsbarROIs_min100_smoothing3_rh.mat')
roi_to_convert = 'MNI_Amygdala_R_roi';
roi_faces = table_of_faces.(roi_to_convert) ~= 0;

% load surface
side = 2;
for side_index = side
    if side_index == 1
        hem = 'lh';
    else
        hem = 'rh';
    end
    mni2fs_dir = '\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master';
    surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);
    fs_surf(side_index) = export(gifti(surf_fn));
end

% convert faces to vertices
vertices = unique(reshape(fs_surf(side).faces(roi_faces, :), [sum(roi_faces) * 3, 1]));

% save roi vertices
save('\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\hnct_iceeg_scripting\brainmap_analysis\roi\MNI_Amygdala_R_roi.mat', 'vertices')
