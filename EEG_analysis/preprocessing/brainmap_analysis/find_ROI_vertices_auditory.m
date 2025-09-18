
%% Make time series of individual ROIs from Fox et al. 2005 paper
rootfolder = 'N:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP Movies';
output_folder = 'E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\23 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy\Wendy Method\ROI';

% vertex_val_folder = 'F:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\17 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy';
vertex_val_folder = 'E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\23 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy\Wendy Method\Only Nearest Vertex';

trial_types = {'CP', 'CNP'};


cd(rootfolder)
xtick_spacing_filename = 'spectrogram_times_57bins.mat';
load(xtick_spacing_filename)
mni2fs_dir = 'E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master';
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
inflationstep = 5;

roi_dir = dir('E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\Formatted ROIs2');
all_folders = cell(length(roi_dir) - 2, 1);
for i = 3 : length(roi_dir)
   all_folders{i - 2} = ['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\Formatted ROIs2\', roi_dir(i).name]; 
end
roi_labels = {'calcarine', 'caudate', 'cuneus', 'dmn_lib', 'dmn', 'frontal_inf',...
              'frontal_med', 'frontal_mid', 'frontal_sup', 'heschl', 'insula',...
              'occipital', 'paracentral_lobule', 'parietal', 'postcentral',...
              'precentral', 'precuneus', 'temporal'};
colors = distinguishable_colors(length(roi_labels));

folders = all_folders([5:6, 7, 8, 33:34, 35, 36, 37:42, 43:44, 45:48, 49:54,...
    57:58, 61:62, 65:70, 77:78, 79:82, 83:84, 85:86, 99:108]);
roi_groups = {5:6, 7, 8, 33:34, 35, 36, 37:42, 43:44, 45:48, 49:54,...
    57:58, 61:62, 65:70, 77:78, 79:82, 83:84, 85:86, 99:108};
    
clusters = {};
output_suffix = '_grouped_ROIs';

for network_index = 1:length(folders)
    folder = folders{network_index};
    splitpath = strsplit(folder, filesep);
    ROI_name = splitpath{end};
    cd(folder);
    if ~exist([ROI_name, '.mat'], 'file')
        load('voxelInfo.mat');
        voxel_cor = xyz.vXYZ;
        cd(rootfolder);
        save([ROI_name, '.mat'],'voxel_cor','-v7.3');
    end
end

folders = [folders(:); clusters(:)]';
start_side = 1;
sides_num = 2;
cd(rootfolder)

for t = 1:length(trial_types)
    
    load([vertex_val_folder filesep trial_types{t} filesep 'vertex_values.mat'])
    vertex_valuesL(vertex_valuesL == 0) = NaN;
    vertex_valuesR(vertex_valuesR == 0) = NaN;
    vertex_values = {vertex_valuesL, vertex_valuesR};
    vertex_side_sums = cell(1, 2);
    for side_index = start_side:sides_num
        vertex_side_sums{side_index} = squeeze(nansum(vertex_values{side_index}, 2));
        divisor = squeeze(sum(~isnan(vertex_values{side_index}), 2));
        divisor(divisor == 0) = 1;
        vertex_side_sums{side_index} = vertex_side_sums{side_index} ./ sqrt(divisor);
    end
    vertex_values_bilateral = cat(1, vertex_valuesL, vertex_valuesR);
    vertex_side_sums_bilateral = squeeze(nansum(vertex_values_bilateral, 2));
    divisor = squeeze(sum(~isnan(vertex_values_bilateral), 2));
    divisor(divisor == 0) = 1;
    vertex_side_sums_bilateral = vertex_side_sums_bilateral ./ sqrt(divisor);

    for side_index = start_side:sides_num
        if side_index == 1
            side = 'L';
            hem = 'lh';
        else
            side = 'R';
            hem = 'rh';
        end
        
        surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);
        fs_surf(side_index) = export(gifti(surf_fn));
        v_num = size(fs_surf(side_index).vertices, 1);
        fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] *...
            Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
        fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
        surfrender_fn = fullfile(mni2fs_dir,['/surf/' hem '.inflated' num2str(inflationstep) '.surf.gii']);
        inflated_surf(side_index) = export(gifti(surfrender_fn));
        inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] *...
            Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
        inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);
    end

    ROI_unique = cell(1, length(folders));
    for network_index = 1:length(folders)
        splitpath = split(folders{network_index}, filesep);
        ROI_name = splitpath{end};
        side_delim = regexp(ROI_name, '_[LR]_roi');
        if ~isempty(side_delim)
            ROI_name_stripped = ROI_name(1:side_delim - 1);
        else
            ROI_name_stripped = ROI_name;
        end
        if ~any(strcmp(ROI_unique, ROI_name_stripped))
            ROI_unique{network_index} = ROI_name_stripped;
        end
    end

    ROI_unique = ROI_unique(~cellfun('isempty', ROI_unique));

    for i = length(ROI_unique):-1:1
        ROIs(i).name = '';
        ROIs(i).side_data = struct([]);
    end
    ROI_unique = cell(1, length(ROI_unique));
    cur_ROI = 1;
    advance = false;
    
    %% Find vertices within 15 mm of each electrode position
    for network_index = 1:length(folders)

        splitpath = split(folders{network_index}, filesep);
        ROI_name = splitpath{end};
        side_delim = regexp(ROI_name, '_[LR]_roi');
        
        if ~isempty(side_delim)
            ROI_name_stripped = ROI_name(1:side_delim - 1);
        else
            ROI_name_stripped = ROI_name;
        end

        load([ROI_name, '.mat']);

        if endsWith(ROI_name, 'L_roi')
            side_index = 1;
            side = 'L';
        elseif endsWith(ROI_name, 'R_roi')
            side_index = 2;
            side = 'R';
        else
            side_index = 1;
            side = 'Bilateral';
        end

        if strcmp(side, 'Bilateral')
            vertex_values = vertcat(vertex_side_sums{:});
            fs_vertices = vertcat(fs_surf.vertices);
            inflated_vertices = vertcat(inflated_surf.vertices);
        else
            vertex_values = vertex_side_sums{side_index};
            fs_vertices = fs_surf(side_index).vertices;
            inflated_vertices = inflated_surf(side_index).vertices;
        end

        % calculate the number of electrode of this network on this side
        % use the 3D coordinates of the projected electrode vertices
        if ~ismember('ROI_vertex_inds', who('-file', [ROI_name, '.mat']))
            ROI_coordinates = voxel_cor';
            ROI_vertex_inds = NaN(1, size(vertex_values, 1));
            Ind = 0;
            r = 1.5;
            if strcmp(ROI_name_stripped, 'MNI_Fusiform')
                for e = 1:size(ROI_coordinates, 1)
                    ROI_voxel_vertex_xy_dist = sqrt((ROI_coordinates(e,1) - fs_vertices(:,1)).^2 + (ROI_coordinates(e,2) - fs_vertices(:,2)).^2);
                    [~, minInd] = min(ROI_voxel_vertex_xy_dist);

                    if ROI_coordinates(e, 3) - fs_vertices(minInd, 3) < 181/3 && ROI_coordinates(e, 3) - fs_vertices(minInd, 3) > 181/60
                        ROI_coordinates(e, 3) = fs_vertices(minInd, 3);
                    end
                end
            end

            for e = 1:size(fs_vertices, 1)%ROIs_num
                % for each of the "electrodes" or ROIs, find the distance
                % between that and each vertex
                if any(all(fs_vertices(e, :) >= ROI_coordinates - r & fs_vertices(e, :) < ROI_coordinates + r, 2))
                    ROI_vertex_inds(e) = e;
                end
            end
            ROI_vertex_inds(isnan(ROI_vertex_inds) | ROI_vertex_inds == 0) = [];
            save([ROI_name, '.mat'], 'ROI_vertex_inds', '-append')
        end
        if ~any(strcmp(ROI_unique, ROI_name_stripped))
            ROIs(cur_ROI).name = ROI_name_stripped;
            ROI_unique{cur_ROI} = ROI_name_stripped;
            for i = 2:-1:1
                ROIs(cur_ROI).side_data(i).side = '';
                ROIs(cur_ROI).side_data(i).vertices = [];
                ROIs(cur_ROI).side_data(i).mean_power_timecourse = [];
            end
        else
            advance = true;
        end

        ROIs(cur_ROI).side_data(side_index).side = side;
        ROIs(cur_ROI).side_data(side_index).vertices = inflated_vertices(ROI_vertex_inds, :);
        ROIs(cur_ROI).side_data(side_index).mean_power_timecourse = nanmean(vertex_values(ROI_vertex_inds,:),1);
        if advance || strcmp(side, 'Bilateral')
            cur_ROI = cur_ROI + 1;
            advance = false;
        end

    end
    
    ROI_name_str = '';
    
    for cur_ROI = 1:length(ROIs)
        ROI_name_str = strcat(ROI_name_str,['_' ROIs(cur_ROI).name]);
    end
    
    %% plot ROIs in functional groups
    for cur_ROI = 1 : length(roi_groups)
        network_name = roi_labels{cur_ROI};
        for sub_roi = 1 : length(roi_groups{cur_ROI})
            if size(ROIs(cur_ROI).side_data, 2) > 1
                %network_vertices = cat(1, ROIs(cur_ROI).side_data(1).vertices,...
                %    ROIs(cur_ROI).side_data(2).vertices);
                timecourse(:, sub_roi) = nanmean(cat(1, ROIs(cur_ROI).side_data(1).mean_power_timecourse, ROIs(cur_ROI).side_data(2).mean_power_timecourse), 1);
            else
                %network_vertices = ROIs(cur_ROI).side_data(1).vertices;
                timecourse(:, sub_roi) = ROIs(cur_ROI).side_data(1).mean_power_timecourse;
            end
        end
        roi_timecourse(:, cur_ROI) = mean(timecourse, 2);
    end
    figure('Position', get(0, 'Screensize'))
    for i = 1 : length(roi_groups)
        plot((T - 1) * 1000, roi_timecourse(:, i), 'Linewidth', 2, 'Color', colors(i, :)); hold on
    end
    hold off
    title('All ROIs', 'Fontsize', 14, 'Interpreter', 'none');
    xlabel('Time (ms)', 'FontSize', 14)
    ylabel('Average z-score gamma power', 'FontSize', 14)
    xlim([(T(1) - 1) * 1000 (T(end) - 1) * 1000]);
    ylim([-4 10])
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(roi_labels, 'Location', 'Northwest', 'Interpreter', 'none')
    saveas(gcf, [output_folder filesep trial_types{t} output_suffix '.tif'])

    %% plot time-series using vertex time-series of each ROI
    %{
    for cur_ROI = 1:length(ROIs)
        network_name = ROIs(cur_ROI).name;
        if size(ROIs(cur_ROI).side_data, 2) > 1
            network_vertices = cat(1, ROIs(cur_ROI).side_data(1).vertices, ROIs(cur_ROI).side_data(2).vertices);
            timecourse = nanmean(cat(1, ROIs(cur_ROI).side_data(1).mean_power_timecourse, ROIs(cur_ROI).side_data(2).mean_power_timecourse), 1);
        else
            network_vertices = ROIs(cur_ROI).side_data(1).vertices;
            timecourse = ROIs(cur_ROI).side_data(1).mean_power_timecourse;
        end
        figure('Position', get(0, 'Screensize'))
        plot((T - 1) * 1000, timecourse, 'Linewidth', 2, 'Color', ColorSet(cur_ROI, :));
        title(network_name, 'Fontsize', 14, 'Interpreter', 'none');
        xlabel('Time (ms)', 'FontSize', 14)
        ylabel('Gamma power z-score', 'FontSize', 14)
        xlim([(T(1) - 1) * 1000 (T(end) - 1) * 1000]);
        ylim([-4 10])
        xL = xlim;
        yL = ylim;
        line([0 0], yL,'Color','k');
        line(xL, [0 0],'Color','k');
    end
    %}
    
    %% plot combined time-series of each network
    %{
    figure('Position', get(0, 'Screensize'))
    names = cell(1, length(ROIs));
    for cur_ROI = 1:length(ROIs)
        network_name = ROIs(cur_ROI).name;
        names{cur_ROI} = network_name;

        if size(ROIs(cur_ROI).side_data, 2) > 1
            network_vertices = cat(1, ROIs(cur_ROI).side_data(1).vertices, ROIs(cur_ROI).side_data(2).vertices);
            timecourse = nanmean(cat(1, ROIs(cur_ROI).side_data(1).mean_power_timecourse, ROIs(cur_ROI).side_data(2).mean_power_timecourse), 1);
        else
            network_vertices = ROIs(cur_ROI).side_data(1).vertices;
            timecourse = ROIs(cur_ROI).side_data(1).mean_power_timecourse;
        end

        plot((T-1) * 1000, timecourse, 'Linewidth', 2, 'Color', ColorSet(cur_ROI, :))
        hold on;
    end

    xlim([(T(1) - 1) * 1000 (T(end) - 1) * 1000])
    ylim([-4 7])
    xL = xlim;
    title(['ROI Timecourses ' trial_types{t} ' Trials'])
    xlabel('Time (ms)','FontSize',14)
    ylabel('Gamma power z-score','FontSize',14)

    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(names,'Location','Northwest', 'Interpreter', 'none')

    set(gca,'Fontsize',14)
    set(gcf, 'color', [1 1 1]);
    
    if ~exist(output_folder, 'dir')
        mkdir(output_folder)
    end
    saveas(gcf, [output_folder filesep trial_types{t} output_suffix '.tif'])
    %}
end
