
% Load clusters to utilize
% This script assumes that that the reduced vertices are equal between CP
% and CNP groups
load(['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\',...
    'Group\Inflated Surface\23 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy\',...
    'Wendy Method\K-Means Clusters\CP_clusters\',...
    'vertex_values_bilat_23pts57bins_recentered_rejoutliers_sounds_restricted_CP.mat'],...
    'vertices_by_cluster_reduced', 'vertices_by_cluster_reduced_200ms', 'vertices_by_cluster_reduced_400ms')

% Vertex values to cluster
load(['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\',...
    '23 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy\Wendy Method\',...
    'CNP\vertex_values.mat'])

% Load relevant patient processing metadata
load('E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info_19.10.03.mat')
patients = patients(1 : end - 1); % remove patient 51

% Data location
data_location = 'E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies';

% Image storage folder. create if does not exist
image_storage_folder = ['E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\',...
    'HNCT ERP movies\Group\Inflated Surface\',...
    '23 Patients Auditory ID 57 Bins Sounds Restricted by CP Accuracy\tmpStorage\',...
    'Wendy Method\CNP_w_CP_clusters'];
if ~isfolder(image_storage_folder)
    mkdir(image_storage_folder)
end

% Define preferences
overlay = 0;
createMontage = 0;
laterality = 0;
views = 4;
inflationstep = 5;
subtraction = false;

% New addition: If IsVisual doesn't exist specifying which data is from the
% visual dataset, assume all data is not visual.
if ~exist('IsVisual', 'var')
    IsVisual = false(1, numel(patients));
end

num = squeeze(sum(vertex_valuesL ~= 0, 2));
num(num == 0) = 1;

% Determine if we want to divide by sqrt(num patients) (Wendys method) per
% vertex or to divide by num patients
wendys_method = 0;
if wendys_method
    if laterality == 0
        vertex_values_sum_L = squeeze(nansum(vertex_valuesL, 2)) ./ sqrt(num);
        num = squeeze(sum(vertex_valuesR ~=0, 2));
        num(num == 0) = 1;
        vertex_values_sum_R = squeeze(nansum(vertex_valuesR, 2)) ./ sqrt(num);
    elseif laterality == 3
        vertex_values_sum = squeeze(nansum(vertex_valuesL, 2)) ./ sqrt(num);
    end
else
    if laterality == 0
        vertex_values_sum_L = squeeze(nansum(vertex_valuesL, 2)) ./ num;
        num = squeeze(sum(vertex_valuesR ~=0, 2));
        num(num == 0) = 1;
        vertex_values_sum_R = squeeze(nansum(vertex_valuesR, 2)) ./ num;
    elseif laterality == 3
        vertex_values_sum = squeeze(nansum(vertex_valuesL, 2)) ./ num;
    end
end

xtick_spacing_filename = 'spectrogram_times_57bins.mat';
mni2fs_dir = 'E:\Chris HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master';
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
alpha = 0.8;

% the colors for the different patients' electrodes
colors = distinguishable_colors(length(patients), .5 * ones(1, 3));

%% Side laterality - which side of the brain to draw
start_side = 1;
if laterality == 0 % bilateral
    sides_num = 2;
elseif laterality == 1 % left only
    sides_num = 1;
elseif laterality == 2 % right only
    start_side = 2;
    sides_num = 2;
elseif laterality == 3 % do both left and right on the left brain
    sides_num = 1;
end

figHandle = gobjects(1, sides_num);

%% previously displayElectrodes_phase2_fig6_inflated 
for side_index = 1 : 2
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
        hem = 'lh';
    else
        side = 'R';
        hem = 'rh';
    end
    
    % hard coded number of frames for each overlay setting, need to change later
    if overlay == 1 || overlay == 0 % arrival times
        total_num_frames = 1;
    elseif overlay == 2 || overlay == 3
        total_num_frames = 57;
    end
    
    % Setup some mesh variables
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
    
    % create a figure of the left or right side of the brain
    figHandle(side_index) = figure('Position',[70,70,700,700]);
    temp_surf = [];
    temp_surf.hem = hem;
    temp_surf.inflationstep = inflationstep;
    temp_surf.decimation = 0;
    temp_surf = mni2fs_brain(temp_surf);
    set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
    set(gcf, 'color', [1 1 1]);
    axis off;
    axis equal;
    vertex_values = zeros(v_num, length(patients));
    
    % if you are plotting arrival times, then each frame is a different set
    % of electrodes
    if overlay == 1
        electrode_mapping_times = 10;
    else
        electrode_mapping_times = 1;
    end
    
    for e = 1 : electrode_mapping_times
        scatter_count = 0;
        for p = 1 : length(patients)
            % Get electrode info for the patient
            % this part plots the electrodes for each patient and the overlay data
            % for each patient
            patient_folder = [data_location, '/', patients{p}];
            
            display(patients{p})
            
            % Now load some data files needed for electrode placement and surface shading
            cd([data_location, '/', patients{p}]);
            load([side '_MontageMap.mat']);
            electrode = eval([side '_MontageMap' ]);
            
            %         electrode_num = size(electrode,1);
            % load the X-flipped montage and add to list of electrodes
            if laterality == 3
                load('R_MontageMap_flipX.mat');
                electrode = cat(1,electrode,R_MontageMap);
            end
            
            % Discard depth electrodes or plot only some electrodes for arrival times
            if IsVisual(p) && ~isempty(electrode)
                % find out which depth electrodes should be discarded from the
                % analysis
                if ~strcmp(patients{p}, '193AF')
                    load('labels.mat')
                    load('labels_depth_most.mat')
                else
                    load('labels_first_surgery.mat')
                    labels_depth = {};
                end
                
                % find the indices of the white matter electrodes in labels
                labels_depth_idx = [];
                for l = 1:length(labels_depth)
                    labels_depth_idx = [labels_depth_idx find(strcmp(labels_depth{l},labels))];
                end
                
                % remove them from the electrode list that you will later access
                [side_labels_depth, depth_idx, electrode_idx] = intersect(labels_depth_idx,electrode(:,1));
                electrode(electrode_idx,:) = [];
                
                % if you are trying to plot arrival times, then you will need
                % to plot only the electrodes in that time period
                if strcmp(overlay, 'Arrival Times')
                    load('arrival_times_by_period_100ms.mat')
                    electrodes_this_period = find(arrival_times_by_period(:,e) == 1);
                    electrode = electrode(ismember(electrode(:,1),electrodes_this_period),:);
                end
            end
            
            % Plot patient's electrodes
            if isempty(electrode) % Some files only have electrodes on the one side
                disp('-- No electrode data!');
                electrodes_present(side_index,p) = 0;
            else
                all_electrodes = electrode(:,1);
                electrodes_present(side_index,p) = 1;
                [elecVertices] = ProjectElectrode2TransSurf_inflated(v_num, fs_surf(side_index).vertices, all_electrodes, electrode);
                scatter_count = scatter_count + 1;
            end
            
            % cycle through all the frames to aggregate zscores for each vertex
            if overlay == 2 || overlay == 3
                for frame_index = 1:total_num_frames
                    if isempty(electrode) % Some files only have electrodes on the left, hence the if statement
                        disp('-- No electrode data!');
                        FaceVertexAlphaData    = zeros(v_num, 1);
                        vertexCdata            = zeros(v_num, 1);
                    else
                        disp([patients{p} ' Aggregating Frame ' num2str(frame_index)])
                        [electrode_vertex_values] = electrode_data_overlay(overlay,frame_index,all_electrodes, subtraction, file_suffixes{p});
                        [vertexCdata]        = CalculateTransparency(v_num, inflated_surf(side_index).vertices, elecVertices, electrode_vertex_values, overlay);
                        vertex_values(:,p,frame_index) = vertexCdata;
                    end
                end
            end
        end
        
        % Frame Generation if saving without overlaying color onto the brain
        
        % if you are doing arrival times, you need to save the figures at
        % the end of every time period and redraw the electrodes
        if overlay == 1
            disp(e)
            % Grab the figure of the side of the brain you are working on
            figure(figHandle(side_index))
            
            cd(image_storage_folder)
            frames_set(side_index,e,:) = savebrainimages(side,e);
            
            % remove the electrodes to make the next one
            children = get(gca,'children');
            delete(children(1:scatter_count,1))
        end
    end
end

% final colors for k-means using confirmed not perceived trials:
% dark blue DMN
% red visual assocation
% green visual cortex
% darker cyan fronto-parietal
ColorSetMaps = [
    1.00 1.00 1.00 % no cluster - white
    0.00 1.00 0.00 % Data 1 - green
    0.00 0.00 1.00 % Data 2 - blue
    1.00 0.00 0.00 % Data 3 - red
    0.00 1.00 1.00 % Data 4 - cyan
    1.00 0.00 1.00 % Data 5 - magenta
    0.75 0.75 0.00 % Data 6 - RGB
    0.25 0.25 0.25 % Data 7
    0.75 0.25 0.25 % Data 8
    0.95 0.95 0.00 % Data 9
    0.25 0.25 0.75 % Data 10
    ];
ColorSet = [
    0.00 1.00 0.00 % Data 1 - green
    0.00 0.00 1.00 % Data 2 - blue
    1.00 0.00 0.00 % Data 3 - red
    0.00 1.00 1.00 % Data 4 - cyan
    1.00 0.00 1.00 % Data 5 - magenta
    0.75 0.75 0.00 % Data 6 - RGB
    0.25 0.25 0.25 % Data 7
    0.75 0.25 0.25 % Data 8
    0.95 0.95 0.00 % Data 9
    0.25 0.25 0.75 % Data 10
    ];

%% previously displayElectrodes_phase3_fig6_inflated

cd(image_storage_folder)

% define the range of number of clusters you want
min_clusters = 2;
max_clusters = 10;

% combine left and right side vertex time series
if laterality == 0
    vertex_values_sum_bilateral = cat(1,vertex_values_sum_L,vertex_values_sum_R);
    poststim = vertex_values_sum_bilateral(:,29:end);
elseif laterality == 3
    poststim = vertex_values_sum(:,29:57);
end

% Only want to do kmeans clustering for post stimulus ERPs
% get rid of the vertices that have very small standard deviations but
% remember them so you can accurately assign colors to the correct
% vertices later
poststim_std = std(poststim, 0, 2);
constant_vertices_idx = poststim_std == 0;
changing_vertices_idx = poststim_std > 0;
vertices_by_cluster = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_200ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_400ms = zeros(numel(changing_vertices_idx), max_clusters);

% take only the vertices that are not constant
poststim_reduced = poststim(changing_vertices_idx,:);
if laterality == 0
    vertex_values_sum_bilateral_reduced = vertex_values_sum_bilateral(changing_vertices_idx,:);
    vertex_values_sum_reduced = vertex_values_sum_bilateral_reduced;
else
    vertex_values_sum_reduced = vertex_values_sum(changing_vertices_idx,:);
end

% redo the above for the first 200 ms (4 frames) post-stimulus
poststim_reduced_200ms = poststim_reduced(:,1:7);
poststim_reduced_400ms = poststim_reduced(:,1:13);

% Initialize vars for sillouette values, average sillouette values and p
% values per cluster
s = NaN(sum(changing_vertices_idx), max_clusters);
s_mean = NaN(max_clusters);
p_val = NaN(max_clusters);
s_200ms = NaN(sum(changing_vertices_idx), max_clusters);
s_200ms_mean = NaN(max_clusters);
p_val_200ms = NaN(max_clusters);
s_400ms = NaN(sum(changing_vertices_idx),max_clusters);
s_400ms_mean = NaN(max_clusters);
p_val_400ms = NaN(max_clusters);

% Initialize vars for something...
meantrace = NaN(max_clusters, max_clusters, 57);
baseline_trace = NaN(max_clusters, 28);
poststim_trace = NaN(max_clusters, 29);
meantrace_200ms = NaN(max_clusters, max_clusters, 57);
baseline_trace_200ms = NaN(max_clusters, 28);
poststim_trace_200ms = NaN(max_clusters, 29);
meantrace_400ms = NaN(max_clusters, max_clusters, 57);
baseline_trace_400ms = NaN(max_clusters, 28);
poststim_trace_400ms = NaN(max_clusters, 29);

load('spectrogram_times_57bins.mat')
for i = min_clusters : max_clusters
    
    % Indicate status
    fprintf('\Computing %d clusters\n', i)
    
    % Graph silhouette values
    disp(['Calculating silhouette values for ' num2str(i) ' clusters'])
    figure;
    cidx = vertices_by_cluster_reduced(:, i);
    [s(:, i), ~] = silhouette(poststim_reduced, cidx, 'correlation');
    
    % calculate the mean of the silhouette values in that cluster
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced(:, i) == j);
        s_mean(i, j) = mean(s(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_' num2str(i) '.tif']);
    
    % Graph cluster time courses
    figure
    legend_str = cell(1, i);
    % Gamma part earlier in the code because of graphing reasons
    for j = 1:i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced(:,i) == j);
        vals = vertex_values_sum_reduced(cluster_vertex_idx, :);
        plot(T-1,nanmean(vals,1), 'Color', ColorSet(j,:),'Linewidth',2); hold on;
        meantrace(i,j,:) = nanmean(vals,1);
        baseline_trace(j,:) = abs(meantrace(i,j,1:28));
        poststim_trace(j,:) = abs(meantrace(i,j,29:57));
        [p_val(i,j)] = ranksum(baseline_trace(j,:),poststim_trace(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val(i,j)) ', s = ' num2str(s_mean(i,j))];
    end
    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power z-score','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim([-.5 1])
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest')
    legend('boxoff');
    set(gca,'Fontsize',14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf,['Gamma_clusters_' num2str(i) '.tif']);
    
    clearvars legend_str
    
    %% K means clustering of the first 200 ms (4 frames) post-stimulus

    % Graph silhouette values
    disp(['Calculating silhouette values for ' num2str(i) ' clusters'])
    figure
    cidx = vertices_by_cluster_reduced_200ms(:, i);
    [s_200ms(:, i), ~] = silhouette(poststim_reduced, cidx,'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_200ms(:, i) == j);
        s_200ms_mean(i, j) = mean(s_200ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_200ms_' num2str(i) '.tif']);
    
    % Gamma part earlier in the code because of graphing reasons
    for j = 1:i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_200ms(:,i) == j);
        % plot that index
        plot(T-1,nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1), 'Color', ColorSet(j,:)); hold on;
        
        % do stats on the prestim and poststim gamma power zscore changes
        meantrace_200ms(i,j,:) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1);
        baseline_trace_200ms(j,:) = abs(meantrace_200ms(i,j,1:28));
        poststim_trace_200ms(j,:) = abs(meantrace_200ms(i,j,29:57));
        [p_val_200ms(i,j)] = ranksum(baseline_trace_200ms(j,:),poststim_trace_200ms(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_200ms(i,j)) ', s = ' num2str(s_200ms_mean(i,j))];
    end
    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim([-.5 1])
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest','AutoUpdate', 'off')
    legend('boxoff');
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T-1,squeeze(meantrace(i,j,:)), 'Color',[.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_200ms_' num2str(i) '.tif']);
    
    clearvars legend_str
    
    %% K means clustering of the first 400 ms (4 frames) post-stimulus   
    
    % Graph silhouette values
    figure
    disp(['Calculating silhouette values for ' num2str(i) ' clusters'])
    cidx = vertices_by_cluster_reduced_400ms(:, i);
    [s_400ms(:, i), ~] = silhouette(poststim_reduced,cidx, 'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_200ms(:, i) == j);
        s_400ms_mean(i, j) = mean(s_400ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_400ms_' num2str(i) '.tif']);
    
    % Gamma part earlier in the code because of graphing reasons
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_400ms(:, i) == j);
        % plot that index
        plot(T-1,nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1), 'Color', ColorSet(j,:)); hold on;
        % do stats on the prestim and poststim gamma power zscore changes
        meantrace_400ms(i,j,:) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1);
        baseline_trace_400ms(j,:) = abs(meantrace_400ms(i,j,1:28));
        poststim_trace_400ms(j,:) = abs(meantrace_400ms(i,j,29:57));
        [p_val_400ms(i,j)] = ranksum(baseline_trace_400ms(j,:),poststim_trace_400ms(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_400ms(i,j)) ', s = ' num2str(s_400ms_mean(i,j))];
    end
    
    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim([-.5 1])
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest', 'AutoUpdate', 'off')
    legend('boxoff');
    % set the size of the legend and axes fonts
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T-1,squeeze(meantrace(i,j,:)), 'Color',[.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_400ms_' num2str(i) '.tif']);
    
    clearvars legend_str
end

% put back the constant vertices so they are in the same order again
vertices_by_cluster(changing_vertices_idx, 1 : size(vertices_by_cluster_reduced, 2)) = vertices_by_cluster_reduced;
vertices_by_cluster_200ms(changing_vertices_idx,:) = vertices_by_cluster_reduced_200ms;
vertices_by_cluster_400ms(changing_vertices_idx,:) = vertices_by_cluster_reduced_400ms;

% separate into left and right vertices again (if you aren't doing both
% sides on the left)
if laterality == 0
    vertices_by_cluster_L = vertices_by_cluster(1:size(vertex_values_sum_L, 1),:);
    vertices_by_cluster_R = vertices_by_cluster(size(vertex_values_sum_L, 1)+1:end,:);
    
    vertices_by_cluster_L_200ms = vertices_by_cluster_200ms(1:size(vertex_values_sum_L, 1),:);
    vertices_by_cluster_R_200ms = vertices_by_cluster_200ms(size(vertex_values_sum_L, 1)+1:end,:);
    
    vertices_by_cluster_L_400ms = vertices_by_cluster_400ms(1:size(vertex_values_sum_L, 1),:);
    vertices_by_cluster_R_400ms = vertices_by_cluster_400ms(size(vertex_values_sum_L, 1)+1:end,:);
end

% graph clusters onto brain maps
for side_index = start_side : sides_num
    
    load ('cmap6to8WRX.mat');
    
    % indicate the figure of the left or right side of the brain
    figure(figHandle(side_index));
    
    % NOTE: now the frame index is actually the number of clusters in each
    % frame rather a frame in terms of a time point
    for frame_index = min_clusters : max_clusters
        
        % overlay the colors
        disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status
        vertex_values_sum = [];
        if laterality == 3
            vertex_values_sum = vertices_by_cluster(:,frame_index);
        elseif side_index == 1 && laterality == 0
            vertex_values_sum = vertices_by_cluster_L(:,frame_index);
        elseif side_index == 2 && laterality == 0
            vertex_values_sum = vertices_by_cluster_R(:,frame_index);
        end
        colormap(ColorSetMaps)
        caxis([-0.5 10.5])
        
        if any(electrodes_present(side_index,:) == 1)
            % for now just overlay the color as number of patients overlapping in
            % this area without any transparency
            h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                'FaceVertexCData', vertex_values_sum, 'FaceColor', 'interp',...
                'LineStyle', 'none', 'FaceAlpha', alpha);
            lightset = [0.6 0.5 0.1];
            material(lightset);
        end
        
        % Save each side of the brain in 4 views
        frames_set(side_index,frame_index,:) = savebrainimages_inflated(side_index, frame_index, views);
        if any(electrodes_present(side_index, :) == 1) && exist('h', 'var')
            % get rid of the color patch so you can lay on the color patch for the
            % next frame
            delete(h)
        end
        
    end
end

% Put the 2 or 4 views into one image
for frame_index = min_clusters : max_clusters
    
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    for side_index = start_side:sides_num
        if side_index == 1
            side = 'L';
        else
            side = 'R';
        end
        for view_index = 1:numel(ViewAngles)
            eval([side '_' ViewAngles{view_index} ' = [];'])
        end
    end
    
    % name the frames for each side
    for side_index = start_side:sides_num
        % Determine what side we are working on
        if(side_index == 1)
            side = 'L';
            ViewAngles = {'Medial'};
        else
            side = 'R';
            ViewAngles = {'Lateral'};
        end
        
        if views == 4
            ViewAngles={'Lateral','Medial','Ventral','Posterior'};
        elseif views == 2
            ViewAngles={'Lateral','Medial'};
        elseif strcmp(views, 'Lateral')
            ViewAngles={'Lateral'};
        elseif strcmp(views, 'Medial')
            ViewAngles={'Medial'};
        end
        
        for view_index=1:length(ViewAngles)
            Frame = frames_set(side_index,frame_index,view_index).cdata;
            % Do cropping here
            Frame(all(all(Frame == 255, 3), 2), :, :) = [];
            Frame(:, all(all(Frame == 255, 3), 1), :) = [];
            eval([side '_' ViewAngles{view_index} ' = Frame;'])
        end
    end
    
    % create the frames for each side
    for side_index = start_side:sides_num
        if side_index == 1
            L_Lateral_Padded = [L_Lateral, ones(size(L_Lateral, 1),...
                max([0, size(L_Medial, 2) - size(L_Lateral, 2)]), 3) * 255];
            
            L_Medial_Padded = [L_Medial, ones(size(L_Medial, 1),...
                max([0, size(L_Lateral, 2) - size(L_Medial, 2)]), 3) * 255];
            
            if laterality ~= 3
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(R_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];
            else
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(L_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];
            end
            
            L_Lateral_Medial = [L_Medial_Padded; L_Lateral_Padded];
            
            if laterality == 0
                R_Posterior_Padded = [R_Posterior, ones(size(R_Posterior, 1),...
                    max([0, size(L_Ventral, 2) - size(R_Posterior, 2)]), 3) * 255];
                
                L_Ventral_R_Posterior = [L_Ventral_Padded; R_Posterior_Padded];
                
                L_Lateral_Medial_Padded = [L_Lateral_Medial; ones(max([0, size(L_Ventral_R_Posterior, 1) - ...
                    size(L_Lateral_Medial, 1)]), size(L_Lateral_Medial, 2), 3) * 255];
                
                L_Ventral_R_Posterior_Padded = [L_Ventral_R_Posterior; ones(max([0, size(L_Lateral_Medial, 1) - ...
                    size(L_Ventral_R_Posterior, 1)]), size(L_Ventral_R_Posterior, 2), 3) * 255];
                
                combined_views_left{frame_index} = [L_Ventral_R_Posterior_Padded, L_Lateral_Medial_Padded];
            elseif laterality == 3
                L_Posterior_Padded = L_Posterior;
                
                L_Ventral_L_Posterior = [L_Ventral_Padded;L_Posterior_Padded];
                
                L_Lateral_Medial_Padded = [L_Lateral_Medial; ones(max([0, size(L_Ventral_L_Posterior, 1) - ...
                    size(L_Lateral_Medial, 1)]), size(L_Lateral_Medial, 2), 3) * 255];
                
                L_Ventral_L_Posterior_Padded = [L_Ventral_L_Posterior; ones(max([0, size(L_Lateral_Medial_Padded, 1) - ...
                    size(L_Ventral_L_Posterior, 1)]), size(L_Ventral_L_Posterior, 2), 3) * 255];
                
                combined_views_left{frame_index} = [L_Ventral_L_Posterior_Padded, L_Lateral_Medial_Padded];
            end
            
            figure
            imshow(combined_views_left{frame_index});
            set(gca,'position',[0 0 1 1],'units','normalized')
            
            % adding time stamp text to the figure
            if laterality == 3 && ismember(overlay,[2 3])
                load(xtick_spacing_filename)
                time_txt = ['Time = ' num2str(round((T(frame_index)-1)*1000)) 'ms'];
                text(400,900,time_txt,'Fontsize',24)
            end
            fn = ['combined_viewsL' '_' num2str(frame_index)];
            print(fn, '-dtiff');
            close
            
        elseif side_index == 2
            
            R_Lateral_Padded = [R_Lateral, ones(size(R_Lateral, 1),...
                max([0, size(R_Medial, 2) - size(R_Lateral, 2)]), 3) * 255];
            
            R_Medial_Padded = [R_Medial, ones(size(R_Medial, 1),...
                max([0, size(R_Lateral, 2) - size(R_Medial, 2)]), 3) * 255];
            
            R_Lateral_Medial = [R_Medial_Padded; R_Lateral_Padded];
            
            
            R_Ventral_Padded = [R_Ventral, ones(size(R_Ventral, 1),...
                max([0, size(L_Posterior, 2) - size(R_Ventral, 2)]), 3) * 255];
            
            if laterality == 0
                L_Posterior_Padded = [L_Posterior, ones(size(L_Posterior, 1),...
                    max([0, size(R_Ventral, 2) - size(L_Posterior, 2)]), 3) * 255];
                
                R_Ventral_L_Posterior = [R_Ventral_Padded; L_Posterior_Padded];
                
                R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(max([0, size(R_Ventral_L_Posterior, 1) - ...
                    size(R_Lateral_Medial, 1)]), size(R_Lateral_Medial, 2), 3) * 255];
                
                R_Ventral_L_Posterior_Padded = [R_Ventral_L_Posterior; ones(max([0, size(R_Lateral_Medial, 1) - ...
                    size(R_Ventral_L_Posterior, 1)]), size(R_Ventral_L_Posterior, 2), 3) * 255];
                
                combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_L_Posterior_Padded];
            elseif laterality == 2
                R_Posterior_Padded = R_Posterior; %(80:495,155:395,:);
                
                R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(116,545,3)*255];
                
                R_Ventral_R_Posterior = [R_Ventral_Padded;R_Posterior_Padded];
                
                combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_R_Posterior];
                combined_views_right{frame_index} = [combined_views_right{frame_index}; 255*ones(15, 776, 3)];
            end
            
            figure
            imshow(combined_views_right{frame_index});
            set(gca,'position',[0 0 1 1],'units','normalized')
            fn = ['combined_viewsR' '_' num2str(frame_index)];
            saveas(gcf,[fn '.tif'])
            close
        end
    end
    
    if laterality == 0
        figure
        combined_views = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left{frame_index}, 1) - ...
            size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
            [combined_views_left{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
            size(combined_views_left{frame_index}, 1)]), size(combined_views_left{frame_index}, 2), 3) * 255]];
        imshow(combined_views);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis tight
        frame(frame_index,:,:,:) = combined_views;
        if exist('j','var')
            fn = ['combined_views_full_' num2str(frame_index) '_' num2str(j)];
        else
            fn = ['combined_views_full_' num2str(frame_index)];
        end
        
        if ismember(overlay,[2 3])
            load(xtick_spacing_filename)
            time_txt = ['Time = ' num2str(round((T(frame_index)-1)*1000)) 'ms'];
            if views == 4
                uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                    'FontSize', 48, 'String', time_txt);
                text(size(combined_views, 2) - 1.5 * uc.Extent(3),...
                    size(combined_views, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 48)
            end
        elseif overlay == 1
            T = [100 200 300 400 500 600 700 800 900 1000];
            time_txt = ['Time = ' num2str(T(frame_index)) ' ms'];
            text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
        end
        
        saveas(gcf,[fn '.tiff'])
    end
end

%% Save data
savename = 'vertex_values_bilat_23pts57bins_recentered_rejoutliers_sounds_restricted_CP_averagePatientMean';
if laterality == 0
    save([savename, '.mat'], 'changing_vertices_idx', 'p_val', 'p_val_200ms',...
        'p_val_400ms', 'poststim', 'poststim_reduced', 'poststim_reduced_200ms', 'poststim_reduced_400ms', 's', 's_200ms', 's_200ms_mean',...
        's_400ms', 's_400ms_mean', 's_mean', 'vertex_valuesL', 'vertex_valuesR', 'vertex_values_sum_L', 'vertex_values_sum_R',...
        'vertex_values_sum_bilateral', 'vertex_values_sum_bilateral_reduced', 'vertices_by_cluster', 'vertices_by_cluster_L',...
        'vertices_by_cluster_R', 'vertices_by_cluster_reduced', 'vertices_by_cluster_reduced_200ms', 'vertices_by_cluster_reduced_400ms')
elseif laterality == 3
    save([savename, '_flipX.mat'], 'changing_vertices_idx', 'p_val', 'p_val_200ms',...
        'p_val_400ms', 'poststim', 'poststim_reduced', 'poststim_reduced_200ms', 'poststim_reduced_400ms', 's', 's_200ms', 's_200ms_mean',...
        's_400ms', 's_400ms_mean', 's_mean', 'vertex_valuesL',...
        'vertex_values_sum', 'vertex_values_sum_reduced', 'vertices_by_cluster',...
        'vertices_by_cluster_reduced', 'vertices_by_cluster_reduced_200ms', 'vertices_by_cluster_reduced_400ms')
end




