
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



