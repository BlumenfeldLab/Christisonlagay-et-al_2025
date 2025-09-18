function electrodeLocationMap_fs(data_location, patients, createMontage, laterality, inflationstep, pre_proj_plot)
%{
Chris Micek Edits

  A trimmed-down version of 'displayElectrodesInflated'. Creates
  interactive figures for the left and right hemispheres of the inflated
  surface with the electrodes for each participant in 'patients' plotted.
  Clicking on an electrode contact (when 3D rotation is disabled) displays
  its associated participant ID and channel label.

INPUTS:
  REQUIRED:
    - data_location: Location of patient data folders for brain surface
                     overlays. On Chris's computer, this is 'E:\Chris
                     HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP
                     movies'.

    - patients: Cell array of patient ID codes. There must be a corresponding
                folder in 'data_location' for each.

   - createMontage: Should be 1 if a never-before analyzed patient is
                    included in 'patients':
                      1 = create each patient's montage and map for EEG
                          channels (i.e. read in and save the electrode names
                          and locations from their 'Map... .xlsx' file)
                      0 = nothing

    - laterality: Which hemisphere(s) to plot overlay data onto:
                    0 = bilateral
                    1 = left only
                    2 = right only
                    3 = project both left and right electrodes on the left
                        brain (flipX)

    - inflationstep: How inflated the brain surface should be; options range
                     from 1 (least inflated) to 6 (most inflated)

    - pre_proj_plot: Boolean indicating whether to additionally plot
                     electrode locations prior to their projection onto the
                     brain surface.
%}

image_storage_folder = [data_location '/Group/Wendy picture storage2'];
mni2fs_dir = 'E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master';
mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;
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

%% Generating all frames for each side of the brain
for side_index = start_side : sides_num    % Generate frames for L and R
    
    % Determine what side we are working on
    if side_index == 1
        side = 'L';
        hem = 'lh';
    else
        side = 'R';
        hem = 'rh';
    end
    
    %% hard coded number of frames for each overlay setting, need to change later
    total_num_frames = 1;
   
    
    %------------------------------------------------------------------
    %   Set-up some local variables
    %------------------------------------------------------------------
    
    surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);
    fs_surf(side_index) = export(gifti(surf_fn));
    v_num = size(fs_surf(side_index).vertices, 1);
    %fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] *...
    %   Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
    %fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
    
    surfrender_fn = fullfile(mni2fs_dir,['/surf/' hem '.inflated' num2str(inflationstep) '.surf.gii']);
    inflated_surf(side_index) = export(gifti(surfrender_fn));
    %inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] *...
    %    Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
    %inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);
    
    %     cortex_filename{side_index}  =    [side '_cortex.txt'];
    %     smooth_filename{side_index}  =    [side '_smooth.txt'];
    %     cd([data_location '\' VTK_folder_name]);
    %     cortex_v_num(side_index)     =    eval(['metadata.' side '_cortex_v_num']);
    %     cortex_f_num(side_index)     =    eval(['metadata.' side '_cortex_f_num']);
    %     smooth_v_num(side_index)     =    eval(['metadata.' side '_smooth_v_num']);
    %     smooth_f_num(side_index)    =    eval(['metadata.' side '_smooth_f_num']);
    %
    %     % Input and interpret the vertice data for the cortex and smooth
    %     % surfces
    %     [cortex(side_index).vertices, cortex(side_index).faces] = FormatMeshSurface(cortex_v_num(side_index), cortex_f_num(side_index), data_location, cortex_filename{side_index});
    %     [smooth(side_index).vertices, smooth(side_index).faces] = FormatMeshSurface(smooth_v_num(side_index), smooth_f_num(side_index), data_location, smooth_filename{side_index});
    
    
    %% create a figure of the left or right side of the brain
    % you will be working with this for different electrodes, overlays,
    % views, and light
    
    figHandle(side_index) = figure('Position',[70,70,700,700]); % Original = [50,50,600,600]Creates a figure graphic object
    %     g=trisurf(cortex.faces,cortex.vertices(:,1),cortex.vertices(:,2),cortex.vertices(:,3));
    
    uc = uicontrol('Style', 'text', 'FontName', 'Helvetica', 'Units', 'normalized',...
                    'Position', [0.12, 0.93, 0.9, 0.05], 'FontSize', 20, 'BackgroundColor', 'w');
    
    temp_surf = [];
    temp_surf.hem = hem; % choose the hemesphere 'lh' or 'rh'
    temp_surf.inflationstep = inflationstep; % 1 no inflation, 6 fully inflated
    temp_surf.decimation = 0;
    temp_surf = mni2fs_brain(temp_surf);
    set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
    
%     g = trisurf(inflated_surf(side_index).faces,inflated_surf(side_index).vertices(:,1),...
%         inflated_surf(side_index).vertices(:,2),inflated_surf(side_index).vertices(:,3));
%     
%     lighting flat;
%     set(g, 'FaceColor', [1 1 1], 'LineStyle', 'none')%, 'EdgeLighting','gouraud')%, 'FaceLighting', 'gouraud');
    set(gcf, 'color', [1 1 1]);
    hold on
    %     xlabel('X');
    %     ylabel('Y');
    %     zlabel('Z');
    
    axis off;
    axis equal;

    electrode_mapping_times = 1;

    for e = 1:electrode_mapping_times
        scatter_count = 0;
        scatter_list = gobjects(2, length(patients));
        e_verts = cell(1, length(patients));
        side_label_inds = cell(1, length(patients));
        labels = cell(1, length(patients));
        for p = 1:length(patients)
            
            %% Get electrode info for the patient
            
            % this part plots the electrodes for each patient and the overlay data
            % for each patient
            patient_folder = [data_location, '/', patients{p}];
            
            %% Make the Montage
            % if the left/right electrode montages have not been created, do so
            % now, indicated by createMontage = 0 or = 1
            if createMontage == 1
                % in this case, just do the normal thing: create the normal
                % (non-X-flipped) montage for left and right
                flipX = 0;
                create_electrode_montage(patient_folder,flipX)
                if  laterality == 3 % for x-flipped in addition to regular
                    % in this case, create the x-flipped montage for the left and
                    % right electrodes but you will later only use the RIGHT
                    % electrodes with the coordinates flipped to the LEFT
                    flipX = 1;
                    create_electrode_montage(patient_folder,flipX)
                end
            end
            
            %         % Combine EEG Channel with Electrode Position
            %         cd(data_location)
            %
            %         [overlay_frame_data] = OrganizeElectrodePosition(data_location, patient_folder, p, overlay, createMontage, overlay_frame_data);
            
            display(patients{p})
            
            % Now load some data files needed for electrode placement and surface shading
            cd([data_location, '/', patients{p}]);
            load([side '_MontageMap.mat']);
            electrode = eval([side '_MontageMap' ]);
            
            if ~isempty(electrode)
                % convert to fs space
                electrode_temp = [electrode(:, 2 : 4), ones(size(electrode, 1), 1)] * ((Tfstovox_rcor' * Trsvoxtomni_rcor' * (Tmni' ^ -1)) ^ -1);
                electrode(:, 2 : 4) = electrode_temp(:, 1 : 3);
            end
            
            %         electrode_num = size(electrode,1);
            % load the X-flipped montage and add to list of electrodes
            if laterality == 3
                load('R_MontageMap_flipX.mat');
                electrode = cat(1,electrode,R_MontageMap);
            end
            
            %% Plot patient's electrodes
            
%             electrode = [electrode(:, 1), [electrode(:, 2:end), ones(size(electrode, 1), 1)] * Tmni'];
%             electrode = electrode(:, 1:4);
            
            % assign this patient's electrodes one of the colors
            colorElectrode = colors(p, :);
            
            side_label_inds{p} = electrode;
            
            if isempty(electrode) % Some files only have electrodes on the one side
                disp('-- No electrode data!');
                electrodes_present(side_index,p) = 0;
            else
                all_electrodes = electrode(:,1);
                electrodes_present(side_index,p) = 1;
                
                if pre_proj_plot
                    scatter3(electrode(:, 2), electrode(:, 3), electrode(:, 4), 40, 'k')
                end
                
                %                 [elecVertices2]= ProjectElectrode2TransSurf(smooth_v_num, smooth(side_index).vertices, all_electrodes, electrode, side_index);
                %                 electrode(:,2:4)                       = smooth(side_index).vertices(elecVertices2,1:3);
                [elecVertices] = ProjectElectrode2TransSurf(v_num, fs_surf(side_index).vertices, electrode);
                
                e_verts{p} = elecVertices;
                
                %Plot the electrodes in a specified color and count them
                %in case you need to remove them later
                try
                    p_data = load([patient_folder, filesep, 'labels_all_gray_matter.mat']);
                catch
                    p_data = load([patient_folder, filesep, 'labels_first_surgery_all_gray_matter.mat']);
                end
                
                labels{p} = p_data.labels;
                
                scatter3(inflated_surf(side_index).vertices(elecVertices,1),...
                    inflated_surf(side_index).vertices(elecVertices,2),...
                    inflated_surf(side_index).vertices(elecVertices,3),40, colorElectrode,'filled');
                
                point_set = scatter3(1, 1, 1, 110, 'r', 'filled', 'MarkerFaceAlpha', 0,'MarkerEdgeAlpha', 0);
                set(point_set, 'XData', [], 'YData', [], 'ZData', [])
                
                scatter_list(2, p) = point_set;
                
                scatter_list(1, p) = scatter3(inflated_surf(side_index).vertices(elecVertices,1),...
                    inflated_surf(side_index).vertices(elecVertices,2),...
                    inflated_surf(side_index).vertices(elecVertices,3), 110, 'r','filled',...
                    'MarkerFaceAlpha', 0,'MarkerEdgeAlpha', 0);
                
                scatter_count = scatter_count + 1;
            end
            
        end
        
        for p = 1:length(patients)
            try
                scatter_list(1, p).ButtonDownFcn = {@display_name, scatter_list(2, :), uc, p, patients{p}, labels{p},...
                    [side_label_inds{p}(:, 1), inflated_surf(side_index).vertices(e_verts{p},1),...
                    inflated_surf(side_index).vertices(e_verts{p},2),...
                    inflated_surf(side_index).vertices(e_verts{p},3)]};
            catch
            end
        end
    end
    %------------------------------------------------------------------
    %   Frame Generation with overlaying color onto the brain surface
    %------------------------------------------------------------------

    % go to the correct folder
    cd(image_storage_folder)
    load('cmap6to8WRX.mat');

    for frame_index = 1:total_num_frames

        %% overlay the colors

        disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status

        lightset = [0.6 0.5 0.1];
        material(lightset);

    end

end

end



%==================================================================
%
%Title:         ProjectElectrode2TransSurf
%Description:
%
%Inputs:
%==================================================================

% Project electrodes to the surface
%Input: number of vertices of the projected surface, matrix of vertices, all
%electrodes, electrodes that is going to be projected, Powerdata, frame index
%Output:electtrode vertices on the new surface, the color/power value of that electrode
function [elecVertices]=ProjectElectrode2TransSurf(vertices_num, vertices, electrode)
elecVertices = zeros(size(electrode,1),1);
for elec_index = 1:size(electrode,1)
    xCoord = ones(vertices_num,1)*electrode(elec_index, 2);
    yCoord = ones(vertices_num,1)*electrode(elec_index, 3);
    zCoord = ones(vertices_num,1)*electrode(elec_index, 4);
    dists = sqrt((vertices(:,1) - xCoord).^2 + (vertices(:,2) - yCoord).^2 + (vertices(:,3) - zCoord).^2);
    [minVal, minInd] = min(dists);
    elecVertices(elec_index, 1) = minInd;
end
% load('labels');
% overlay_frame_data = overlay_frame_data';
% electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);
end


%==================================================================
%
%Title:         create_electrode_montage
%Description:   This function takes a comma delimited sheet with all the
%               names of the electrodes and coordinates and reads it into
%               Matlab. It then saves the electrodes according to left and
%               right side with the (x,y,z) coordinates and the index of
%               that electrode in 'labels'. Essentially the electrode names
%               in 'labels' comes from the EEG .edf file and the coordinates
%               come from the .mgrid file so this reconciles the naming of
%               electrodes between the two
%
%Inputs:        "patient_folder"    [str] location of this patient's files
%               "flipX"             [int] 0 means it's the regular montage
%                                         1 means use the X-flipped montage
%==================================================================

function create_electrode_montage(patient_folder,flipX)
cd (patient_folder);

% if you are creating a montage for the x-flipped electrodes, you want to
% open and then save everything with a "_flipX" attached to it
if flipX == 0
    flipX_str = '';
elseif flipX == 1
    flipX_str = '_flipX';
end

% Read in Electrode locations (Map.xls or Map_flipX.xls) and Montage file (Montage.xls).
% load('labels.mat');
try
    load('labels_all_gray_matter.mat');
catch
    load('labels_first_surgery_all_gray_matter.mat');
end
if size(labels, 1) == 1
    labels = labels';
end
b1 = linspace(1, length(labels), length(labels))';
b1 = num2cell(b1);
b2 = labels(:,1);
montage = [b1 b2];
% [a1, a2, map]=xlsread('Map.xls');
[a1, a2, map]=xlsread('Map_all_gray_matter.xlsx');
% Convert the long channel names (i.e "A_L_Most_Mid_Frontal_Polar_12" into the abbreviated names found in the montage (i.e. "A12")
try
    for i=1:size(map,1)
        str             =   map{i,1};
        str_delim       =   regexp(str,'_','split');
        if length(str_delim)==1 || length(str_delim) == 2
            if startsWith(str_delim{1,1}, {'L', 'R'})
                map{i,5}=strtrim(str_delim{1,1}(1));
                map{i,1}=strtrim(str_delim{1,2});
                map{i,1}(map{i,1} == ';') = [];
            else
                map{i,5}=strtrim(str_delim{1,1}(end));
            end
        else
            firstLetter     =   str_delim{1};
            positionLetter  =   strtrim(str_delim{2});
            digitStr        =   str(isstrprop(str,'digit'));
            replaceStr      =   strcat(firstLetter, digitStr);
            map{i,1}        =   replaceStr; % i.e. "A12"
            map{i,5}        =   positionLetter; % Either L or R
        end
    end
catch
    disp(['Error reading Map file! There was a problem on line ' num2str(i) '. Verify that file is in the correct format (i.e each row is Letter_Side*Number. Unacceptable entries: J_Sup_Parietal8,PEG1)']);
    error('MATLAB:myCode:dimensions', '');
end

Montage=montage;
%%
%     try
% Combine montage with 3-d electrode locations
for i=1:size(montage,1)
    MontageMap(i,1)=montage{i,1}; % Abbreviated Channel ("A12")
    for j=1:size(map,1)
        if strcmp(strtrim(montage{i,2}), map{j,1})
            if i == 113
                leah = [];
            end
            MontageMap(i,2)=map{j,2}; % X Coordinates
            MontageMap(i,3)=map{j,3}; % Y Coordinates
            MontageMap(i,4)=map{j,4}; % Z Coordinates
            %                     MontageMap(i,5) = newTimes{i,1}; % Arrival Time
            Position(i,1)=map{j,5}; % Side (R or L)
        end
    end
end
%         % catch
%         %     disp(['Error on MontageMap write! There was a problem on line ' num2str(i) '. Verify that file is in the correct format (i.e each row is Letter_Side*Number. Unacceptable entries: J_Sup_Parietal8)']);
%         %     error('MATLAB:myCode:dimensions', '');
%     end
j=1;k=1;

% check MontageMap to see if there are any errors that result in
% (0,0,0) coordinates for (x,y,z)
for i = 1:size(MontageMap,2)
    if MontageMap(i,1) == 0 && MontageMap(i,2) == 0 && MontageMap(i,3) == 0
        display(MontageMap{i,1});
        error('This electrode name in Map or Map_flipX is not matching, double check the original Map sheet for misnamed electrodes compared to labels');
    end
end


% Break the montage into L and R sides
for i=1: length(Position)
    if (Position(i)=='L')
        L_MontageMap(j,1:4)=MontageMap(i,1:4);
        j=j+1;
    else
        R_MontageMap(k,1:4)=MontageMap(i,1:4);
        k=k+1;
    end
end

if (j==1)
    L_MontageMap=[];
elseif (k==1)
    R_MontageMap=[];
end
% Save data into mat files
MontageMap=[L_MontageMap; R_MontageMap];
MontageMap=sortrows(MontageMap,1);
save(['Montage' flipX_str '.mat'], 'Montage');
save(['MontageMap' flipX_str '.mat'], 'MontageMap');
save(['L_MontageMap' flipX_str '.mat'], 'L_MontageMap');
save(['R_MontageMap' flipX_str '.mat'], 'R_MontageMap');

end

function display_name(~, event, point_set, uc, patient_ind, patient_name, labels, MontageMap)

name = labels{MontageMap(ismember(MontageMap(:, 2:end), event.IntersectionPoint, 'rows'), 1)};

for pset = 1:length(point_set)
    try
        set(point_set(pset), 'MarkerFaceAlpha', 0,'MarkerEdgeAlpha', 0)
    catch
    end
end

set(point_set(patient_ind), 'XData', event.IntersectionPoint(1), 'YData', event.IntersectionPoint(2),...
    'ZData', event.IntersectionPoint(3), 'MarkerFaceAlpha', 1,'MarkerEdgeAlpha', 1)

uc.String = ['Subject ', patient_name, ', Electrode ', name];

end


