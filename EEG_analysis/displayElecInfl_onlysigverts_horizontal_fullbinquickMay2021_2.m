
set(gcf, 'Renderer', 'OpenGL')

addpath('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\MultipleTestingToolbox\MultipleTestingToolbox')
addpath('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\Mass_Univariate_ERP_Toolbox-master')
mappeddatadrive='V:';
time_offset=2;
lengthofbaseline=[500];  %%baseline 500 or 1000 ms prestim
flip=[0]; %%bilater=0; flipX=3;

%%laterality: not flipped
laterality=0;
flipxmuliplier=1;
laterality_label=[''];
display_color_bar=0;

%% common baseline
baseline_name='common';
baselinetype='_common_r';

%% baseline start
baseline_start_real=59;
baseline_start=19;
baselinesuffixes='baseline5to0';

%% include baseline in stats
statsinclusion=['_baselinestats_'];
statsinclusionverbose=['including baseline statistics'];
analysis_start=41;

statstype=2;
%%
if statstype==1
    statstiming='early';
    finalbin=97;
    total_num_frames2=length(baseline_start_real:finalbin);
    
    
else
    statstiming='late';
                    finalbin=137;
    total_num_frames2=length(baseline_start_real:finalbin);
    
end

%% otherstuff
baseline_end_real=78;
baseline_end=38;

xtick_spacing_filename = ['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\filter_times_win100ms_shift25ms_full.mat'];
load (xtick_spacing_filename)
final_bin=length(T); %first 500 ms
runduration=[num2str(4000) 'ms'];

%% number of parcels
nparcels=[ 80 ];



addpath('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\MultipleTestingToolbox\MultipleTestingToolbox')
addpath('D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\Mass_Univariate_ERP_Toolbox-master')
mappeddatadrive='V:';

lengthofbaseline=[500];  %%baseline 500 or 1000 ms prestim
flip=[0]; %%bilater=0; flipX=3;

%%laterality: not flipped
laterality=0;
flipxmuliplier=1;
laterality_label=[''];

%% common baseline
baseline_name='common';
baselinetype='_common_r';

%% baseline start
baseline_start_real=59;
baseline_start=19;
baselinesuffixes='baseline5to0';

%% include baseline in stats
statsinclusion=['_baselinestats_'];
statsinclusionverbose=['including baseline statistics'];
analysis_start=41;

%% otherstuff
baseline_end_real=78;
baseline_end=38;

load (xtick_spacing_filename)
final_bin=length(T); %first 500 ms
timeline_stuff='late';

%% number of parcels
nparcels=[ 80 ];


for overlays=[3]
    if isequal(overlays, 1)
        overlay='CP';
    elseif isequal(overlays, 2)
        overlay='CNP';
    else
        overlay='Subtraction';
    end
    
    inflationstep=5;
    createMontage = 0;
    views = 4;
    plot_electrodes=[];
    parcel_flavor=['kmeans_mask_refitted']
    
    data_location=[mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies'];
    data_location2=[mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior'];
    save_location_1 = ['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\'];
    save_location_2_1= [num2str(lengthofbaseline) 'ms_' baseline_name '_baseline_' statsinclusion];
    save_location_2_2=['100msbins75msoverlap'];
    save_location_2_3=['filtered_2to2s'];
    save_location_2=[save_location_2_1 save_location_2_3];
    save_location=[save_location_1 save_location_2_2 '\' save_location_2  '\' num2str(nparcels) '_parcels_' parcel_flavor  '\' laterality_label overlay runduration];
    
    if ~exist(save_location, 'dir')
        mkdir(save_location)
    end
    
    verticesdata=['\\172.23.254.106\Data25\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 157 Bins Sounds Restricted by CP Accuracy 100ms\filter_baseline5to0_common_r\Gamma_power_' overlay '\vertex_values.mat'];
    
        load(verticesdata)
    
    rights=[]; lefts=[];
    if isfile([save_location '\' statstiming 'significantvertices_' save_location_2 '_' num2str(nparcels) 'parcels.mat'])
        disp(['Loading significant vertices for ' laterality_label ' ' overlay ' using with ' num2str(nparcels) ' parcels, with a ' num2str(lengthofbaseline) 'ms ' baseline_name ' baseline '  statsinclusionverbose])
        
        load([save_location '\' statstiming 'significantvertices_' save_location_2 '_' num2str(nparcels) 'parcels.mat'])
    else
        disp(['Calculating significant vertices for  ' laterality_label ' ' overlay ' using with ' num2str(nparcels) ' parcels, with a ' num2str(lengthofbaseline) 'ms ' baseline_name ' baseline '  statsinclusionverbose])
        
        for a=1:100
            
            if isequal(a/10, ceil(a/10))
                disp(['Loading map ' num2str(a) ])
            end
            
            load ([save_location '\' statstiming 'significant_vertices_map' num2str(a) '.mat'])
            
            rights=cat(3, rights, significantvertex_valuesR);
            lefts=cat(3, lefts, significantvertex_valuesL);
            
        end
        
        
        fraction_significant_right=sum(rights,3)/size(rights,3);
        sign_right=ge(fraction_significant_right,.5);
        
        fraction_significant_left=sum(lefts,3)/size(lefts,3);
        sign_left=ge(fraction_significant_left,.5);
        
        %                             if isequal(signal_on, 39)
        %                                 sign_left2=zeros(size(sign_left,1),38);
        %                                 sign_left=[sign_left2, sign_left];
        %
        %                                 sign_right2=zeros(size(sign_right,1),38);
        %                                 sign_right=[sign_right2, sign_right];
        %
        %                             end
        
        save([save_location '\' statstiming 'significantvertices_' save_location_2 '_' num2str(nparcels) 'parcels.mat'], 'sign_left', 'sign_right');
        
        
        
    end
    
    
    
    %%%create file for later
    daterun = date;
    save([save_location '\' statstiming 'vertexsum_' save_location_2 '_' num2str(nparcels) 'parcels.mat'], 'daterun' )
                      
    
    overlay_labels = {'Electrode Distribution', 'Arrival Times', 'CP', 'CNP', 'Electrode Density', 'Subtraction', 'Comparison'};
    
    % If overlay numeric code, replace with the label
    if isnumeric(overlay)
        overlay = overlay_labels{overlay};
    end
    
    %% Define directory path to which images will be saved
    if laterality == 3; save_location = [save_location, '/flipX']; end
    
    new_color=[''];
    
    suffix='_horizontal3views';
    
    
    switch overlay
        case 'Electrode Distribution'
            image_storage_folder = [save_location '/Electrode_Distribution'];
        case 'CP'
            image_storage_folder = [save_location '/' statstiming 'Gamma_power_CP' suffix new_color];
        case 'CNP'
            image_storage_folder = [save_location '/' statstiming 'Gamma_power_CNP' suffix new_color];
        case 'Electrode Density'
            image_storage_folder = [save_location '/Electrode_Density'];
        case 'Subtraction'
            image_storage_folder = [save_location '/' statstiming 'Gamma_power_Subtraction' suffix new_color];
        case 'Comparison'
            image_storage_folder = [save_location '/Gamma_power_Comparison'];
    end
    
    if ~exist(image_storage_folder, 'dir')
        mkdir(image_storage_folder)
    end
    
    % Load the mni template and transform
    mni2fs_dir = 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\mni2fs-master';
    mni_nifti_path = fullfile(mni2fs_dir, 'MNI_T1_1mm_stripped.nii');
    load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
    mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
    Tmni = mnit1_1mm_stripped.transform;
    
    
    
    % Define the face transparency value
    alpha = 0.7;
    
    % Side laterality - which side of the brain to draw
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
    
    % Create figure
    figHandle = gobjects(1, sides_num);
    
    total_num_frames = numel(T);
    
    
    for side_index=1:2
        if side_index == 1
            vertex_values = vertex_valuesL;
            side = 'L';
            hem = 'lh';
        elseif side_index == 2
            vertex_values = vertex_valuesR;
            side = 'R';
            hem = 'rh';
        end
        electrodes_present(side_index,:) = 1;
        
        
        %%%set up brain!
        
        % Setup surface variables
        surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);  %%need directory
        fs_surf(side_index) = export(gifti(surf_fn));  %%we need here through line 472 for parcellation
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
        figHandle(side_index) = figure('Position', [70, 70, 700, 700]);
        
        % Initialize a temporary surface variable
        temp_surf = [];
        temp_surf.hem = hem; % choose the hemesphere 'lh' or 'rh'
        temp_surf.inflationstep = inflationstep; % 1 no inflation, 6 fully inflated
        temp_surf.decimation = 0;
        temp_surf = mni2fs_brain(temp_surf);
        set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
        set(gcf, 'color', [1 1 1]);
        hold on
        axis off;
        axis equal;
        
        %         figure('Color','k','position',[20 72 800 600])
        %         hold on
        % Load and Render the FreeSurfer surface
        S = [];
        
        
        
        %
        %
        %         mni2fs_dir = [mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\Functions for Inflated Brain Display\mni2fs-master'];
        %         mni_nifti_path = fullfile(mni2fs_dir, '..', 'MNI_T1_1mm_stripped.nii');
        %         mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
        %         S.hem = hem; % choose the hemesphere 'lh' or 'rh'
        %         S.inflationstep = 5; % 1 no inflation, 6 fully inflated
        %         S.plotsurf = 'inflated';
        %         S.lookupsurf = 'mid';
        %         S.decimation = false; % Decimate the surface for speed. (Use FALSE for publishable quality figures).
        %         S = mni2fs_brain(S);
        %
        hold on;
        
        
        % Setup surface variables
        surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);  %%need directory
        fs_surf(side_index) = export(gifti(surf_fn));  %%we need here through line 472 for parcellation
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
        figHandle(side_index) = figure('Position', [70, 70, 700, 700]);
        
        % Initialize a temporary surface variable
        temp_surf = [];
        temp_surf.hem = hem; % choose the hemesphere 'lh' or 'rh'
        temp_surf.inflationstep = inflationstep; % 1 no inflation, 6 fully inflated
        temp_surf.decimation = 0;
        temp_surf = mni2fs_brain(temp_surf);
        set(temp_surf.p, 'Faces', inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices)
        set(gcf, 'color', [1 1 1]);
        hold on
        axis off;
        axis equal;
        
        if ~strcmp(overlay, 'Arrival Times')
            
            % go to the correct folder
            color_map_outer_bounds = [-15 15];
            color_map_inner_bounds = [-2 2];
            %                         color_map = [mappeddatadrive ':\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\cmap6to8WRX.mat'];
            %                         color_map = 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\verysmallinnerboundscolor.mat';
            color_map = 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\extendedcolorbar.mat';
            
            cd(image_storage_folder)
            load(color_map);
            
            if ~any(strcmp(overlay, {'Electrode Density', 'Comparison'}))
                
                vertex_values_sum = squeeze(nansum(vertex_values, 2));
                divisor = squeeze(sum(vertex_values ~= 0, 2));
                divisor(divisor == 0) = 1;
                
                % Take weighted mean across patients
                vertex_values_sum = vertex_values_sum ./ sqrt(divisor);
                
                % load the colormap and set axis of colormap
                colormap(mycmap)
                caxis(color_map_outer_bounds)
                
            end
            
            
            
            if isequal(side_index,1)
                leftverts=vertex_values_sum;
                save([save_location '\' statstiming 'vertexsum_' save_location_2 '_' num2str(nparcels) 'parcels.mat'], 'leftverts', '-append');
            elseif isequal(side_index,2)
                rightverts=vertex_values_sum;
                save([save_location '\' statstiming 'vertexsum_' save_location_2 '_' num2str(nparcels) 'parcels.mat'], 'rightverts', '-append');
            end
            
            
            
            for frame_index = baseline_start_real: finalbin
                frame_index2=frame_index-baseline_start_real+1;
                %% overlay the colors
                disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status
                
                % sum the count of vertices across patients
                if (any(strcmp(overlay, {'CP', 'CNP', 'Subtraction'})) && ndims(vertex_values) == 3) || ...
                        any(strcmp(overlay, {'Electrode Density', 'Comparison'}))% the surface is colored according to zscores of power
                    
                    
                    if isequal(side_index,1)
                        sign_vertices=sign_left;
                    elseif isequal(side_index,2)
                        sign_vertices=sign_right;
                    end
                    
                    h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                        'FaceVertexCData',(vertex_values_sum(:, frame_index).*sign_vertices(:,frame_index2)),'FaceColor','interp',...
                        'LineStyle', 'none', 'FaceAlpha', 'flat');
                    
                    alpha_data = ones(numel(h.FaceVertexCData), 1) * alpha;
                    %
                    %
                    withinbounds=find(h.FaceVertexCData >= color_map_inner_bounds(1) & h.FaceVertexCData <= color_map_inner_bounds(2));
                    alpha_data(h.FaceVertexCData >= color_map_inner_bounds(1) & h.FaceVertexCData <= color_map_inner_bounds(2)) = (((abs(h.FaceVertexCData(withinbounds)))/2).^4)*alpha;
                    
                    lightset = [0.6 0.5 0.1];
                    material(lightset);
                    
                    if ~isequal(max(sign_vertices(:,frame_index2)),0)
                        set(h, 'FaceVertexAlphaData', alpha_data)
                    else
                        set(h, 'FaceVertexAlphaData', alpha_data)
                    end
                    
                end
         glah=1;       
                
                %% Save each side of the brain in 3 views
                frames_set(side_index,frame_index, :) = savebrainimages(side_index, frame_index2, views, figHandle);
                if any(electrodes_present(side_index,:) == 1) && exist('h','var')
                    % get rid of the color patch so you can lay on the color patch for the
                    % next frame
                    delete(h)
                end
            end
        end
    end
    
    
    %% Put the 3 views into one image
    
    if strcmp(overlay, 'Arrival Times')
        all_frames_num = electrode_mapping_times;
    else
        all_frames_num = total_num_frames;
    end
    
    for frame_index =  baseline_start_real: finalbin %1 : all_frames_num
        ViewAngles={'Lateral','Medial','Ventral','Posterior'};
        for side_index = start_side:sides_num
            if side_index == 1
                side = 'L';
            else
                side = 'R';
            end
            for view_index = 1 : numel(ViewAngles)
                eval([side '_' ViewAngles{view_index} ' = [];'])
            end
        end
        
        %% name the frames for each side
        for side_index = start_side : sides_num
            
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
        
        %% create the frames for each side
        for side_index = start_side : sides_num
            if side_index == 1
                
                L_Ventral = cat(1, L_Ventral, ones(5, size(L_Ventral, 2), 3) * 255);
                L_Ventral = cat(1, ones(5, size(L_Ventral, 2), 3) * 255, L_Ventral);
                
                L_Lateral = cat(1, L_Lateral, ones(floor((size(L_Ventral, 1) - size(L_Lateral, 1)) / 2), size(L_Lateral, 2), 3) * 255);
                L_Lateral = cat(1, ones(size(L_Ventral, 1) - size(L_Lateral, 1), size(L_Lateral, 2), 3) * 255, L_Lateral);
                
                L_Medial = cat(1, L_Medial, ones(floor((size(L_Ventral, 1) - size(L_Medial, 1))/ 2), size(L_Medial, 2), 3) * 255);
                L_Medial = cat(1, ones(size(L_Ventral, 1) - size(L_Medial, 1), size(L_Medial, 2), 3) * 255, L_Medial);
                
                L_Posterior = cat(1, L_Posterior, ones(floor((size(L_Ventral, 1) - size(L_Posterior, 1))/ 2), size(L_Posterior, 2), 3) * 255);
                L_Posterior = cat(1, ones(size(L_Ventral, 1) - size(L_Posterior, 1), size(L_Posterior, 2), 3) * 255, L_Posterior);
                
                %combined_views_left{frame_index} = cat(2, L_Lateral, L_Medial, L_Ventral, L_Posterior);
                %                                     combined_views_left{frame_index} = cat(2, L_Ventral, L_Lateral, L_Medial, L_Posterior);
                buffer=ones(size(L_Ventral, 1), 5, 3)*255;
                combined_views_left{frame_index} = cat(2, L_Ventral, buffer, L_Lateral, buffer, L_Medial, buffer);
                
                % add a colorbar
                if display_color_bar
                    color_bar_image = imread('L_1_Ventral.tif');
                    color_bar_image = color_bar_image(:, 850 : 1024, :);
                    vert_dim = size(combined_views_left{frame_index}, 1);
                    horz_dim = floor((vert_dim / size(color_bar_image, 1)) * size(color_bar_image, 2));
                    for r = 1 : 3
                        color_bar_image_resize(:, :, r) = imresize(color_bar_image(:, :, r), [vert_dim, horz_dim]);
                    end
                    combined_views_left{frame_index} = [combined_views_left{frame_index}, color_bar_image_resize];
                end
                
                figure
                imshow(combined_views_left{frame_index});
                set(gca, 'position', [0 0 1 1], 'units', 'normalized')
                
                % adding time stamp text to the figure
                
                % save figure
                saveas(gcf,['combined_viewsL' '_' num2str(frame_index), '.tif'])
                close
                
            elseif side_index == 2
                
                R_Ventral = cat(1, R_Ventral, ones(5, size(R_Ventral, 2), 3) * 255);
                R_Ventral = cat(1, ones(5, size(R_Ventral, 2), 3) * 255, R_Ventral);
                
                R_Lateral = cat(1, R_Lateral, ones(floor((size(R_Ventral, 1) - size(R_Lateral, 1)) / 2), size(R_Lateral, 2), 3) * 255);
                R_Lateral = cat(1, ones(size(R_Ventral, 1) - size(R_Lateral, 1), size(R_Lateral, 2), 3) * 255, R_Lateral);
                
                R_Medial = cat(1, R_Medial, ones(floor((size(R_Ventral, 1) - size(R_Medial, 1))/ 2), size(R_Medial, 2), 3) * 255);
                R_Medial = cat(1, ones(size(R_Ventral, 1) - size(R_Medial, 1), size(R_Medial, 2), 3) * 255, R_Medial);
                
                R_Posterior = cat(1, R_Posterior, ones(floor((size(R_Ventral, 1) - size(R_Posterior, 1))/ 2), size(R_Posterior, 2), 3) * 255);
                R_Posterior = cat(1, ones(size(R_Ventral, 1) - size(R_Posterior, 1), size(R_Posterior, 2), 3) * 255, R_Posterior);
                
                %combined_views_right{frame_index} = cat(2, R_Posterior, R_Ventral, R_Medial, R_Lateral);
                %                                     combined_views_right{frame_index} = cat(2, R_Posterior, R_Medial, R_Lateral, R_Ventral);
                
                buffer=ones(size(R_Ventral, 1), 5, 3)*255;

                combined_views_right{frame_index} = cat(2, buffer,R_Medial, buffer, R_Lateral, buffer, R_Ventral);
                
                figure
                imshow(combined_views_right{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                saveas(gcf,['combined_viewsR' '_' num2str(frame_index), '.tif'])
            end
        end
        
        %------------------------------------------------------------------
        %   Movie Assembly Procedure
        %------------------------------------------------------------------
        % Add all 6 components of each movie frame together onto one plot after
        % cropping. Then after assembling the frame, add to an avi movie buffer
        % output
        if laterality == 0
            figure
            combined_views = cat(2, combined_views_right{frame_index}, combined_views_left{frame_index});
            
            imshow(combined_views);
            set(gca,'position',[0 0 1 1],'units','normalized')
            axis tight
            %                                 frame(frame_index,:,:,:) = combined_views;
            if exist('j','var')
                fn = ['combined_views_full_' num2str(frame_index) '_' num2str(j)];
            else
                fn = [overlay '_' num2str(frame_index)];
            end
            
            if any(strcmp(overlay, {'CP', 'CNP', 'Subtraction', 'Comparison'}))
                if strcmp(overlay, 'Comparison') && strcmp(params.CompareMethod, 'XCorr')
                    time_txt = ['Lag = ' num2str(round((T2(frame_index))*1000)) 'ms'];
                else
                    t_num = round((T(frame_index) - time_offset) * 1000);
                    
                    xtick_spacing_filename_1=strfind(xtick_spacing_filename, '\');
                    xtick_spacing_filename_2=xtick_spacing_filename(xtick_spacing_filename_1(end)+1:end);
                    if strcmp(xtick_spacing_filename_2, 'spectrogram_times_57bins.mat')
                        time_txt = [num2str(t_num - 125) ' to ' num2str(t_num + 125) ' ms'];
                    elseif strcmp(xtick_spacing_filename_2, 'spectrogram_times_40bins.mat')
                        time_txt = [num2str(t_num - 25) ' to ' num2str(t_num + 25) ' ms'];
                    elseif strcmp(xtick_spacing_filename_2, 'filter_times_win100ms_shift50ms.mat') ||contains(xtick_spacing_filename_2, 'filter_times_win100ms_shift25ms')
                        time_txt = [num2str(t_num - 50) ' to ' num2str(t_num + 50) ' ms'];
                    end
                end
                if views == 4
                    
                    uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                        'FontSize', 24, 'String', time_txt);
                    
                    text(size(combined_views, 2) - uc.Extent(3) - 250,...
                        size(combined_views, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 18)
                end
            elseif strcmp(overlay, 'Arrival Times')
                T = [100 200 300 400 500 600 700 800 900 1000];
                time_txt = ['Time = ' num2str(T(frame_index)) ' ms'];
                text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
            end
            
            saveas(gcf, [fn '.tif'])
            %                                 close
        end
    end
    
    % if you want to create a montage of evenly space post-stimulus frames,
    % then uncomment this code
    if views == 2
        
        %indices of frames to capture in montage
        frame_idx = [33 37 41 45 49 53 57];
        
        large_image = [];
        for i = 1:length(frame_idx)
            large_image = cat(1,squeeze(frame(frame_idx(i), :, :, :)),large_image);
        end
        figure
        imshow(large_image)
        set(gca,'position',[0 0 1 1],'units','normalized')
        saveas(gcf,'frames_montage.tiff')
        saveas(gcf,'frames_montage.eps')
    end
    close all;
end


function frames_set = savebrainimages(side_index,frame_index,views, figHandle)

figure(figHandle(side_index))
zoomin=1.24; %%used to be 1.25
if views == 4
    % For each frame, take a snapshot for each of these views
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    % These arrays define the viewing and light perspective of the frame
    zoom_factors = [zoomin, zoomin, zoomin, zoomin];
    if side_index == 1
        side = 'L';
        Light_Pos = [1200 -200 500; -700 -800 -100;-200, -100, -800; 200, 800, 0];
        View_Pos = [90 0;-90, 0;0 -90;180 0];
    elseif side_index == 2
        side = 'R';
        Light_Pos = [-1200 -300 700; 1000 -1100 -200; -200, -100, -800;-200, 800, 0];
        View_Pos = [-90 0;90, 0;0 -90;180 0];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif views == 2
    zoom_factors = [zoomin, zoomin];
    if side_index == 1
        side = 'L';
        ViewAngles = {'Lateral', 'Medial'};
        Light_Pos = [1200 -200 500; -700 -800 -100];
        View_Pos = [90 0;-90, 0];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral', 'Medial'};
        Light_Pos = [-1200 -300 700; 1000 -1100 -200];
        View_Pos = [-90 0;90, 0];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif strcmp(views, 'Lateral')
    zoom_factors = zoomin;
    if side_index == 1
        side = 'L';
        ViewAngles = {'Lateral'};
        Light_Pos = [1200 -200 500];
        View_Pos = [90 0];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral'};
        Light_Pos = [-1200 -300 700];
        View_Pos = [-90 0];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif strcmp(views, 'Medial')
    zoom_factors = zoomin;
    if side_index == 1
        side = 'L';
        ViewAngles = {'Medial'};
        Light_Pos = [-700 -800 -100];
        View_Pos = [-90, 0];
        %         Light_Pos = [-700 -800 -300];
        %         View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Medial'};
        Light_Pos = [1000 -1100 -200];
        View_Pos = [90, 0];
    else
        Light_Pos=0;
        View_Pos=0;
    end
end

for i = 1:length(ViewAngles)
    lightHandle = light('Position', Light_Pos(i,:));
    view(View_Pos(i, :));
    if i == 1
        % with this command, matlab will freeze the aspect ratio of
        % your image instead of changing the object size to fill the
        % window each time
        axis vis3d
    end
    camzoom(zoom_factors(i))
    drawnow
    fn = [side '_' num2str(frame_index) '_' ViewAngles{i}];
    frames_set(i) = getframe;
    
    % before you save this picture, add a colorbar
    if frame_index == 1
        colorbar
        print(gcf, fn, '-dtiff');
    end
    
    camzoom(1/zoom_factors(i))
    
    % you have to turn the light off before applying the new light for
    % the next position. The light will be the most recent axes child
    axes_children = get(gca,'children');
    delete(axes_children(1))
    
    % also take off the colorbar so that the next image will not
    % have a colorbar until right before you save it
    if frame_index == 1
        delete(colorbar)
    end
end
end
