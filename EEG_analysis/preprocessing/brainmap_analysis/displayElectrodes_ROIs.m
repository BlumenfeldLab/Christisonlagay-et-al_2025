% Inputs: overlay_frame_data: 0 = nothing, just draw electrodes on blank brain
%                             1 = arrival_times
%                             2 = zScores confirmed perceived
%                             3 = zscores confirmed not perceived
%         createMontage: 1 = create each patient's montage and map for EEG channels, 0 = nothing
%         laterality: collect data and draw on: 0 = bilateral, 1 = left only, 2 = right only, 3 = project both left and right electrodes on the left brain
%         views: 4 = 4 views (lateral, medial, ventral, and posterior)
%                2 = views (Left medial and right lateral)

data_location='D:\Human CPT icEEG Project\icEEG Analysis\Group Analysis\CPT ERP movies';
patients ={'193AF','154NP','200JW','201MU','210VG','212DD','213SP','215JF','217CG','168KG','278JA'};
addpath('D:\Human CPT icEEG Project\icEEG Analysis\Group Analysis\ROI Analysis')


ROI_Path='D:\Human CPT icEEG Project\icEEG Analysis\Group Analysis\ROI Analysis\FOrmatted ROI\Formatted ROIs2\MNI_';

mni_nifti_path = 'custom_medial_frontal.nii.gz';
mni_nifti=gunzip(mni_nifti_path);
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti{1});
medialFrontal=mnit1_1mm_stripped.img;

ROIs={'Calcarine','Cuneus','Lingual','Occipital_Sup','Occipital_Mid','Occipital_Inf','Fusiform','Frontal_Mid','Frontal_Sup_Medial','Parahippocampal','Parietal_Sup','Angular','Parietal_Inf','SupraMarginal','Supp_Motor_Area','Precentral'};

ROIs_titles={'Occipital','Fusiform','Caudal Middle Frontal','Superior Frontal','Parahippocampal','Parietal','Supplementary Motor Area','Precentral'};

ROIcmap =[
    227 71 77 %red
    228 157 70 %orange
    236 88 247 %magenta
    44 155 100%green;
    166 117 178 %purple
    250 217 120 %gold
    72 97 166  %blue
    70 163 201 %cyan
    ]/255;

ROIcmap= [0.905882358551025	0.905882358551025	0.905882358551025;ROIcmap;0.905882358551025	0.905882358551025	0.905882358551025];

overlay = 2;
createMontage = 0;
laterality = 0;
views = 4;

image_storage_folder = [data_location '\Wendy picture storage2'];
xtick_spacing_filename = 'spectrogram_times_64bins_2s.mat';

% list of patients with only right sided contacts

VTK_folder_name = 'VTK Leah';

for ROI_counter=1:length(ROIs_titles)
    plot(10,10,'s','LineWidth',1.5,'MarkerSize',16,'MarkerEdgeColor','k','MarkerFaceColor',ROIcmap(ROI_counter+1,:)),hold on
end

lgd=legend(ROIs_titles,'FontSize',14);
title(lgd,'Regions of Interest','FontSize',16)
legend('boxoff')
metadata = CollectEEGMetadata(data_location,VTK_folder_name);

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
ROI_vertices_fusiform=[];

%% Generating all frames for each side of the brain
for side_index = start_side:sides_num    % Generate frames for L and R
    if(side_index == 1)
        side = 'L';
    else
        side = 'R';
    end
    %------------------------------------------------------------------
    %   Set-up some local variables
    %------------------------------------------------------------------
    cortex_filename{side_index}  =    [side '_cortex.txt'];
    smooth_filename{side_index}  =    [side '_smooth.txt'];
    cd([data_location '\' VTK_folder_name]);
    cortex_v_num(side_index)     =    eval(['metadata.' side '_cortex_v_num']);
    cortex_f_num(side_index)     =    eval(['metadata.' side '_cortex_f_num']);
    smooth_v_num(side_index)     =    eval(['metadata.' side '_smooth_v_num']);
    smooth_f_num(side_index)    =    eval(['metadata.' side '_smooth_f_num']);
    
    % Input and interpret the vertice data for the cortex and smooth
    % surfaces
    [cortex(side_index).vertices, cortex(side_index).faces] = FormatMeshSurface(cortex_v_num(side_index), cortex_f_num(side_index), data_location, cortex_filename{side_index});
    [smooth(side_index).vertices, smooth(side_index).faces] = FormatMeshSurface(smooth_v_num(side_index), smooth_f_num(side_index), data_location, smooth_filename{side_index});
    
    
    %% create a figure of the left or right side of the brain
    % you will be working with this for different electrodes, overlays,
    % views, and light
    
    figHandle(side_index) = figure('Position',[70,70,700,700]);
    g=trisurf(cortex(side_index).faces,cortex(side_index).vertices(:,1),cortex(side_index).vertices(:,2),cortex(side_index).vertices(:,3));
    
    lighting flat;
    set(g, 'FaceColor', [1 1 1], 'LineStyle', 'none')
    set(gcf, 'color', [1 1 1]);
    hold on
    axis off;
    axis equal;
    
    %% hard coded number of frames for each overlay setting, need to change later
    if overlay == 1 || overlay == 0
        total_num_frames = 1;
    elseif overlay == 2 || overlay == 3 || overlay == 4 || overlay == 5 || overlay == 6%
        total_num_frames = 64;
    end
    
    
    vertex_values = nan(cortex_v_num(side_index),1);
    counter=0;
    for ROI_counter=[1 7 8  9:11 15 16]
        counter= counter+1;
        if(side_index == 1)
            addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter},'_L_roi')))
        else
            addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter},'_R_roi')))
        end
        vertexCdata = [];
        if overlay == 1
            electrode_mapping_times = 10;
        else
            electrode_mapping_times = 1;
        end
        
        for e = 1:electrode_mapping_times
            
            scatter_count = 0;
            ROI_electrode_values=[];
            load('voxelInfo.mat')
            voxel_cor=xyz.vXYZ;
            
            switch ROI_counter
                case 1
                    
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+1},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+1},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor2=xyz.vXYZ;
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+2},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+2},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor3=xyz.vXYZ;
                    voxel_cor3 = (union(voxel_cor3',(ROI_voxel_projection_Z(cortex(side_index).vertices,voxel_cor3)), 'rows'))';
                    voxel_cor =  (union(voxel_cor3',union(voxel_cor', voxel_cor2', 'rows'), 'rows'))';
                    
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat('Fusiform','_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat('Fusiform','_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor_fusiform=xyz.vXYZ;
                    voxel_cor_fusiform = ROI_voxel_projection_Z(cortex(side_index).vertices,voxel_cor_fusiform);
                    voxel_cor_MedialOccipital=(setdiff(voxel_cor',voxel_cor_fusiform,'rows'))';
                    
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+3},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+3},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor=xyz.vXYZ;
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+4},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+4},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor2=xyz.vXYZ;
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+5},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+5},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor3=xyz.vXYZ;
                    
                    voxel_cor_LateralOccipital =  (union(voxel_cor3',union(voxel_cor', voxel_cor2', 'rows'), 'rows'))';
                    voxel_cor=  union(voxel_cor_MedialOccipital', voxel_cor_LateralOccipital', 'rows')';
                case 8
                    [idx,center]=kmeans(voxel_cor',2);
                    voxel_cor=voxel_cor(:,find(idx==(find(center(:,2)==max(center(:,2))))));
                    
                case 11
                    
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+1},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+1},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor2=xyz.vXYZ;
                    voxel_cor =(union(voxel_cor', voxel_cor2', 'rows'))';
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+2},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+2},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor3=xyz.vXYZ;
                    voxel_cor =(union(voxel_cor', voxel_cor3', 'rows'))';
                    
                    if(side_index == 1)
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+3},'_L_roi')))
                    else
                        addpath(strcat(ROI_Path,strcat(ROIs{ROI_counter+3},'_R_roi')))
                    end
                    load('voxelInfo.mat')
                    voxel_cor4=xyz.vXYZ;
                    voxel_cor =(union(voxel_cor', voxel_cor4', 'rows'))';
            end
            
            if(strcmp(ROIs{ROI_counter},'Fusiform')) %  For   Fusiform
                voxel_cor = ROI_voxel_projection_Z(cortex(side_index).vertices,voxel_cor);
                ROI_vertices = ROI_mapping_voxel2vetices(cortex(side_index).vertices,voxel_cor');
            elseif(strcmp(ROIs{ROI_counter},'Frontal_Sup_Medial')) %  For   Fusiform
                ROI_vertices = ROI_mapping_voxel2vetices3(cortex(side_index).vertices,voxel_cor,side_index);
            elseif(strcmp(ROIs{ROI_counter},'Frontal_Mid')) %  For   Fusiform
                ROI_vertices = ROI_mapping_voxel2vetices4(cortex(side_index).vertices,voxel_cor,side_index);
            elseif(strcmp(ROIs{ROI_counter},'Frontal_Sup_Orb'))
                voxel_cor = ROI_voxel_projection_Z2(cortex(side_index).vertices,voxel_cor);
                ROI_vertices = ROI_mapping_voxel2vetices(cortex(side_index).vertices,voxel_cor');
                ROI_vertices=ROI_vertices(find((cortex(side_index).vertices(ROI_vertices,2))<=maxy));
                ROI_vertices=setdiff(ROI_vertices,ROI_vertices_tmp);
            elseif (strcmp(ROIs{ROI_counter},'Frontal_Med_Orb'))
                ROI_vertices_custom = ROI_mapping_voxel2vetices2(cortex(side_index).vertices,voxel_cor_custom);
                voxel_cor_Z = ROI_voxel_projection_Z2(cortex(side_index).vertices,voxel_cor);
                voxel_cor_Y = ROI_voxel_projection_Y(cortex(side_index).vertices,voxel_cor');
                ROI_vertices=union(ROI_mapping_voxel2vetices(cortex(side_index).vertices,voxel_cor_Z'),ROI_mapping_voxel2vetices(cortex(side_index).vertices,voxel_cor));
                maxy=max(cortex(side_index).vertices(ROI_vertices,2));
                maxx=max(cortex(side_index).vertices(ROI_vertices,1));
                minx=min(cortex(side_index).vertices(ROI_vertices,1));
                ROI_vertices=union(ROI_vertices,ROI_vertices_custom);
                if(side_index == 1)
                    ROI_vertices=ROI_vertices(find((cortex(side_index).vertices(ROI_vertices,1))<=(maxx-6)));
                else
                    ROI_vertices=ROI_vertices(find((cortex(side_index).vertices(ROI_vertices,1))>=(minx+6)));
                end
                ROI_vertices_tmp=ROI_vertices;
            else
                ROI_vertices = ROI_mapping_voxel2vetices(cortex(side_index).vertices,voxel_cor);
                ROI_vertices=setdiff(ROI_vertices,ROI_vertices_fusiform);
            end
            
            for p = 1
                
                patient_folder = [data_location, '/', patients{p}];
                
                %% Make the Montage
                % if the left/right electrode montages have not been created, do so
                % now, indicated by createMontage = 0 or = 1
                if createMontage == 1
                    if  laterality == 3
                        flipX = 1;
                        create_electrode_montage(patient_folder,flipX)
                    else
                        flipX = 0;
                        create_electrode_montage(patient_folder,flipX)
                    end
                end
                display(patients{p})
                
                % Now load some data files needed for electrode placement and surface shading
                cd([data_location, '/', patients{p}]);
                load([side '_MontageMap.mat']);
                electrode = eval([side '_MontageMap' ]);
                
                % load the X-flipped montage and add to list of electrodes
                if laterality == 3
                    load('R_MontageMap_flipX.mat');
                    electrode = cat(1,electrode,R_MontageMap);
                end
                
                %% Discard depth electrodes or plot only some electrodes for arrival times
                if ~isempty(electrode)
                    % find out which depth electrodes should be discarded from the
                    % analysis
                    load('labels.mat')
                    load('labels_depth_most.mat')
                    
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
                    if overlay == 1
                        load('arrival_times_by_period_100ms.mat')
                        electrodes_this_period = find(arrival_times_by_period(:,e) == 1);
                        electrode = electrode(ismember(electrode(:,1),electrodes_this_period),:);
                    end
                end
                
                %% Plot patient's electrodes
                if isempty(electrode)
                    disp('-- No electrode data!');
                    electrodes_present(side_index,p) = 0;
                else
                    all_electrodes = electrode(:,1);
                    electrodes_present(side_index,p) = 1;
                    
                    [elecVertices2]= ProjectElectrode2TransSurf(smooth_v_num, smooth(side_index).vertices, all_electrodes, electrode, side_index);
                    electrode(:,2:4)                       = smooth(side_index).vertices(elecVertices2,1:3);
                    [elecVertices] = ProjectElectrode2TransSurf(cortex_v_num, cortex(side_index).vertices, all_electrodes, electrode, side_index);
                    scatter_count = scatter_count + 1;
                end
            end
            [Lia, Locb] = ismember((1:length(cortex(side_index).vertices))',ROI_vertices','rows');
            vertex_values(find(Lia))=counter;
            ROI_vertices=[];
        end
    end
    
    %% Frame Generation if saving without overlaying color onto the brain
    if overlay == 1
        display(e)
        figure(figHandle(side_index))
        cd(image_storage_folder)
        frames_set(side_index,e,:) = savebrainimages(side,e);
        children = get(gca,'children');
        delete(children(1:scatter_count,1))
    end
    
    if ~ismember(overlay,1)

        load ([image_storage_folder '\cmap6to8WRX.mat']);
        
        for frame_index = 1
            
            %% overlay the colors
            
            disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status
            
            % sum the count of vertices across patients
            vertex_values_sum = nan(cortex_v_num(side_index),1);
            if overlay == 2 || overlay == 3 || overlay == 4 || overlay == 5 || overlay == 6 % the surface is colored according to zscores of power
                vertex_values_sum=vertex_values;
                % load the colormap and set axis of colormap
                colormap(ROIcmap)
                caxis([0 length(ROIs_titles)+1]);
                %caxis([-15 20]);
                if any(electrodes_present(side_index,:) == 1)
                    
                    h = patch('Faces',cortex(side_index).faces, 'Vertices', cortex(side_index).vertices, 'FaceVertexCData',vertex_values_sum,'FaceColor','interp',...
                        'FaceLighting','gouraud','SpecularStrength',0.2,'LineStyle', 'none',  'EdgeLighting', 'gouraud');
                    alpha = 0.7;
                    alpha_data = ones(numel(h.FaceVertexCData), 1) * alpha;
                    alpha_data(abs(h.FaceVertexCData) <= 2) = 0;
                    %                         lightset = [0.6 0.5 0.1];
                    lightset = [0.5 0.5 0.1];
                    material(lightset);
                    set(h, 'FaceVertexAlphaData', alpha_data)
                end
            end
            
            %% Save each side of the brain in 3 views
            frames_set(side_index,frame_index,:) = savebrainimages(side_index,frame_index,views);
            if any(electrodes_present(side_index,:) == 1) && exist('h','var')
                delete(h)
            end
        end
    end
    
    if side_index == 1
        vertex_valuesL = vertex_values;
    elseif side_index == 2
        vertex_valuesR = vertex_values;
    end
end
save matlab.mat;
if laterality == 0
    save('vertex_values_allevents_64bins(filter).mat','vertex_valuesL','vertex_valuesR','-v7.3')
    save('vertex_values_allevents_64bins.mat','vertex_valuesL','vertex_valuesR')
    cd(data_location);
    save('vertex_values_allevents_64bins(filter).mat','vertex_valuesL','vertex_valuesR','-v7.3')
    save('vertex_values_allevents_64bins.mat','vertex_valuesL','vertex_valuesR')
elseif laterality == 3 || laterality == 1
    save('vertex_values.mat','vertex_valuesL','-v7.3')
elseif laterality == 2
    save('vertex_values.mat','vertex_valuesR','-v7.3')
end

% Put the 3 views into one image
if overlay == 1
    all_frames_num = electrode_mapping_times;
else
    all_frames_num = total_num_frames;
end

for frame_index = 1 : 1
    
    %% name the frames for each side
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
        end
        for view_index=1:length(ViewAngles)
            Frame = frames_set(side_index,frame_index,view_index).cdata;
            eval([side '_' ViewAngles{view_index} ' = Frame;'])
        end
    end
    
    %% creat the frames for each side
    if views == 4
        for side_index = start_side:sides_num
            if side_index == 1
                L_Lateral_Crop = L_Lateral;
                %             imshow(L_Lateral_Crop)
                
                L_Medial_Crop = [L_Medial;ones(5,544,3)*255];
                %             figure
                %             imshow(L_Medial_Crop)
                
                %L_Ventral_Crop = L_Ventral(:,1:225,:);
                %L_Ventral_Crop = [L_Ventral;ones(17,230,3)*255];
                
                
                %Aya
                L_Ventral_Crop = [L_Ventral;ones(17,231,3)*255];
                
                
                
                %             figure
                %             imshow(L_Ventral_Crop)
                
                %             L_Ventral_padded = [ones(15,1090,3)*255; L_Ventral_Crop,ones(572,839,3)*255];
                %             L_Lateral_Medial = [L_Medial_Crop, L_Lateral_Crop];
                %             combined_views_left{frame_index} = [L_Lateral_Medial; L_Ventral_padded];
                L_Lateral_Medial = [L_Medial_Crop; L_Lateral_Crop];
                %                 L_Lateral_Medial_padded = [L_Lateral_Medial; ones(177,544,3)*255];
                
                %Aya
                L_Lateral_Medial_padded = [L_Lateral_Medial; ones(175,544,3)*255];
                
                
                %             figure
                %             imshow(L_Lateral_Medial)
                
                if laterality == 0
                    R_Posterior_Crop = [R_Posterior, ones(413,5,3)*255];
                    %               figure
                    %               imshow(R_Posterior_Crop)
                    
                    %                     L_Ventral_R_Posterior = [L_Ventral_Crop;R_Posterior_Crop];
                    
                    
                    %Aya
                    L_Ventral_R_Posterior = [L_Ventral_Crop(:,1:230,:);R_Posterior_Crop];
                    
                    
                    %             figure
                    %             imshow(L_Ventral_R_Posterior)
                    L_Lateral_Medial_padded_L_Ventral_R_Posterior = [L_Ventral_R_Posterior, L_Lateral_Medial_padded];
                    combined_views_left{frame_index} = [ones(974,11,3)*255,L_Lateral_Medial_padded_L_Ventral_R_Posterior];
                elseif laterality == 3
                    %                     L_Posterior_Crop = [L_Posterior(1:396,1:225,:)];
                    L_Posterior_Crop = L_Posterior;
                    %                     figure
                    %                     imshow(L_Posterior_Crop)
                    
                    L_Ventral_L_Posterior = [L_Ventral_Crop; L_Posterior_Crop];
                    %                     figure
                    %                     imshow(L_Ventral_L_Posterior)
                    
                    L_Lateral_Medial_padded= L_Lateral_Medial_padded(1:958,:,:);
                    combined_views_left{frame_index} = [L_Ventral_L_Posterior, L_Lateral_Medial_padded];
                end
                
                figure
                imshow(combined_views_left{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                fn = ['combined_viewsL' '_' num2str(frame_index)];
                saveas(gcf,['D:\Human CPT icEEG Project\Analysis Results\v1.0' '\' fn '.tiff'])
                
                %             export_fig(fn, '-tiff');
                print(fn, '-dtiff');
                close
            elseif side_index == 2
                
                %                 R_Lateral_Crop = [ones(5,544,3)*255;R_Lateral];
                
                %                 Aya
                R_Lateral_Crop = [ones(5,545,3)*255;R_Lateral];
                
                %             figure
                %             imshow(R_Lateral_Crop)
                
                R_Medial_Crop = R_Medial(1:413,:,:);
                %             figure
                %             imshow(R_Medial_Crop)
                R_Lateral_Medial = [R_Medial_Crop; R_Lateral_Crop];
                %                 R_Lateral_Medial_padded = [R_Lateral_Medial; ones(126,544,3)*255];
                %                 Aya
                R_Lateral_Medial_padded = [R_Lateral_Medial; ones(126,545,3)*255];
                
                R_Lateral_Medial_padded = [R_Lateral_Medial_padded,ones(957,11,3)*255];
                %             figure
                %             imshow(L_Lateral_Medial_padded)
                
                %                 R_Ventral_Crop = [R_Ventral, ones(544,5,3)*255];
                
                
                R_Ventral_Crop = [R_Ventral, ones(545,5,3)*255];
                
                %             figure
                %             imshow(R_Ventral_Crop)
                
                %             R_Ventral_padded = [ones(572,849,3)*255, R_Ventral_Crop];
                %             R_Lateral_Medial = [R_Lateral_Crop, R_Medial_Crop];
                %             combined_views_right{frame_index} = [R_Lateral_Medial; R_Ventral_padded];
                
                if laterality == 0
                    %                     L_Posterior_Crop = [ones(17,230,3)*255; L_Posterior];
                    %Aya
                    L_Posterior_Crop = [ones(15,230,3)*255; L_Posterior(:,1:230,:)];
                    %               figure
                    %               imshow(L_Posterior_Crop)
                    
                    R_Ventral_L_Posterior = [R_Ventral_Crop;L_Posterior_Crop];
                    %               figure
                    %               imshow(R_Ventral_R_Posterior)
                    R_Lateral_Medial_padded_R_Ventral_L_Posterior = [R_Lateral_Medial_padded, R_Ventral_L_Posterior];
                    
                    %                     combined_views_right{frame_index} = [R_Lateral_Medial_padded_R_Ventral_L_Posterior;ones(17,785,3)*255];
                    combined_views_right{frame_index} = [R_Lateral_Medial_padded_R_Ventral_L_Posterior;ones(17,786,3)*255];
                    
                elseif laterality == 2
                    R_Posterior_Crop = R_Posterior(10:405,1:225,:);
                    %               figure
                    %               imshow(R_Posterior_Crop)
                    
                    R_Ventral_R_Posterior = [R_Ventral_Crop;R_Posterior_Crop];
                    %               figure
                    %               imshow(R_Ventral_R_Posterior)
                    
                    combined_views_right{frame_index} = [R_Lateral_Medial_padded, R_Ventral_R_Posterior];
                end
                
                figure
                imshow(combined_views_right{frame_index});
                set(gca,'position',[0 0 1 1],'units','normalized')
                fn = ['combined_viewsR' '_' num2str(frame_index)];
                %             export_fig(fn, '-tiff');
                saveas(gcf,[fn '.tiff'])
                %             print(fn, '-dtiff');
                close
            end
        end
        txt_pos_x = 1100;
        txt_pos_y = 900;
    elseif views == 2
        txt_pos_x = 200;
        txt_pos_y = 200;
        L_Medial_Crop = L_Medial(1:396,:,:);
        R_Lateral_Crop = R_Lateral(10:405,:,:);
    end
    
    %------------------------------------------------------------------
    %   Movie Assembly Procedure
    %------------------------------------------------------------------
    % Add all 6 components of each movie frame together onto one plot after
    % cropping. Then after assembling the frame, add to an avi movie buffer
    % output
    
    if laterality == 0
        figure
        if views == 4
            combined_views = [combined_views_right{frame_index} combined_views_left{frame_index}];
        elseif views == 2
            combined_views = [L_Medial_Crop R_Lateral_Crop];
        end
        imshow( combined_views );
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis tight
        frame(frame_index,:,:,:) = combined_views;
        if exist('j','var')
            fn = ['combined_views_full_' num2str(frame_index) '_' num2str(j)];
        else
            fn = ['combined_views_full_' num2str(frame_index)];
        end
        
        if ismember(overlay,[2 3 4 5 6])
            load(xtick_spacing_filename)
            time_txt = ['Time = ' num2str((T(frame_index))) 's'];
        elseif overlay == 1
            T = [100 200 300 400 500 600 700 800 900 1000];
            time_txt = ['Time = ' num2str(T(frame_index)) 's'];
        end
        %         text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
        colormap(ROIcmap)
        caxis([1 length(ROIs_titles)]);
        %caxis([-15 15]);
        %         colorbar
        %         saveas(gcf,[fn '.tiff'])
        
        
        
        %         saveas(gcf,[data_location '\' patients{1} '\' 'events' '\' fn '.tiff'])
        
        saveas(gcf,['D:\Human CPT icEEG Project\Analysis Results\v30.0' '\' fn '.tiff'])
        
        
        %         export_fig(fn, '-tiff');
        %         print(fn, '-dtiff');
        % clear combined_views_with_colorbar; clear combined_views_left; clear combined_colorbars;
        % clear Beta_L_Lateral_Crop; clear Beta_L_Medial_Crop; clear Beta_L_Ventral_Crop;
        % clear Beta_R_Lateral_Crop; clear Beta_R_Medial_Crop; clear Beta__Ventral_Crop;
        % clear combined_views;
        
        close
    end
end
% if you want to create a montage of evenly space post-stimulus frames,
% then uncomment this code

if views == 2
    %indices of frames to capture in montage
    %frame_idx = [33 37 41 45 49 53 57];
    %     frame_idx = [57 53 49 45 41 37 33];
    frame_idx = [60 56 52 48 44 40 36];large_image = [];
    for i = 1:length(frame_idx)
        large_image = cat(1,squeeze(frame(frame_idx(i),:,:,:)),large_image);
    end
    figure
    imshow(large_image)
    set(gca,'position',[0 0 1 1],'units','normalized')
    saveas(gcf,'frames_montage.tiff')
    saveas(gcf,'frames_montage.eps')
end

% end
close all

%% Below are functions called upon within this function


%==================================================================
%
%Title:         savebrainimages
%Description:   This function takes a brain you have prepared and shines a
%               light on it in 3 views and saves them
%
%Inputs:        "side"          [str] left = 'L', right = 'R'
%               "frame_index"   [int] the "frame" of the "movie" you are
%                               on, which can be #clusters, time period, or
%                               bin number for power analysis
%==================================================================

function frames_set = savebrainimages(side_index,frame_index,views)

if views == 4
    % For each frame, take a snapshot for each of these views
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    % These arrays define the viewing and light perspective of the frame
    if side_index == 1
        side = 'L';
        Light_Pos = [1200 -200 500; -700 -800 -100;-200, -100, -800; 200, 800, 0];
        View_Pos = [90 0;-90, 0;0 -90;180 0];
        %         View_Pos_xyz = [200 0 0;-200 0 0;0 0 -200;0 1000 0];
    elseif side_index == 2
        side = 'R';
        Light_Pos = [-1200 -300 700; 1000 -1100 -200; -200, -100, -800;-200, 800, 0];
        View_Pos = [-90 0;90, 0;0 -90;180 0];
    else
        Light_Pos=0;
        View_Pos=0;
    end
elseif views == 2
    if side_index == 1
        side = 'L';
        ViewAngles = {'Medial'};
        Light_Pos = [-700 -800 -300];
        View_Pos = [-90 -30];
    elseif side_index == 2
        side = 'R';
        ViewAngles = {'Lateral'};
        Light_Pos = [-1200 -300 -200];
        View_Pos = [-90 -15];
    else
        Light_Pos=0;
        View_Pos=0;
    end
end

for i = 1:length(ViewAngles)
    lightHandle = light('Position', Light_Pos(i,:));
    
    view(View_Pos(i,:));
    if i == 1
        % with this command, matlab will freeze the aspect ratio of
        % your image instead of changing the object size to fill the
        % window each time
        axis vis3d
    end
    drawnow
    %set(gcf, 'Position', [0 0 545 572]);
    fn = [side '_' num2str(frame_index) '_' ViewAngles{i}];
    frames_set(i) = getframe;
    %         frames_set(side_index,frame_index,i) = getframe;
    
    % before you save this picture, add a colorbar
    %if frame_index == 1
    %             colorbar
    print(gcf, fn, '-dtiff');
    %             saveas(gcf, fn);
    %end
    
    % you have to turn the light off before applying the new light for
    % the next position. The light will be the most recent axes child
    axes_children = get(gca,'children');
    delete(axes_children(1))
    
    % also take off the colorbar so that the next image will not
    % have a colorbar until right before you save it
    %if frame_index == 1
    delete(colorbar)
    %end
end
end

%==================================================================
%
%Title:         FormatMeshSurface
%Description:   This function loads a surface txt file output from Bioimage
%               Suite into Matlab, and arranges it into a readable format.
%               Repeat twice (once for cortex and once for smooth). The
%               output is a reorganized matrix of vertices.
%
%Inputs:        "vertices_num"  [int] The number of vertices in the file
%               "faces_num"     [int] The number of faces in the file
%               "data_location" [string] should point to the root folder
%                               for the seizure
%                               i.e. 'C:\IcEEGProject_ProcessedData\25\25_1'
%               "filename"      [string] Location of the vertice data
%==================================================================

function [mesh_vertices, mesh_faces]=FormatMeshSurface(vertices_num, faces_num, data_location, filename)
% cd(patient_folder);
% cd ..
cd(data_location)
data_location_previous=pwd;
vertices_range_str=strcat('A1..I',num2str(ceil(vertices_num/3)));
raw_mesh_vertices=dlmread(fullfile( data_location_previous, filename),'',vertices_range_str);
faces_range_str=strcat('B',num2str(ceil(vertices_num/3)+1),'..D',num2str(ceil(vertices_num/3)+faces_num));
raw_mesh_faces=dlmread(fullfile( data_location_previous, filename),'',faces_range_str );
mesh_faces=raw_mesh_faces+1;

vertices_range=1:ceil(vertices_num/3);
mesh_vertices=zeros(ceil(vertices_num/3)*3,3);
mesh_vertices(3*vertices_range-2,1:3)=raw_mesh_vertices(vertices_range, 1:3);
mesh_vertices(3*vertices_range-1,1:3)=raw_mesh_vertices(vertices_range, 4:6);
mesh_vertices(3*vertices_range,1:3)=raw_mesh_vertices(vertices_range, 7:9);
mesh_vertices(ismember(mesh_vertices,zeros(1,3),'rows'),:)=[];
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
function [elecVertices]=ProjectElectrode2TransSurf(vertices_num, vertices, all_electrodes, electrode, side_index)
elecVertices = zeros(size(electrode,1),1);
for elec_index = 1:size(electrode,1)
    xCoord = ones(vertices_num(side_index),1)*electrode(elec_index,2);
    yCoord = ones(vertices_num(side_index),1)*electrode(elec_index,3);
    zCoord = ones(vertices_num(side_index),1)*electrode(elec_index,4);
    [minVal, minInd] = min(sqrt((vertices(:,1) - xCoord).^2 + (vertices(:,2) - yCoord).^2 + (vertices(:,3) - zCoord).^2));
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
load('labels_reduced.mat');
labels = labels';
b1 = linspace(1, length(labels), length(labels))';
b1 = num2cell(b1);
b2 = labels(:,1);
montage = [b1 b2];
[a1, a2, map]=xlsread(['Map_reduced.xlsx']);
% Convert the long channel names (i.e "A_L_Most_Mid_Frontal_Polar_12" into the abbreviated names found in the montage (i.e. "A12")
try
    for i=1:size(map,1)
        str             =   map{i,1};
        str_delim       =   regexp(str,'_','split');
        if length(str_delim)==1 || length(str_delim) == 2
            map{i,5}=strtrim(str_delim{1,1}(end));
        else
            firstLetter     =   str_delim{1};
            positionLetter  =   strtrim(str_delim{2});
            digitStr        =   str(isstrprop(str,'digit'));
            replaceStr      =   strcat(positionLetter, digitStr);
            map{i,1}        =   replaceStr; % i.e. "A12"
            map{i,5}        =  'L'; % Either L or R
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
    if (Position(i)=='L');
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

%==================================================================
%
%Title:         electrode_data_overlay
%Description:   This function loads a surface txt file output from Bioimage
%               Suite into Matlab, and arranges it into a readable format.
%               Repeat twice (once for cortex and once for smooth). The
%               output is a reorganized matrix of vertices.
%
%Inputs:        "vertices_num"  [int] The number of vertices in the file
%               "faces_num"     [int] The number of faces in the file
%               "data_location" [string] should point to the root folder
%                               for the seizure
%                               i.e. 'C:\IcEEGProject_ProcessedData\25\25_1'
%               "filename"      [string] Location of the vertice data
%==================================================================

function [electrode_vertex_values] = electrode_data_overlay(overlay,frame_index,all_electrodes)
overlay_frame_data = [];
if overlay == 0 || overlay == 1
    %% Give every overlay value a 1
    % once aggregated across patients it gives the density map of
    % electrodes at each vertex
    load('labels.mat');
    overlay_frame_data = ones(length(labels),1);
    % elseif overlay == 1
    %     %% First arrival times - need to add code for 2nd and 3rd------------------
    %     load('arrival_time.mat')%, 'half_max_times');
    %     newTimes = [];
    %     for i = 1:size(half_max_times,2)
    %         arrivalOne = half_max_times{1,i};
    %         if ~isempty(half_max_times{1,i})
    %             firstTime = arrivalOne(1);
    %         else
    %             firstTime = 0;
    %         end
    %         newTimes = [newTimes, firstTime];
    %     end
    %     newTimes = newTimes';
    %     overlay_frame_data = newTimes;
elseif overlay == 2 || overlay == 3 || overlay == 4 || overlay == 5 || overlay == 6
    %% getting power zscores from mean power of each electrode
    %     load('meanpower_traces_zscore_event_4s_64bins.mat')
    %     load('meanpower_traces_zscore_allevents_64bins.mat')
    
    
    %      load('meanpower_traces_events.mat')
    load('meanpower_traces_events_Zscores_Corrected.mat')
    %     load('meanpower_traces_events_Zscores.mat')
    if overlay == 2
        trialtype = 1; % block onset transient
    elseif overlay == 3
        trialtype = 2; % block offset transient
    elseif overlay == 4
        trialtype = 3; % button
    elseif overlay == 5 % sum
        trialtype = 4;
    elseif overlay == 6 % Substraction
        trialtype = 5;
        %         meanpower_traces(:,5,:) = meanpower_traces(:,1,:) - meanpower_traces(:,3,:);
        meanpower_traces(:,5,:) = (meanpower_traces(:,1,:) + meanpower_traces(:,3,:))/2;
        
    end
    overlay_frame_data = squeeze(meanpower_traces(:,trialtype,:));
end

% once you have the overlay data, you need to match it to an electrode from
% labels

overlay_frame_data = overlay_frame_data';
electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);

end


%% Determine the transparency for every face/vertex


%==================================================================
%
%Title:         CalculateTransparency
%Description:   This function takes the input values from overlay frame
%               data for each electrode and gives it to all the vertices
%               within radius2 distance from each electrode's vertex
%
%Inputs:        "vertices_num"  The number of vertices in the file
%               "elecVertices"  vertex indices of each electrode to be
%                               plotted on this side of the brain
%               "electrode_vertex_values"
%                               the overlay value assigned to each of those
%                               vertices in elecVertices
%               "overlay"      	[int] this number is an input to the parent
%                               function. If overlay = 0 then it does not
%                               do any weighing based on other electrodes
%                               within a radius2 proximity. All vertices
%                               within radius 2 receive a 1
%
%Outputs:       "vertexCdata"   matrix the length of the number of vertices
%                               with the electrode_vertex_values assigned
%                               to all vertices within radius2 of each
%                               electrode
%==================================================================

%This function should be combined with OverlayColor to speed up processing
%time...

function [ vertexCdata]= CalculateTransparency(vertices_num,vertices, elecVertices, electrode_vertex_values, overlay)
maxAlpha    = 1;
radius1     = 1;
radius2     = 15;
electrodes_num      = length(elecVertices);
electrode_vector    = vertices(elecVertices,:);
vertexCdata = zeros( vertices_num, 1 );
j=1;
k=1;

for vertIndex = 1 : vertices_num
    sharedTransparency = [];
    sharedElectrode = [];
    
    xCoord = ones(electrodes_num,1)*vertices(vertIndex, 1);
    yCoord = ones(electrodes_num,1)*vertices(vertIndex, 2);
    zCoord = ones(electrodes_num,1)*vertices(vertIndex, 3);
    
    % Calculate "Transparency" in other words the "weight" based on
    % distance from the vertex
    [distanceToElectrode, electrodeIndice] = sort(sqrt((electrode_vector(:,1) - xCoord).^2 + (electrode_vector(:,2) - yCoord).^2 + (electrode_vector(:,3) - zCoord).^2));
    
    if overlay == 0
        if any(distanceToElectrode <= radius2)
            vertexCdata(vertIndex,1) = 1;
        end
    else
        
        for n = 1:length(electrodeIndice)
            if distanceToElectrode(n) <= radius1
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                sharedTransparency = [sharedTransparency; maxAlpha];
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
            elseif distanceToElectrode(n) <= radius2
                % Calculates the transparencies for each electrode within a
                % distance less than 15mm from this vertex
                sharedTransparency = [sharedTransparency; (maxAlpha - ((distanceToElectrode(n) - radius1)/(radius2-radius1))*maxAlpha)];
                % Saves the indices/index of the electrodes/electrode that contribute transparency to this vertex
                sharedElectrode = [sharedElectrode; electrodeIndice(n)];
            end
        end
        % Grabs the zScore values associated with the electrodes that
        % contribute power to this vertex
        zIn = electrode_vertex_values(sharedElectrode)';
        weightedZ = []; %transparencySum = [];
        % Calculates the weight of the z score according to the transparencies and zScores associated with this vertex
        for h = 1:length(zIn)
            %                 weightedZ = [weightedZ zIn(h)*(sharedTransparency(h).^2)]
            weightedZ = [weightedZ zIn(h)*(sharedTransparency(h))];
            %                 transparencySum = [transparencySum sharedTransparency(h)];
        end
        
        
        %             indDelete = [];
        %             for y = 1:length(zIn)
        %                 if isnan(weightedZ(y))
        %                     indDelete = [indDelete y];
        %                 end
        %             end
        %             weightedZ(indDelete) = [];
        %             transparencySum(indDelete) = [];
        
        
        weightedZ = nansum(weightedZ);
        %             transparencySum = sum(transparencySum);
        %             weightedZ = weightedZ/transparencySum;
        % If there is no weighted zScore that means the vertex was greater than
        % 15mm from eevery electrode and has no power or zScore to display in
        % this frame
        if isempty(zIn)
            vertexCdata(vertIndex,1) = 0;
        elseif isnan(zIn)
            vertexCdata(vertIndex,1) = 0;
        else
            vertexCdata(vertIndex,1) = weightedZ;
        end
    end
end
end
