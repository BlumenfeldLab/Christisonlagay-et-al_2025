mappeddatadrive='V:';

mni2fs_dir='D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\mni2fs-master';
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);



mni_nifti_path = fullfile(mni2fs_dir, '\MNI_T1_1mm_stripped.nii');

mnit1_1mm_stripped = niftiread(mni_nifti_path);

load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

inflationstep=5;




% VTK_folder_name = 'VTK Leah';  %%surface
% data_location='D:\Human CPT icEEG Project\icEEG Analysis\Group Analysis\CPT ERP movies';  %%where data lives
%addpath('D:\Human CPT icEEG Project\icEEG Analysis\Analysis Codes\nifti_tools\nifti_tools')  %%%read MRI, ignore here
%file= gunzip('D:\ch2better_WhiteMatter2.nii.gz');  %%this is a white
%matter mask
%nii = mni2fs_load_nii(file{1});
%mask=nii.img;
storage_location_1=['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\'];

start_side=1;
sides_num=2;
nparcels=90;
ROIcolors=[distinguishable_colors(nparcels)];
% a = makeColorMap([40, 0, 115]/255,[90, 0, 255]/255, 20);
% a2 = makeColorMap([20, 20, 90]/255,[90, 90, 255]/255, 20);
% a3 = makeColorMap([0, 0, 128]/255,[0, 240, 255]/255, 30);
% a4 = makeColorMap([0, 80, 100]/255,[0, 230, 255]/255, [115, 255, 240]/255, 50);
% a1 = makeColorMap([255, 170, 0]/255,[255, 240, 100]/255, 30);
% ROIcolors_1=[a; a1; a2; a3; a4];
% permuted_colors=randperm(length(ROIcolors_1));
% ROIcolors=ROIcolors_1(permuted_colors, :);
color_map=colormap(ROIcolors);
caxis([1 nparcels]);
verticesIndex2=[];
for MapIndex=1:100
    storage_location_2=['parcellations\parcellation_' num2str(nparcels) '_kmeans_nowhitemattermask\store_' num2str(MapIndex)];
    storage_location=[storage_location_1 storage_location_2];
    
    if ~exist(storage_location, 'dir')
        mkdir(storage_location);
    end

    %metadata = CollectEEGMetadata(data_location,VTK_folder_name);
    for side_index = start_side:sides_num    % Generate frames for L and R
        %[r,c,v] = ind2sub(size(mask),find(mask == 2));
        %         if(side_index == 1)
        %         voxel_cor=[r'+1;c';v'];
        %         else
        %         voxel_cor=[r'-1;c';v'];
        %         end
        
        
        % Determine what side we are working on
        if(side_index == 1)
            side = 'L';
            hem = 'lh';
            
        else
            side = 'R';
            hem = 'rh';
            
        end
        
        surf_fn = fullfile(mni2fs_dir,['\surf\' hem '.surf.gii']);  %%need directory
        fs_surf(side_index) = export(gifti(surf_fn));  %%we need here through line 472 for parcellation
        v_num = size(fs_surf(side_index).vertices, 1);
        fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] *...
            Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
        fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
        surfrender_fn = fullfile(mni2fs_dir,['\surf\' hem '.inflated' num2str(inflationstep) '.surf.gii']);
        inflated_surf(side_index) = export(gifti(surfrender_fn));
        inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] *...
            Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
        inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);
        
        %     cortex_filename{side_index}  =    [side '_cortex.txt'];
        %     smooth_filename{side_index}  =    [side '_smooth.txt'];
        %     cd([data_location '\' VTK_folder_name]);
        %     cortex_v_num(side_index)     =    eval(['metadata.' side '_cortex_v_num']);
        %     cortex_f_num(side_index)     =    eval(['metadata.' side '_cortex_f_num']);
        %     smooth_v_num(side_index)     =    eval(['metadata.' side '_smooth_v_num']);
        %     smooth_f_num(side_index)    =    eval(['metadata.' side '_smooth_f_num']);
        
        % Input and interpret the vertice data for the cortex and smooth
        % surfces
        %     [cortex(side_index).vertices, cortex(side_index).faces] = FormatMeshSurface(cortex_v_num(side_index), cortex_f_num(side_index), data_location, cortex_filename{side_index});
        %     [smooth(side_index).vertices, smooth(side_index).faces] = FormatMeshSurface(smooth_v_num(side_index), smooth_f_num(side_index), data_location, smooth_filename{side_index});
        %
        
        
        %whiteMatter_vertices = ROI_mapping_voxel2vetices2(cortex(side_index).vertices,voxel_cor);
        
        %vertciesIndex=refineBorders2(cortex(side_index),whiteMatter_vertices,nparcels/2+1);
        verticesIndex=refineBorders2(inflated_surf(side_index),nparcels/2);
        
        
        if(side_index==1)

            save([storage_location '\verticesIndex_L_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)], 'verticesIndex')
            verticesIndex2=[verticesIndex2; verticesIndex];
            vertexvaluesL=verticesIndex;
        else

            save([storage_location '\verticesIndex_R_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)], 'verticesIndex')
            verticesIndexR=verticesIndex+nparcels/2;
            %         verticesIndexR=verticesIndex;
            verticesIndex2=[verticesIndex2; verticesIndexR];
            vertexvaluesR=verticesIndexR;
        end
        
    end
    

    
    
    displayElectrodesInflated_v2_parcellation([mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies'],...
        [storage_location] ,...
        {'543MS'},...
        8,...
        1,...
        0,...
        color_map,...
        0,...
        0,...
        4,...
        inflationstep,...
        0,...
        {'_sounds_restricted_CPonly'},...
        vertexvaluesL,...
        vertexvaluesR,...
        nparcels)
     close all
    
end




% function [mesh_vertices, mesh_faces]=FormatMeshSurface(vertices_num, faces_num, data_location, filename)
% % cd(patient_folder);
% % cd ..
% cd(data_location)
% data_location_previous=pwd;
% vertices_range_str=strcat('A1..I',num2str(ceil(vertices_num/3)));
% raw_mesh_vertices=dlmread(fullfile( data_location_previous, filename),'',vertices_range_str);
% faces_range_str=strcat('B',num2str(ceil(vertices_num/3)+1),'..D',num2str(ceil(vertices_num/3)+faces_num));
% raw_mesh_faces=dlmread(fullfile( data_location_previous, filename),'',faces_range_str );
% mesh_faces=raw_mesh_faces+1;
%
% vertices_range=1:ceil(vertices_num/3);
% mesh_vertices=zeros(ceil(vertices_num/3)*3,3);
% mesh_vertices(3*vertices_range-2,1:3)=raw_mesh_vertices(vertices_range, 1:3);
% mesh_vertices(3*vertices_range-1,1:3)=raw_mesh_vertices(vertices_range, 4:6);
% mesh_vertices(3*vertices_range,1:3)=raw_mesh_vertices(vertices_range, 7:9);
% mesh_vertices(ismember(mesh_vertices,zeros(1,3),'rows'),:)=[];
% end

