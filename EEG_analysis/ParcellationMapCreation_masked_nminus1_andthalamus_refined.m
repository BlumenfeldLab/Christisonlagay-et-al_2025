mni2fs_dir='D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\mni2fs-master';
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mappeddatadrive='V:';


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
load allmni2fsfaces.mat

start_side=1;
sides_num=2;
nparcels=80;
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
    storage_location_2=['parcellations\parcellation_' num2str(nparcels) '_kmeans_mask_refined\store_' num2str(MapIndex)];
        %storage_location_2=['parcellations\parcellation_'
        %num2str(nparcels) '_kmeans_mask_refitted\store_'
        %num2str(MapIndex)];  %%folder previously named 'refitted'; updated
        %folder location as to not write over data in case of code being
        %run again.

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
            load HOArefinedtable_min200_smoothing6_lh_15-Apr-2021_version6.mat
            allfaces=allfaces_lh;
            allvertices=allvertices_lh;
            


        else
            side = 'R';
            hem = 'rh';
            load HOArefinedtable_min200_smoothing6_rh_09-Mar-2021_version3.mat
            allfaces=allfaces_rh;
            allvertices=allvertices_rh;
        end
        
        
        whitematter_faces_1=refinedtable.background;
        whitematter_faces=find(whitematter_faces_1(:,1)~=0);
        
        
%         if (side_index == 2) %%%Thalamus only include for right hemisphere
%             thalamic_faces_1=refinedtable.MNI_Thalamus_R_roi;
%             thalamic_faces=find(thalamic_faces_1(:,1)~=0);
%             sharedfaces=find(ismember(whitematter_faces,thalamic_faces));
%             whitematter_faces(sharedfaces)=[];
%         end
        
        face_to_vertex1=allfaces(whitematter_faces,:);
        whitematter_vertices=unique(face_to_vertex1);
        
        
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
        
        %         inflated_surf(side_index).vertices(whitematter_vertices,:)=[];
        
        %         verticesIndex=refineBorders2(inflated_surf(side_index),nparcels/2);
        if (side_index==1)
            verticesIndex=refineBorders3(inflated_surf(side_index),whitematter_vertices,nparcels/2+1);
        else
            verticesIndex=refineBorders3(inflated_surf(side_index),whitematter_vertices,nparcels/2);
            
            thalamic_faces_1=refinedtable.Thalamus;
            thalamic_faces=find(thalamic_faces_1(:,1)~=0);
            face_to_vertex1=allfaces(thalamic_faces,:);
            thalamic_vertices=unique(face_to_vertex1);
            max_vertexIndex=max(unique(verticesIndex));
            verticesIndex(thalamic_vertices)=max_vertexIndex+1;

        end
        
        if(side_index==1)

            save([storage_location '\verticesIndex_L_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)], 'verticesIndex')
            verticesIndex2=[verticesIndex2; verticesIndex];
            vertexvaluesL=verticesIndex;
        else

            save([storage_location '\verticesIndex_R_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)], 'verticesIndex')
            verticesIndexR=verticesIndex+nparcels/2;
            verticesIndexR(verticesIndexR==nparcels/2)=0;
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




