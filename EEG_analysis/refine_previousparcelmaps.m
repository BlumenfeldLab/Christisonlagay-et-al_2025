mni2fs_dir='D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\mni2fs-master';
load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mappeddatadrive='V:';


mni_nifti_path = fullfile(mni2fs_dir, '\MNI_T1_1mm_stripped.nii');

mnit1_1mm_stripped = niftiread(mni_nifti_path);

load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
Tmni = mnit1_1mm_stripped.transform;

inflationstep=5;


storage_location_1=['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\'];
load allmni2fsfaces.mat

start_side=1;
sides_num=2;
nparcels=80;
ROIcolors=[distinguishable_colors(nparcels)];


color_map=colormap(ROIcolors);
caxis([1 nparcels]);
verticesIndex2=[];
for MapIndex=1:100
    storage_location_2=['parcellations\parcellation_' num2str(nparcels) '_kmeans_mask_refitted\store_' num2str(MapIndex)];
    storage_location=[storage_location_1 storage_location_2];
    
    previousstorage_location=[storage_location_1 'parcellations\parcellation_' num2str(nparcels) '_kmeans_mask\store_' num2str(MapIndex)];

    if ~exist(storage_location, 'dir')
        mkdir(storage_location);
    end

    for side_index =start_side:sides_num    % Generate frames for L and R
    
        
        % Determine what side we are working on
        if(side_index == 1)
            side = 'L';
            hem = 'lh';
            load HOArefinedtable_min200_smoothing6_lh_15-Apr-2021_version6.mat
            allfaces=allfaces_lh;
            allvertices=allvertices_lh;
            load ([previousstorage_location '\verticesIndex_L_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)])
 


        else
            side = 'R';
            hem = 'rh';
            load HOArefinedtable_min200_smoothing6_rh_09-Mar-2021_version3.mat
            allfaces=allfaces_rh;
            allvertices=allvertices_rh;
            load ([previousstorage_location '\verticesIndex_R_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)])

        end
        
        
        whitematter_faces=find(refinedtable.background~=0);
        whitematter_verts1=allfaces(whitematter_faces,:);
        whitematter_verts=unique(whitematter_verts1);
        
      
        if (side_index == 2) %%%Thalamus only include for right hemisphere
            thalamic_faces=find(refinedtable.Thalamus~=0);
            thalamic_verts1=allfaces(thalamic_faces,:);
            thalamic_verts=unique(thalamic_verts1);
        end
        
        
        if (side_index==1)
            verticesIndex(whitematter_verts)=0;
        else
            allthalamicparcelIDs=unique(verticesIndex(thalamic_verts));

            verticesIndex(whitematter_verts)=0;
            verticesIndex(thalamic_verts)=max(allthalamicparcelIDs);

        end
        
        if(side_index==1)

            save([storage_location '\verticesIndex_L_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)], 'verticesIndex')
            vertexvaluesL=verticesIndex;
        else

            save([storage_location '\verticesIndex_R_' num2str(nparcels) 'ParcellationMap_' num2str(MapIndex)], 'verticesIndex')
            vertexvaluesR=verticesIndex;
        end
%         
    end
    

    
    
    displayElectrodesInflated_parcellation_June2021update([mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies'],...
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




