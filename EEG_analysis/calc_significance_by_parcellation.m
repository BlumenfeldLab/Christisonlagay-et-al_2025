%%% runs after extract_vertices_by_parcellation
%%Written by K Christison-Lagay April 2021

%%
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

%% otherstuff
baseline_end_real=78;
baseline_end=38;

xtick_spacing_filename = ['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\filter_times_win100ms_shift25ms_full.mat'];
load (xtick_spacing_filename)
final_bin=length(T); %first 500 ms
runduration=[num2str(4000) 'ms'];

%% number of parcels
nparcels=[ 80 ];


for overlays=[1:3]
    if isequal(overlays, 1)
        overlay='CP';
    elseif isequal(overlays, 2)
        overlay='CNP';
    else
        overlay='Subtraction';
    end
    
    start_pt=1;
    
    for maps=start_pt:100
        
        
        inflationstep=5;
        createMontage = 0;
        views = 4;
        MapIndex=maps;
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
        
        load ([save_location '\allverticevalues_' num2str(maps) '.mat'])
        
        
        significant_vertices_filename=[save_location '\significant_vertices_map' num2str(MapIndex) '.mat'];
        disp(['Running stats for ' laterality_label ' ' overlay ' using map ' num2str(MapIndex) ' with ' num2str(nparcels) ' parcels, with a ' num2str(lengthofbaseline) 'ms ' baseline_name ' baseline '  statsinclusionverbose])
        
        
        
        
        % change to local version
        mni2fs_dir = 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\mni2fs-master';
        mni_nifti_path = fullfile(mni2fs_dir, 'MNI_T1_1mm_stripped.nii');
        load([mni2fs_dir, filesep, 'surf', filesep, 'transmats.mat']);
        mnit1_1mm_stripped = mni2fs_load_nii(mni_nifti_path);
        Tmni = mnit1_1mm_stripped.transform;
        
        
        
        
        sides_num = 2;
        for statstype=2
            total_num_frames = numel(T);

            if statstype==1
                statstiming='early';
                finalbin=97;
                total_num_frames2=length(baseline_start_real:finalbin);
                

            else
                statstiming='late';
                finalbin=137;

                total_num_frames2=length(baseline_start_real:finalbin);

            end
                
            significant_vertices_filename=[save_location '\' statstiming 'significant_vertices_map' num2str(MapIndex) '.mat'];
            
            
            
            %% Generating all frames for each side of the brain
            for side_index = 1 %:2    % Generate frames for L and R
                
                % Determine what side we are working on
                if side_index == 1
                    side = 'L';
                    hem = 'lh';
                    verticesIndex=verticesIndexL;
                    chan_hood=chan_hoodL;
                    ROI_electrode_values= vertex_valuesL;
                else
                    side = 'R';
                    hem = 'rh';
                    verticesIndex=verticesIndexR;
                    chan_hood=chan_hoodR;
                    ROI_electrode_values= vertex_valuesR;
                    
                end
                
                % Setup surface variables
                surf_fn = fullfile(mni2fs_dir,['/surf/' hem '.surf.gii']);  %%need directory
                fs_surf(side_index) = export(gifti(surf_fn));  %%we need here through line 472 for parcellation
                v_num = size(fs_surf(side_index).vertices, 1);
                f_num = size(fs_surf(side_index).faces, 1);
                
                fs_surf(side_index).vertices = [fs_surf(side_index).vertices, ones(v_num, 1)] *...
                    Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
                fs_surf(side_index).vertices = fs_surf(side_index).vertices(:, 1:3);
                surfrender_fn = fullfile(mni2fs_dir,['/surf/' hem '.inflated' num2str(inflationstep) '.surf.gii']);
                inflated_surf(side_index) = export(gifti(surfrender_fn));
                inflated_surf(side_index).vertices = [inflated_surf(side_index).vertices, ones(v_num, 1)] *...
                    Tfstovox_rcor' * Trsvoxtomni_rcor' / Tmni';
                inflated_surf(side_index).vertices = inflated_surf(side_index).vertices(:, 1:3);
                
                
                
                
                % Initialize empty vertex_values variable
                
                
                
                
                
                
                sigclusters=zeros(size(inflated_surf(side_index).vertices,1),total_num_frames2);
                side_index
                disp(['>> Calculating significant parcels for ' side ' Cortex']);  % Output Status
                
                % sum the count of vertices across patients
                vertex_values_sum = zeros(size(inflated_surf(side_index).vertices,1),1);
                
                
                Baseline=repmat(nanmean(ROI_electrode_values(:,(baseline_start_real:baseline_end_real),:),2),1,total_num_frames2,1);
                
                [pval, t_orig, clust_info, seed_state, est_alpha]=clust_perm1_iceeg_sumt(ROI_electrode_values(:,baseline_start_real:finalbin,:)-Baseline,chan_hood);
                pos_clust_id=clust_info.pos_clust_ids;
                neg_clust_id=clust_info.neg_clust_ids;
                
                %
                
                posn=find(clust_info.pos_clust_pval<=0.05);
                negn=find(clust_info.neg_clust_pval<=0.05);
                
                for frame_index=1:total_num_frames2
                    for pcc=1:length(posn)
                        idpos=find(pos_clust_id(:,frame_index)==posn(pcc))
                        
                        if(isempty(idpos))
                        else
                            
                            for nnn=1:length(idpos)
                                
                                vertex_values_sum(find(verticesIndex==idpos(nnn)))=8;
                                sigclusters(find(verticesIndex==idpos(nnn)),frame_index)=1;
                            end
                        end
                        
                    end
                end
                
                for pcc=1:length(negn)
                    idneg=find(neg_clust_id(:,frame_index)==negn(pcc))
                    if(isempty(idneg))
                    else
                        for nnn=1:length(idneg)
                            vertex_values_sum(find(verticesIndex==idneg(nnn)))=-4;
                            sigclusters(find(verticesIndex==idneg(nnn)),frame_index)=1;
                            
                            
                        end
                        
                        
                    end
                end
                
            end
            
            
            
            
            
        
            if side_index == 1
                significantvertex_valuesL= sigclusters;
                size(sigclusters)
                
                if exist(significant_vertices_filename, 'file')
                    save(significant_vertices_filename,'significantvertex_valuesL','-append')
                else
                    save(significant_vertices_filename,'significantvertex_valuesL', '-append')
                end
                
            elseif side_index == 2
                significantvertex_valuesR= sigclusters;
                size(sigclusters)
                
                if exist(significant_vertices_filename, 'file')
                    save(significant_vertices_filename,'significantvertex_valuesR','-append')
                else
                    save(significant_vertices_filename,'significantvertex_valuesR', '-append')
                end
            end
            
            
            
        end
        
        
    end
end
