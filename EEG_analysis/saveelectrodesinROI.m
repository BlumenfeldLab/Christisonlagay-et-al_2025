mappeddatadrive='V:';
load ([ mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\31 Patients Auditory ID 157 Bins Sounds Restricted by CP Accuracy 100ms\filter_baseline5to0_common_r\Gamma_power_Subtraction\group_elec_vertices.mat']);




load([mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_21_04_15.mat']);


sorted_suffixes=eval('sorted_trials_2to2_r_suffixes');
noisytrials_suffixes=eval(['composite_rejection_recentered_r_suffixes']);


frequency_bands = [0 4; 4 12; 40 115];
%frequency_bands = [40 115];
load('A:\Human CPT icEEG Project\icEEG Analysis\Group Analysis\CPT ERP movies\spectrogram_times_64bins_2s.mat')

data_location=[mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies'];
data_location2=[mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior'];
save_location_1 = ['D:\HNCT\HNCT Auditory\Analysis Code\spectrogram'];
save_location_2_1= ['500ms_common_baseline_'];
save_location_2_2=['100msbins75msoverlap'];
save_location_2_3=['filtered_2to2s'];
save_location_2=[save_location_2_1 save_location_2_3];
save_location=[save_location_1 save_location_2_2 '\' save_location_2];


if ~exist(save_location, 'dir')
    mkdir(save_location)
end

%%load the vertices for each patient's vertices
load([ 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\ROIvertices_forauditory.mat'])

for p_index =[1:5, 7:31] %5, 7:(length(patients)-1)]
    display(p_index)
    %% aggregate mean power zscore for each patient
    patient_folder3= [data_location2, '/HNCT ERP movies/', patients{p_index}];
    cd(patient_folder3)
    
    load R_MontageMap.mat
    load L_MontageMap.mat
    load MontageMap.mat
    
    left_electrodelocations=cell2mat(group_elec_vertices(p_index,1));
    right_electrodelocations=cell2mat(group_elec_vertices(p_index,2));
    
    
    
    %%are left electrodes in an ROI
    Supramar_lh=find(ismember(left_electrodelocations, Supramarginal_lh));
    RostIns_lh=find(ismember(left_electrodelocations, RostralInsula_lh));
    Aud_lh=find(ismember(left_electrodelocations, ComboAuditory_lh));
    MFG_lh=find(ismember(left_electrodelocations, CaudalMFG_lh));
    ParsOp_lh=find(ismember(left_electrodelocations, ParsOperculum_lh));
    Heschlsgyrus_lh=find(ismember(left_electrodelocations, Heschls_lh));

    
    
    
    
    %%are right electrodes in an ROI
    
    Supramar_rh=find(ismember(right_electrodelocations, Supramarginal_rh));
    RostIns_rh=find(ismember(right_electrodelocations, RostralInsula_rh));
    Aud_rh=find(ismember(right_electrodelocations, ComboAuditory_rh));
    MFG_rh=find(ismember(right_electrodelocations, CaudalMFG_rh));
    ParsOp_rh=find(ismember(right_electrodelocations, ParsOperculum_rh));
        Heschlsgyrus_rh=find(ismember(right_electrodelocations, Heschls_rh));

    
    %%what are the electrodes' position in MontageMap?
    if ~isempty(L_MontageMap)&~isempty(R_MontageMap)
        Supramarginal(p_index)={[L_MontageMap(Supramar_lh,1); R_MontageMap(Supramar_rh,1)]};
        RostralInsula(p_index)={[L_MontageMap(RostIns_lh,1); R_MontageMap(RostIns_rh,1)]};
        ComboAuditory(p_index)={[L_MontageMap(Aud_lh,1); R_MontageMap(Aud_rh,1)]};
        CaudalMFG(p_index)={[L_MontageMap(MFG_lh,1); R_MontageMap(MFG_rh,1)]};
        ParsOpercularis(p_index)={[L_MontageMap(ParsOp_lh,1); R_MontageMap(ParsOp_rh,1)]};
        Heschls(p_index)={[L_MontageMap(Heschlsgyrus_lh,1); R_MontageMap(Heschlsgyrus_rh,1)]};

        
        
    elseif  ~isempty(L_MontageMap)&isempty(R_MontageMap)
        Supramarginal(p_index)={[L_MontageMap(Supramar_lh,1)]};
        RostralInsula(p_index)={[L_MontageMap(RostIns_lh,1)]};
        ComboAuditory(p_index)={[L_MontageMap(Aud_lh,1)]};
        CaudalMFG(p_index)={[L_MontageMap(MFG_lh,1)]};
        ParsOpercularis(p_index)={[L_MontageMap(ParsOp_lh,1)]};
        Heschls(p_index)={[L_MontageMap(Heschlsgyrus_lh,1)]};

    elseif isempty(L_MontageMap)&~isempty(R_MontageMap)
        Supramarginal(p_index)={[R_MontageMap(Supramar_rh,1)]};
        RostralInsula(p_index)={[R_MontageMap(RostIns_rh,1)]};
        ComboAuditory(p_index)={[R_MontageMap(Aud_rh,1)]};
        CaudalMFG(p_index)={[R_MontageMap(MFG_rh,1)]};
        ParsOpercularis(p_index)={[R_MontageMap(ParsOp_rh,1)]};
        Heschls(p_index)={[R_MontageMap(Heschlsgyrus_rh,1)]};

    end
    
end

% filename='electrodesinROIs2.mat';
% save(filename, 'Supramarginal', 'RostralInsula', 'ComboAuditory','CaudalMFG', 'ParsOpercularis', 'Heschls');