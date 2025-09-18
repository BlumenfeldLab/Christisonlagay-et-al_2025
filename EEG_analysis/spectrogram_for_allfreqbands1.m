clear all


mappeddatadrive='G:\mnt\Data25\';


%%       load patients, and suffixes
load([mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_21_04_15.mat']);


sorted_suffixes=eval('sorted_trials_2to2_r_suffixes');
noisytrials_suffixes=eval(['composite_rejection_recentered_r_suffixes']);


xtick_spacing_filename = ['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\filter_times_win100ms_shift25ms_full.mat'];

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
        
trialtypes = {'BC_trials','WG_trials', 'blank_TN_trials'};

for p =[1:(length(patients)-1)] 
    display(p)
    %% aggregate mean power zscore for each patient
    patient_folder2= [data_location2, '/', patients{p} '/icEEG/HNCT Auditory ID Combined'];
    cd(patient_folder2)
    
  
        load('labels_all_gray_matter.mat')
   

    
    load(['sorted_trials' sorted_suffixes{p}])
    load(['composite_rejection_recentered' noisytrials_suffixes{p}])
    fs1=fs(p);
  
    
    

    for i = 1:length(labels)
        disp(['Calculating channel ' num2str(i) ' of ' num2str(length(labels)) ' for patient ' num2str(p)])
        %try
        
        %channel_epochs_alltrialtype=[];
        for j = 1:2
           
            eval(['channel_epochs = squeeze(' trialtypes{j} '(i,:,:));']);
            
            %% This section is only for subject 193
            % figure out the number of trials of that type

            
            % set the indicies of rejected trials
            eval(['noisy_trials_idx = noisy_' trialtypes{j} '_idx{i};']);
            
            % change rejected trials NaNs
            channel_epochs(:,noisy_trials_idx) = [];
            % end
            if size(channel_epochs,2) > 1
                S = [];
                for k = 1:size(channel_epochs,2)
                    
              
                    data = squeeze(channel_epochs(:,k));
                    
                    %downsample so we're all on the same frequency
                    datadown =downsample(data, ceil(fs1/1024));
                    
                    %cut it down to -1.5 to 1.5 second
                    datadown =datadown(513:3584);

                    %data = [data_or;data_or];
                    [wt, f] = cwt(datadown, 'amor', 1024, 'VoicesPerOctave', 20);
                    psd = abs(wt.^2)./f;
                    for d_data=1:size(psd,1)
                       
                        if ~isequal(size(psd(d_data,:),2)/1024, ceil(size(psd(d_data,:),2)/1024))
                            psd_prep1=ceil(size(psd(d_data,:),2)/1024);
                            psd_prep2=NaN*ones(1, psd_prep1*1024);
                            psd_prep2(1,1:size(psd,2))=psd(d_data,:);
                            psd_1=psd_prep2;
                        else
                            psd_1=psd(d_data,:);
                        end
                         size1=32;
                        size2=size(psd_1, 2)/size1;
                        data1(d_data,:)=nanmean(reshape(psd_1,[size1,size2]),1);
                    end
                    
                    S(:,:,k) = data1(:,(512/32)+1:size2-(512/32));
                end
                meanpower_average(i,j,:,:) = nanmean(S,3);
                
                
                %%i=# electrodes
                %%j=# conditions
                %%3rd dimension=power at frequency
                %%4th dimension=time bin
            
                
                
            end
        end
    end
    
    
   
    save('meanpower_average_wavelet_nonoverlaped2.mat','meanpower_average','-v7.3');
    clearvars channel_epochs data1 S meanpower_average
end
