%%%updated by KCL from Jiajia's EVENTS_6rois-Frequency.m

% Description: this code takes the confirmed perceived/not percieved trials. It gets rid of the
% rejected trials and blank trials and also the trials >20mean/Devidation then finds the mean power increases
% for the frequency bands specified and plots them. Once complete it saves
% the mean power trace for later use.

% Inputs
% statistic method label- SL = 1: pairred ttest; 2: Nonparamettric Wilcoxon signed rank test
% BC_trials_left or BC_trials_right = EEG data of confirmed perceived
%           trials in a 3D matrix: #electrodes x #samples x #trials
% noisy_BC_trials_left_idx (or right) = array containing indices of trials
%           for each electrode that should be rejected
% window_size = size of window to perform fourier transform on
% sliding_window = how many units to move the window by incrementally
% frequency band = matrix of dimensions #bands x cutoffs that you want to
%           analyze (ie. [0 4; 4 12; 40 115])

% Output
% meanpower_traces = matrix containing the values for the average power
%           zscores trace of each electrode organized by laterality,
%           dimensions: #electrodes x #frequency bands x #bins x contra(1) vs ipsi(2)



%These are the three frequency bands we will be analyzing. 40-115 =
%gamma band

mappeddatadrive='G:\mnt\Data25\';

%These are the three frequency bands we will be analyzing. 40-115 =
%gamma band

%%       load patients, and suffixes
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




%baseline_sample = [1024-512:1024];
baseline_sample = [32-16:32];

%%load the vertices for each patient's vertices
load([mappeddatadrive '\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\Group\Inflated Surface\tmp\Electrode_Distribution\group_elec_vertices.mat'])
%%                        set the intials for meanpower for L+R/L/R     %%%

% meanpower_frequency=NaN(length(ROI_network),length(patient),183,120);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%                        start the loop from each ROI
load([ 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\ROIvertices.mat'])
load([ 'D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\electrodesinROIs.mat'])

log_or_linear='log'; %'linear' 'log'

ROIareas=[{'ParsOpercularis'} ] %{'ComboAuditory'}, {'RostralInsula'}, {'CaudalMFG'}, {'Supramarginal'}, {'ParsOpercularis'}];
for ROIs=1:length(ROIareas)
    currentarea=char(ROIareas(ROIs));
    eval(['currentareaelectrodes=' currentarea ';'])
    for type=1:2
        meanpower_frequency=[];
        for p_index =[1:31] 
            display(p_index)
            %% aggregate mean power zscore for each patient
            patient_folder2= [data_location2, '/', patients{p_index} '/icEEG/HNCT Auditory ID Combined'];
            cd(patient_folder2)
            
            pt_electrodesROI=cell2mat(currentareaelectrodes(p_index));
          
            load('meanpower_average_wavelet_nonoverlaped2.mat')
            
                       
            
            
            test2=[];
            test2=reshape((meanpower_average(pt_electrodesROI,type,:,:)),[size(pt_electrodesROI,1), size(meanpower_average,3), size(meanpower_average,4)]);
            test2commonbackground_1=nanmean(meanpower_average(pt_electrodesROI,:,:,:),2);
            test2commonbackground=reshape(test2commonbackground_1, [size(pt_electrodesROI,1), size(meanpower_average,3), size(meanpower_average,4)]);
            %         test2CP=squeeze(nanmean(meanpower_average(:,1,:,:),2));
            %         test2CNP=squeeze(nanmean(meanpower_average(:,2,:,:),2));
            
            %                 testbaseline=squeeze(nanmean(meanpower_average(:,:,:,:),2));
            
            %test2=squeeze(nanmean(meanpower_average(i,:,:,:),2));
            %  test2=squeeze(meanpower_average(:,:,:,:));
            %    test2 = squeeze(nanmean(test2,1));
            test2 = log(test2);
            test2commonbackground= log(test2commonbackground);
            for i = 1:size(test2,1)
                for j = 1:size(test2,2)
                    %baseline_sample = [1:65-1,193+1:256];
                                    test2(i,j,:) = squeeze(test2(i,j,:)) - nanmean(squeeze(test2(i,j,baseline_sample)));
                    
%                     test2(i,j,:) = squeeze(test2(i,j,:)) - nanmean(squeeze(test2commonbackground(i,j,baseline_sample)));
                    
                    %test2(i,j,:) = squeeze(test2(i,j,:)) - nanmean(squeeze(test2(i,j,1:64)));
                    
                end
            end
            % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
            %                  smooth the loaded meanpower azores in each patient  %%%
            
            %% mean
%             meanpower_frequency(p_index,:,:)=nanmean(test2(:,:,:),1);
            meanpower_frequency=[meanpower_frequency; test2];

            % meanpower_task(l,p_index)=squeeze(nanmean(nanmean(meanpower_traces(ROI_CH_Index{l,p_index},1,task_sample),1),3));
            
        end
        
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        %%                        plotting
        
        load('A:\Human CPT icEEG Project\icEEG Analysis\Group Analysis\CPT ERP movies\spectrogram_times_64bins_2s.mat')
        load([mappeddatadrive 'HNCT\icEEG Analysis\Analysis\EEG_behavior\404LL\icEEG\HNCT Auditory ID Combined\frequencies_forallfreqs.mat'])
        
        %T=1/1024-1:1/1024:2048/1024-1;
        %     load('D:\Human CPT icEEG Project\icEEG Analysis\Group Analysis\ROI Analysis\meanpower_frequency_Event_nonoverlap');
        
        
        %%64=number of timebins included
        %%first # here is how many frequencies to be included
        
        meanpower_frequency_interpolate=NaN(175,64);
        interpolate_point=[1:1:175];
        
        for i= 1:1
            meanpower_frequency_interpolate=NaN(175,64);  
            
            
            %meanpower_frequency dimensions: # electrodes; # of f; 
            test1=squeeze(nanmean(meanpower_frequency(:,:,:),1));
           switch log_or_linear
               case 'linear'

            for  j=1:size(test1,2)
                test2=test1(:,j);
                
                 
                %interp1 here does a linear interpolation of power by rescaling f to be
                %linear. 
                
                vq1 = interp1(f,test1(:,j),(interpolate_point),'spline');

                meanpower_frequency_interpolate(:,j)=vq1;
            end
               case 'log'     

                    for  j=1:size(test1,2)
                test2=test1(:,j);
               

          %skip this to keep it in octave/log space
                vq1 =test1(:,j);

                meanpower_frequency_interpolate(:,j)=vq1;
                    end
           end
        
        end
        
        
        
        
      switch log_or_linear
          
          case 'linear'
        %%use this one if it's in linear frequency space
%         test1=meanpower_frequency_interpolate(1:175,:);
figure;
        set(gcf, 'color', [1 1 1]);
        test1=meanpower_frequency_interpolate(1:175,:);
        imagesc(T,[1:175],test1);
        xlabel('Time (s)','FontSize',16)
        ylabel('Frequency (Hz)','FontSize',16)
        set(gca,'YDir','normal')

        if type==1
            condition='Perceived';
        else
            condition='Not Perceived';
        end

        title(['Log 0-140 Hz Power  ', currentarea,'  ' condition]); hold on;
        color_axis = caxis;

        colormap jet
        colorbar
        caxis([-0.2 .2])
        xlim([-.5 1])
        ylim([0 140])
        xline(0)
        xline(0.075)

          
          
          case 'log'
              
        %%use this one if it's still in octave space
       test1=flipud(meanpower_frequency_interpolate(1:175,:));
figure;
        set(gcf, 'color', [1 1 1]);
       
        imagesc(T,1:175,test1);
        xlabel('Time (s)','FontSize',16)
        ylabel('Frequency (Hz)','FontSize',16)
        set(gca,'YDir','normal')
        set(gca, 'YScale', 'log')

        if type==1
            condition='Perceived';
        else
            condition='Not Perceived';
        end

        title(['Log 0-140 Hz Power  ', currentarea,'  ' condition]); hold on;
        color_axis = caxis;

        colormap jet
        colorbar
        caxis([-0.2 .2])
        xlim([-.5 1])
%         ymin=min(log(f));
%         ymax=max(log(f));
        ylim([0 140])
        xline(0)
        xline(0.075)

 ticks = [1  10  100 140];

set(gca,'ytick',ticks)
      
      end
        
    end
    
end


