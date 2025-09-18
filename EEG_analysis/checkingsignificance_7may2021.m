
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
%%just a quick script to see if anything is showing up sig on the right
%%hemisphere

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


for overlays=[2]
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
end