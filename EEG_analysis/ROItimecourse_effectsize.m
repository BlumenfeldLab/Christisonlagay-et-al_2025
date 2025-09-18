%%%written by Kate Christison-Lagay in April 2021; modified from previous
%%%versions that showed significance in various ways. This one shows it as
%%%a filled in area

p=[];
parcelnumber=80;
analysistypes=[{'CP'}, {'CNP'}, {'Subtraction'}];
%analysistypes=[{'Subtraction'}];
timeoffset=2;

counter=1;
set(0, 'DefaultFigureRenderer', 'painters');

% plottype=[{'all'}, {'collapsed_hemisphere'}];
plottype=['all'];
suffix='4000ms';

%%load timing info
xtick_spacing_filename = ['filter_times_win100ms_shift25ms_full.mat'];
%             xtick_spacing_filename = ['D:\EEG_behavior\hnct_iceeg_scripting\eeg_analysis\filter_times_win100ms_shift25ms.mat'];

load (xtick_spacing_filename)

T_final=(T-timeoffset)*1000;
stimon=find(T_final==0);
[stupidroundingerror firstbin]=min(abs(T_final+500));
[stupidroundingerror neg250bin]=min(abs(T_final+250));
[stupidroundingerror pos250bin]=min(abs(T_final-250));
[stupidroundingerror pos500bin]=min(abs(T_final-500));
[stupidroundingerror pos750bin]=min(abs(T_final-750));
[stupidroundingerror finalbin]=min(abs(T_final-950));




for type=1:length(analysistypes)
    analysistype=char(analysistypes(type));
    path_1=['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_separate_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask_refitted\' analysistype suffix '\'];
    
    path_2=findstr(path_1, '\stats\');
    analysisname_1=path_1(path_2+7:end);
    
    analysisname_2=strrep(analysisname_1, '_', ' ');
    analysisname_2=strrep(analysisname_2, '\', ', ');
    
    
    mappeddatadrive='G:\mnt\Data25\';
    
    %%grab significant vertices
    
    load (['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_common_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask_refitted\' analysistype  suffix '\latesignificantvertices_500ms_common_baseline__baselinestats_filtered_2to2s_' num2str(parcelnumber) 'parcels.mat'])
    load (['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_common_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask_refitted\' analysistype suffix '\latevertexsum_500ms_common_baseline__baselinestats_filtered_2to2s_' num2str(parcelnumber) 'parcels.mat'])
    
    %%load in table of left hemisphere HOAs
    load HOArefinedtable_min200_smoothing6_lh_24-May-2021_version1.mat
    refinedtable_lh=refinedtable;
    
    %%load in table of right hemisphere HOAs
  
        load HOArefinedtable_min200_smoothing6_rh_24-May-2021_version1.mat
    refinedtable_rh=refinedtable;
    
    %%find all the areas that have ROIs
    [C,ia,ib] = intersect(refinedtable_lh.Properties.VariableNames,refinedtable_rh.Properties.VariableNames, 'stable');
    
    if type==1
        plottedareas=listdlg('PromptString', {'Which areas to plot?'}, 'ListString', C);
    end
    
    %%load mesh
    load allmni2fsfaces.mat
   
    
    for a=1:length(plottedareas)
        
        currentarea=char(C(plottedareas(a)));
        if isempty(findobj( 'Type', 'Figure', 'Name', currentarea))
            f=figure('Name', currentarea);
        else
            fh=findobj( 'Type', 'Figure', 'Name', currentarea);
            figure(fh.Number)
        end
        %         h = legend('show','location','best');
        
        hold on
        
        for b=1:2 %%%side stuff
            if isequal(b,1)
                side='lh';
                side_long='left';
                currentarealocation=find(strcmp(refinedtable_lh.Properties.VariableNames, currentarea));
                table_of_faces=refinedtable_lh;
                allfaces=allfaces_lh;
                allvertices=allvertices_lh;
            else
                side='rh';
                side_long='right';
                
                currentarealocation=find(strcmp(refinedtable_rh.Properties.VariableNames, currentarea));
                table_of_faces=refinedtable_rh;
                allfaces=allfaces_rh;
                allvertices=allvertices_rh;
                
            end
            
            currentarea_array=table2array(table_of_faces(:,currentarealocation));
            areas_faces=find(currentarea_array(:,1)~=0);
            
            
            face_to_vertex1=allfaces(areas_faces,:);
            vertices_in_region=unique(face_to_vertex1);
            
            if isequal(side, 'lh')
                
                position_hist=[0.1300    0.2-.15   0.7750     0.08];
                widthfactor=.8;
                
                sign_or_no=find(sum(sign_left(vertices_in_region,:),2)~=0);
                sign_times=find(sum(sign_left(vertices_in_region,:),1)~=0);
                %proportion_sig=sum(sign_left(vertices_in_region,:),1)/size(vertices_in_region,1);
                proportion_sig_1=sum(sign_left(vertices_in_region,:),1)/size(sign_or_no,1);
                total_sig_proportion=size(sign_or_no,1)/size(vertices_in_region,1);
                
                signverts_in_region=vertices_in_region(sign_or_no);
                activity_insignverts_inregion=leftverts(signverts_in_region,:);
                
                switch analysistype
                    case 'CP'
                        colorline=[0 0 1];
                        linestart=24;
                    case 'CNP'
                        
                        colorline=[0 1 1];
                        linestart=22;
                        
                    case 'Subtraction'
                        colorline=[0 0 .6];
                        linestart=22;
                        signstar=1.2;
                end
                
                
                
            else
                position_hist=[0.1300    0.2-.15    0.7750     0.08];
                widthfactor=.5;
                %                         position_hist=[0.1300    0.3369-.25    0.7750    0.08];
                
                sign_or_no=find(sum(sign_right(vertices_in_region,:),2)~=0);
                sign_times=find(sum(sign_right(vertices_in_region,:),1)~=0);
                %proportion_sig=sum(sign_right(vertices_in_region,:),1)/size(vertices_in_region,1);
                proportion_sig_1=sum(sign_right(vertices_in_region,:),1)/size(sign_or_no,1);
                total_sig_proportion=size(sign_or_no,1)/size(vertices_in_region,1);
                
                signverts_in_region=vertices_in_region(sign_or_no);
                activity_insignverts_inregion=rightverts(signverts_in_region,:);
                switch analysistype
                    case 'CP'
                        colorline=[1 0 0];
                        linestart=23;
                        
                    case 'CNP'
                        
                        colorline=[1 0 1];
                        linestart=21;
                        
                        
                    case 'Subtraction'
                        colorline=[.6 0 0];
                        linestart=21;
                        
                end
            end
%             currentarea
%             min(find(proportion_sig_1>0))
            proportion_sig=[0*ones(1,firstbin) proportion_sig_1];
            
            %                     position1=[0.1300    0.3369    0.7750    0.5881];
            position1=[0.1300    0.2    0.7750    0.7];
                                   signstar=1.2;
 
            %                     subplot('Position', position1);
            hold on;
            
            %             linenumber=a*type+(type-1+b-1);
            linenumber=counter;
            zero_rows=sum(activity_insignverts_inregion,2)==0;
            activity_insignverts_inregion(zero_rows,:)=[];
            avactivity=nanmean(activity_insignverts_inregion,1);
            
            avactivity2=avactivity(:, firstbin:finalbin);
            
            avactivity_baseline=mean(avactivity(:,firstbin:stimon-1));
            
            
            [peakpoststimactivity peakpositionbin]=max(avactivity(:,stimon:finalbin));
            
            peaktime=T_final(stimon+peakpositionbin-1);
            
            effectsize=peakpoststimactivity-avactivity_baseline;
            A = exist('table2');
            if isequal(A, 0)
                
                table2=table;
                table2.Analysis=analysistype;
                table2.Area=currentarea;
                table2.Side=side;
                table2.EffectSize=effectsize;
                table2.PeakTime=peaktime;
                
            else
                dataforeffectsize={analysistype, currentarea, side, effectsize, peaktime};
                
                table2=[table2;dataforeffectsize];
            end
        end
        
    end
    
end

writetable(table2,'effectsize.csv')