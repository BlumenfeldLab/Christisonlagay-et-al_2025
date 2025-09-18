% significantvertices_filename=
% vertexvalue_filename=
p=[];
parcelnumber=80;
% analysistypes=[{'CP'}, {'CNP'}, {'Subtraction'}];
analysistypes=[ {'Subtraction'}];


counter=1;
set(0, 'DefaultFigureRenderer', 'painters');

% plottype=[{'all'}, {'collapsed_hemisphere'}];
plottype=['all'];

switch plottype
    case 'all'
        
        for type=1:length(analysistypes)
            analysistype=char(analysistypes(type));
            
            path_1=['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_separate_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask\' analysistype '\'];
            
            path_2=findstr(path_1, '\stats\');
            analysisname_1=path_1(path_2+7:end);
            
            analysisname_2=strrep(analysisname_1, '_', ' ');
            analysisname_2=strrep(analysisname_2, '\', ', ');
            
            
    mappeddatadrive='G:\mnt\Data25\';
            
            %%grab significant vertices
            load (['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_common_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask\' analysistype '\significantvertices_500ms_common_baseline__baselinestats_filtered_2to2s_' num2str(parcelnumber) 'parcels.mat'])
            load (['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_common_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask\' analysistype '\vertexsum_500ms_common_baseline__baselinestats_filtered_2to2s_' num2str(parcelnumber) 'parcels.mat'])
            
            %%load in table of left hemisphere HOAs
            load HOArefinedtable_min200_smoothing6_lh_15-Apr-2021_version6.mat
            refinedtable_lh=refinedtable;
            
            %%load in table of right hemisphere HOAs
            %             load HOArefinedtable_min200_smoothing6_rh.mat
            %             load HOArefinedtable_min200_smoothing6_rh_05-Feb-2021.mat
            load HOArefinedtable_min200_smoothing6_rh_09-Mar-2021_version3.mat
            
            refinedtable_rh=refinedtable;
            
            %%find all the areas that have ROIs
            [C,ia,ib] = intersect(refinedtable_lh.Properties.VariableNames,refinedtable_rh.Properties.VariableNames, 'stable');
            
            if type==1
                plottedareas=listdlg('PromptString', {'Which areas to plot?'}, 'ListString', C);
            end
            
            %%load mesh
            load allmni2fsfaces.mat
            
            %%load timing info
                        xtick_spacing_filename = [mappeddatadrive 'HNCT\icEEG Analysis\Analysis\EEG_behavior\hnct_iceeg_scripting\eeg_analysis\filter_times_win100ms_shift25ms.mat'];
%             xtick_spacing_filename = ['D:\EEG_behavior\hnct_iceeg_scripting\eeg_analysis\filter_times_win100ms_shift25ms.mat'];
            
            load (xtick_spacing_filename)
            
            
            
            for a=1:length(plottedareas) %[1:8, 10:size(table_of_faces,2)]
                
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
                        proportion_sig=sum(sign_left(vertices_in_region,:),1)/size(sign_or_no,1);
                        
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
                                signstar=1.1;
                        end
                        
                        
                        
                    else
                        position_hist=[0.1300    0.2-.15    0.7750     0.08];
                        widthfactor=.5;
                        %                         position_hist=[0.1300    0.3369-.25    0.7750    0.08];
                        
                        sign_or_no=find(sum(sign_right(vertices_in_region,:),2)~=0);
                        sign_times=find(sum(sign_right(vertices_in_region,:),1)~=0);
                        %proportion_sig=sum(sign_right(vertices_in_region,:),1)/size(vertices_in_region,1);
                        proportion_sig=sum(sign_right(vertices_in_region,:),1)/size(sign_or_no,1);
                        
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
                                                    signstar=1.1;

                    %                     position1=[0.1300    0.3369    0.7750    0.5881];
                    position1=[0.1300    0.2    0.7750    0.7];
                    
                    subplot('Position', position1); hold on;
                    
                    %             linenumber=a*type+(type-1+b-1);
                    linenumber=counter;
                    zero_rows=sum(activity_insignverts_inregion,2)==0;
                    activity_insignverts_inregion(zero_rows,:)=[];
                    avactivity=nanmean(activity_insignverts_inregion,1);
                    if ~isequal(sum(~isnan(avactivity)),0)
                        
                        stderrorofmean= std( activity_insignverts_inregion ) / sqrt( length( activity_insignverts_inregion ));
                        %         errorbar(1:77, avactivity, stderrorofmean); hold on
                        
                        options.error='sem';
                        options.color_area=colorline;
                        
                        %                         scatter(sign_times,linestart*ones(size(sign_times)), 'filled', 'MarkerFaceColor', colorline);
                        
                        
                        plot_areaerrorbar(activity_insignverts_inregion, options);
                        h(linenumber)= plot(avactivity, 'Color', colorline, 'DisplayName', [currentarea ' ' analysistype  ' ' side_long]);
                        times=(T-1)*1000;
                        xticks([19 29 39 49 59 69]);
                        xticklabels({'-500', '-250', '0', '250', '500', '750'});
                        %                         xlabel('ms')
                        %                         axis([19 77.5 -3 33])
%                         axis([19 77.5 -10 18])
                        
%                                                 axis([19 77.5 -8 20])
                                                axis([19 77.5 -20 35])

                        xline(42, '--')
                        xline(39, '--')
                        ylabel('weighted zscore')
                        
                        
                        
                        
                        subplot('Position', position_hist); hold on
                        
                        
                        b=bar(proportion_sig, widthfactor, 'FaceColor', colorline);
                        
                        
                        breaches50=find(proportion_sig>.5);
                        reindexed_ind=min(find(breaches50>39));
                        first50plus=breaches50(reindexed_ind);
                        if ~isempty(first50plus)
                            scatter(first50plus, signstar, '*', 'MarkerEdgeColor', colorline)
                        end
                        
                        axis([19 77.5 0 1.2])
                        
                        
                        %                         yline(total_sig_proportion,'Color',colorline)
                        xline(42, '--')
                        xline(39, '--')

                        times=(T-1)*1000;
                        xticks([19 29 39 49 59 69]);
                        xticklabels({'-500', '-250', '0', '250', '500', '750'});
                        xlabel('ms')
                        yticks([0 1]);
                        
                        counter=counter+1;
                        
                    end
                    
                    %             legend()
                    %
                    %             if b==2
                    % %                 fh=findobj( 'Type', 'Figure', 'Name', currentarea);
                    % %                 figure(fh.Number);
                    %
                    %                 legappend([p(1) p(2)], [analysistype ' Left'], [analysistype ' Right'])
                    %
                    % %                 legend([p(1) p(2)], [analysistype ' Left'], [analysistype ' Right'])
                    %             end
                end
                
            end
            
        end
        
        
        for a=1:length(plottedareas) %[1:8, 10:size(table_of_faces,2)]
            
            currentarea=char(C(plottedareas(a)));
            
            fh=findobj( 'Type', 'Figure', 'Name', currentarea);
            figure(fh.Number)
            allnames1=[];
            for b=1:size(h,2)
                
                linename={h(b).DisplayName};
                allnames1=[allnames1; linename];
                
            end
            
            all_lines_area1=contains(allnames1, currentarea);
            linenumbers=find(all_lines_area1);
            currentlinenames=allnames1(linenumbers);
            
            legend(h(linenumbers), currentlinenames)
        end
        
    case 'collapsed_hemisphere'
        for type=1:length(analysistypes)
            analysistype=char(analysistypes(type));
            
            path_1=['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_separate_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask\' analysistype '\'];
            
            path_2=findstr(path_1, '\stats\');
            analysisname_1=path_1(path_2+7:end);
            
            analysisname_2=strrep(analysisname_1, '_', ' ');
            analysisname_2=strrep(analysisname_2, '\', ', ');
            
            
            mappeddatadrive='G:\mnt\Data25\';
            
            %%grab significant vertices
            load (['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_separate_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask\' analysistype '\significantvertices_500ms_separate_baseline__baselinestats_filtered_2to2s_' num2str(parcelnumber) 'parcels.mat'])
            load (['D:\HNCT\HNCT Auditory\Analysis Code\stats_and_parcellation\stats\100msbins75msoverlap\500ms_separate_baseline__baselinestats_filtered_2to2s\' num2str(parcelnumber) '_parcels_kmeans_mask\' analysistype '\vertexsum_500ms_separate_baseline__baselinestats_filtered_2to2s_' num2str(parcelnumber) 'parcels.mat'])
            
            %%load in table of left hemisphere HOAs
            load HOArefinedtable_min200_smoothing6_lh.mat
            refinedtable_lh=refinedtable;
            
            %%load in table of right hemisphere HOAs
            load HOArefinedtable_min200_smoothing6_rh.mat
            refinedtable_rh=refinedtable;
            
            %%find all the areas that have ROIs
            [C,ia,ib] = intersect(refinedtable_lh.Properties.VariableNames,refinedtable_rh.Properties.VariableNames, 'stable');
            
            %%load mesh
            load allmni2fsfaces.mat
            
            %%load timing info
            xtick_spacing_filename = [mappeddatadrive 'HNCT\icEEG Analysis\Analysis\EEG_behavior\hnct_iceeg_scripting\eeg_analysis\filter_times_win100ms_shift25ms.mat'];
            load (xtick_spacing_filename)
            
            
            
            for a=1:length(C)
                
                currentarea=char(C(a));
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
                        currentarealocation_left=find(contains(refinedtable_lh.Properties.VariableNames, currentarea));
                        table_of_faces=refinedtable_lh;
                        allfaces=allfaces_lh;
                        allvertices=allvertices_lh;
                    else
                        side='rh';
                        side_long='right';
                        
                        currentarealocation=find(contains(refinedtable_lh.Properties.VariableNames, currentarea));
                        table_of_faces=refinedtable_rh;
                        allfaces=allfaces_rh;
                        allvertices=allvertices_rh;
                        
                    end
                    
                    currentarea_array=table2array(table_of_faces(:,a));
                    areas_faces=find(currentarea_array(:,1)~=0);
                    
                    
                    face_to_vertex1=allfaces(areas_faces,:);
                    vertices_in_region=unique(face_to_vertex1);
                    
                    if isequal(side, 'lh')
                        sign_or_no=find(sum(sign_left(vertices_in_region,:),2)~=0);
                        signverts_in_region=vertices_in_region(sign_or_no);
                        activity_insignverts_inregion_left=leftverts(signverts_in_region,:);
                        
                    else
                        sign_or_no=find(sum(sign_right(vertices_in_region,:),2)~=0);
                        signverts_in_region=vertices_in_region(sign_or_no);
                        activity_insignverts_inregion_right=rightverts(signverts_in_region,:);
                        
                    end
                end
                
                
                
                
                switch analysistype
                    case 'CP'
                        colorline=[0 0 1];
                    case 'CNP'
                        
                        colorline=[1 0 0];
                        
                    case 'Subtraction'
                        colorline=[1 0 1];
                end
                
                linenumber=counter;
                
                activity_insignverts_inregion=[activity_insignverts_inregion_left; activity_insignverts_inregion_right];
                zero_rows=sum(activity_insignverts_inregion,2)==0;
                activity_insignverts_inregion(zero_rows,:)=[];
                avactivity=nanmean(activity_insignverts_inregion,1);
                
                if ~isequal(sum(~isnan(avactivity)),0)
                    
                    stderrorofmean= std( activity_insignverts_inregion ) / sqrt( length( activity_insignverts_inregion ));
                    options.error='sem';
                    options.color_area=colorline;
                    
                    
                    plot_areaerrorbar(activity_insignverts_inregion, options);
                    h(linenumber)= plot(avactivity, 'Color', colorline, 'DisplayName', [currentarea ' ' analysistype  ' both hemispheres']);
                    
                    axis([19 77 -8 20])
                    title( strrep(currentarea, '_', ' '))
                    times=(T-1)*1000;
                    xticks([19 29 39 49 59 69]);
                    xticklabels({'-500', '-250', '0', '250', '500', '750'});
                    xlabel('ms')
                    ylabel('weighted zscore')
                    
                    
                    xline(39, '--')
                    counter=counter+1;
                    
                end
                
                
                
                
            end
            
        end
        
        
        for a=1:length(C)
            
            currentarea=char(C(a));
            
            fh=findobj( 'Type', 'Figure', 'Name', currentarea);
            figure(fh.Number)
            allnames1=[];
            for b=1:size(h,2)
                
                linename={h(b).DisplayName};
                allnames1=[allnames1; linename];
                
            end
            
            all_lines_area1=contains(allnames1, currentarea);
            linenumbers=find(all_lines_area1);
            currentlinenames=allnames1(linenumbers);
            
            legend(h(linenumbers), currentlinenames)
        end
        
end