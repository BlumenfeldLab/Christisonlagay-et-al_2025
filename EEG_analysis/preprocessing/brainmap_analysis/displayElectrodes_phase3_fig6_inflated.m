
% Clustering
cd(image_storage_folder)

% define the range of number of clusters you want
min_clusters = 2;
max_clusters = 10;

% combine left and right side vertex time series
if laterality == 0
    vertex_values_sum_bilateral = cat(1,vertex_values_sum_L,vertex_values_sum_R);
    poststim = vertex_values_sum_bilateral(:,29:end);
elseif laterality == 3
    poststim = vertex_values_sum(:,29:57);
end

% Only want to do kmeans clustering for post stimulus ERPs
% get rid of the vertices that have very small standard deviations but
% remember them so you can accurately assign colors to the correct
% vertices later
poststim_std = std(poststim, 0, 2);
constant_vertices_idx = poststim_std == 0;
changing_vertices_idx = poststim_std > 0;
vertices_by_cluster = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_200ms = zeros(numel(changing_vertices_idx), max_clusters);
vertices_by_cluster_400ms = zeros(numel(changing_vertices_idx), max_clusters);

% take only the vertices that are not constant
poststim_reduced = poststim(changing_vertices_idx,:);
if laterality == 0
    vertex_values_sum_bilateral_reduced = vertex_values_sum_bilateral(changing_vertices_idx,:);
    vertex_values_sum_reduced = vertex_values_sum_bilateral_reduced;
else
    vertex_values_sum_reduced = vertex_values_sum(changing_vertices_idx,:);
end

% redo the above for the first 200 ms (4 frames) post-stimulus
poststim_reduced_200ms = poststim_reduced(:,1:7);
poststim_reduced_400ms = poststim_reduced(:,1:13);

% Intialize vertices x clusters NaN matrices...
vertices_by_cluster_reduced = NaN(sum(changing_vertices_idx),max_clusters);
vertices_by_cluster_reduced_200ms = NaN(sum(changing_vertices_idx),max_clusters);
vertices_by_cluster_reduced_400ms = NaN(sum(changing_vertices_idx),max_clusters);

% Initialize vars for sillouette values, average sillouette values and p
% values per cluster
s = NaN(sum(changing_vertices_idx), max_clusters);
s_mean = NaN(max_clusters);
p_val = NaN(max_clusters);
s_200ms = NaN(sum(changing_vertices_idx), max_clusters);
s_200ms_mean = NaN(max_clusters);
p_val_200ms = NaN(max_clusters);
s_400ms = NaN(sum(changing_vertices_idx),max_clusters);
s_400ms_mean = NaN(max_clusters);
p_val_400ms = NaN(max_clusters);

% Initialize vars for something...
meantrace = NaN(max_clusters, max_clusters, 57);
baseline_trace = NaN(max_clusters, 28);
poststim_trace = NaN(max_clusters, 29);
meantrace_200ms = NaN(max_clusters, max_clusters, 57);
baseline_trace_200ms = NaN(max_clusters, 28);
poststim_trace_200ms = NaN(max_clusters, 29);
meantrace_400ms = NaN(max_clusters, max_clusters, 57);
baseline_trace_400ms = NaN(max_clusters, 28);
poststim_trace_400ms = NaN(max_clusters, 29);

load('spectrogram_times_57bins.mat')
for i = min_clusters : max_clusters
    
    % Indicate status
    fprintf('\Computing %d clusters\n', i)
    
    % kmeans clustering of poststimulus activity
    [vertices_by_cluster_reduced(:,i), ~]= kmeans(poststim_reduced, i, 'dist', 'correlation', 'rep',...
        20, 'disp', 'final');
    
    % Graph silhouette values
    disp(['Calculating silhouette values for ' num2str(i) ' clusters'])
    figure;
    cidx = vertices_by_cluster_reduced(:, i);
    [s(:, i), ~] = silhouette(poststim_reduced, cidx, 'correlation');
    
    % calculate the mean of the silhouette values in that cluster
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced(:, i) == j);
        s_mean(i, j) = mean(s(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_' num2str(i) '.tif']);
    
    % Graph cluster time courses
    figure
    legend_str = cell(1, i);
    % Gamma part earlier in the code because of graphing reasons
    for j = 1:i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced(:,i) == j);
        vals = vertex_values_sum_reduced(cluster_vertex_idx, :);
        plot(T-1,nanmean(vals,1), 'Color', ColorSet(j,:),'Linewidth',2); hold on;
        meantrace(i,j,:) = nanmean(vals,1);
        baseline_trace(j,:) = abs(meantrace(i,j,1:28));
        poststim_trace(j,:) = abs(meantrace(i,j,29:57));
        [p_val(i,j)] = ranksum(baseline_trace(j,:),poststim_trace(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val(i,j)) ', s = ' num2str(s_mean(i,j))];
    end
    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power z-score','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim([-.5 1])
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest')
    legend('boxoff');
    set(gca,'Fontsize',14)
    set(gcf, 'color', [1 1 1]);
    saveas(gcf,['Gamma_clusters_' num2str(i) '.tif']);
    
    clearvars legend_str
    
    %% K means clustering of the first 200 ms (4 frames) post-stimulus
    
    [vertices_by_cluster_reduced_200ms(:,i), ~]= kmeans(poststim_reduced_200ms, i, 'dist', 'correlation',...
        'rep', 20, 'disp', 'final');
    
    % Graph clusters
    
    
    % Graph silhouette values
    disp(['Calculating silhouette values for ' num2str(i) ' clusters'])
    figure
    cidx = vertices_by_cluster_reduced_200ms(:, i);
    [s_200ms(:, i), ~] = silhouette(poststim_reduced, cidx,'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_200ms(:, i) == j);
        s_200ms_mean(i, j) = mean(s_200ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_200ms_' num2str(i) '.tif']);
    
    % Gamma part earlier in the code because of graphing reasons
    for j = 1:i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_200ms(:,i) == j);
        % plot that index
        plot(T-1,nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1), 'Color', ColorSet(j,:)); hold on;
        
        % do stats on the prestim and poststim gamma power zscore changes
        meantrace_200ms(i,j,:) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1);
        baseline_trace_200ms(j,:) = abs(meantrace_200ms(i,j,1:28));
        poststim_trace_200ms(j,:) = abs(meantrace_200ms(i,j,29:57));
        [p_val_200ms(i,j)] = ranksum(baseline_trace_200ms(j,:),poststim_trace_200ms(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_200ms(i,j)) ', s = ' num2str(s_200ms_mean(i,j))];
    end
    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim([-.5 1])
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest','AutoUpdate', 'off')
    legend('boxoff');
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T-1,squeeze(meantrace(i,j,:)), 'Color',[.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_200ms_' num2str(i) '.tif']);
    
    clearvars legend_str
    
    % K means clustering of the first 400 ms (4 frames) post-stimulus
    
    [vertices_by_cluster_reduced_400ms(:,i), centroids]= kmeans(poststim_reduced_400ms, i, 'dist',...
        'correlation', 'rep', 20, 'disp', 'final');
    
    
    
    
    % Graph silhouette values
    figure
    disp(['Calculating silhouette values for ' num2str(i) ' clusters'])
    cidx = vertices_by_cluster_reduced_400ms(:, i);
    [s_400ms(:, i), ~] = silhouette(poststim_reduced,cidx, 'correlation');
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_200ms(:, i) == j);
        s_400ms_mean(i, j) = mean(s_400ms(cluster_vertex_idx, i));
    end
    title(['Silhouette values for ' num2str(i) ' clusters'])
    saveas(gcf,['Silhouettes_400ms_' num2str(i) '.tif']);
    
    % Gamma part earlier in the code because of graphing reasons
    for j = 1 : i
        % find the indices of vertices both left and right that are part of
        % that cluster
        cluster_vertex_idx = find(vertices_by_cluster_reduced_400ms(:, i) == j);
        % plot that index
        plot(T-1,nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1), 'Color', ColorSet(j,:)); hold on;
        % do stats on the prestim and poststim gamma power zscore changes
        meantrace_400ms(i,j,:) = nanmean(vertex_values_sum_reduced(cluster_vertex_idx,:),1);
        baseline_trace_400ms(j,:) = abs(meantrace_400ms(i,j,1:28));
        poststim_trace_400ms(j,:) = abs(meantrace_400ms(i,j,29:57));
        [p_val_400ms(i,j)] = ranksum(baseline_trace_400ms(j,:),poststim_trace_400ms(j,:));
        legend_str{j} = ['Cluster ' num2str(j) ', p = ' num2str(p_val_400ms(i,j)) ', s = ' num2str(s_400ms_mean(i,j))];
    end
    
    title(['Gamma Power Zscore, # Clusters = ' num2str(i) ],'Fontsize',14);
    xlabel('Time (s)','FontSize',14)
    ylabel('Power','FontSize',14)
    xlim([T(1)-1 T(end)-1]);
    ylim([-.5 1])
    xL = xlim;
    yL = ylim;
    line([0 0], yL,'Color','k');
    line(xL, [0 0],'Color','k');
    legend(legend_str,'Location','Northwest', 'AutoUpdate', 'off')
    legend('boxoff');
    % set the size of the legend and axes fonts
    set(gca,'Fontsize',14)
    % insert shadows to the lines
    for j = 1:i
        shadow_line(j) = plot(T-1,squeeze(meantrace(i,j,:)), 'Color',[.8 .8 .8],'Linewidth',4); hold on;
        uistack(shadow_line(j), 'bottom')
    end
    saveas(gcf,['Gamma_clusters_400ms_' num2str(i) '.tif']);
    
    clearvars legend_str
end


%% put back the constant vertices so they are in the same order again
vertices_by_cluster(changing_vertices_idx, 1 : size(vertices_by_cluster_reduced, 2)) = vertices_by_cluster_reduced;
vertices_by_cluster_200ms(changing_vertices_idx,:) = vertices_by_cluster_reduced_200ms;
vertices_by_cluster_400ms(changing_vertices_idx,:) = vertices_by_cluster_reduced_400ms;

% separate into left and right vertices again (if you aren't doing both
% sides on the left)
if laterality == 0
    vertices_by_cluster_L = vertices_by_cluster(1:size(vertex_values_sum_L, 1),:);
    vertices_by_cluster_R = vertices_by_cluster(size(vertex_values_sum_L, 1)+1:end,:);
    
    vertices_by_cluster_L_200ms = vertices_by_cluster_200ms(1:size(vertex_values_sum_L, 1),:);
    vertices_by_cluster_R_200ms = vertices_by_cluster_200ms(size(vertex_values_sum_L, 1)+1:end,:);
    
    vertices_by_cluster_L_400ms = vertices_by_cluster_400ms(1:size(vertex_values_sum_L, 1),:);
    vertices_by_cluster_R_400ms = vertices_by_cluster_400ms(size(vertex_values_sum_L, 1)+1:end,:);
end

%%graph clusters onto brain maps

for side_index = start_side : sides_num
    
    load ('cmap6to8WRX.mat');
    
    % indicate the figure of the left or right side of the brain
    figure(figHandle(side_index));
    
    % NOTE: now the frame index is actually the number of clusters in each
    % frame rather a frame in terms of a time point
    for frame_index = min_clusters : max_clusters
        
        % overlay the colors
        disp(['>> Rendering ' side ' Cortex, Frame ' num2str(frame_index)]);  % Output Status
        vertex_values_sum = [];
        if laterality == 3
            vertex_values_sum = vertices_by_cluster(:,frame_index);
        elseif side_index == 1 && laterality == 0
            vertex_values_sum = vertices_by_cluster_L(:,frame_index);
        elseif side_index == 2 && laterality == 0
            vertex_values_sum = vertices_by_cluster_R(:,frame_index);
        end
        colormap(ColorSetMaps)
        caxis([-0.5 10.5])
        
        if any(electrodes_present(side_index,:) == 1)
            % for now just overlay the color as number of patients overlapping in
            % this area without any transparency
            h = patch('Faces',inflated_surf(side_index).faces, 'Vertices', inflated_surf(side_index).vertices,...
                'FaceVertexCData', vertex_values_sum, 'FaceColor', 'interp',...
                'LineStyle', 'none', 'FaceAlpha', alpha);
            lightset = [0.6 0.5 0.1];
            material(lightset);
        end
        
        % Save each side of the brain in 4 views
        frames_set(side_index,frame_index,:) = savebrainimages_inflated(side_index, frame_index, views);
        if any(electrodes_present(side_index, :) == 1) && exist('h', 'var')
            % get rid of the color patch so you can lay on the color patch for the
            % next frame
            delete(h)
        end
        
    end
end

% Put the 2 or 4 views into one image
for frame_index = min_clusters : max_clusters
    
    ViewAngles={'Lateral','Medial','Ventral','Posterior'};
    for side_index = start_side:sides_num
        if side_index == 1
            side = 'L';
        else
            side = 'R';
        end
        for view_index = 1:numel(ViewAngles)
            eval([side '_' ViewAngles{view_index} ' = [];'])
        end
    end
    
    % name the frames for each side
    for side_index = start_side:sides_num
        % Determine what side we are working on
        if(side_index == 1)
            side = 'L';
            ViewAngles = {'Medial'};
        else
            side = 'R';
            ViewAngles = {'Lateral'};
        end
        
        if views == 4
            ViewAngles={'Lateral','Medial','Ventral','Posterior'};
        elseif views == 2
            ViewAngles={'Lateral','Medial'};
        elseif strcmp(views, 'Lateral')
            ViewAngles={'Lateral'};
        elseif strcmp(views, 'Medial')
            ViewAngles={'Medial'};
        end
        
        for view_index=1:length(ViewAngles)
            Frame = frames_set(side_index,frame_index,view_index).cdata;
            % Do cropping here
            Frame(all(all(Frame == 255, 3), 2), :, :) = [];
            Frame(:, all(all(Frame == 255, 3), 1), :) = [];
            eval([side '_' ViewAngles{view_index} ' = Frame;'])
        end
    end
    
    
    %% create the frames for each side
    for side_index = start_side:sides_num
        if side_index == 1
            
            L_Lateral_Padded = [L_Lateral, ones(size(L_Lateral, 1),...
                max([0, size(L_Medial, 2) - size(L_Lateral, 2)]), 3) * 255];
            
            L_Medial_Padded = [L_Medial, ones(size(L_Medial, 1),...
                max([0, size(L_Lateral, 2) - size(L_Medial, 2)]), 3) * 255];
            
            if laterality ~= 3
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(R_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];
            else
                L_Ventral_Padded = [L_Ventral, ones(size(L_Ventral, 1),...
                    max([0, size(L_Posterior, 2) - size(L_Ventral, 2)]), 3) * 255];
            end
            
            L_Lateral_Medial = [L_Medial_Padded; L_Lateral_Padded];
            
            
            if laterality == 0
                R_Posterior_Padded = [R_Posterior, ones(size(R_Posterior, 1),...
                    max([0, size(L_Ventral, 2) - size(R_Posterior, 2)]), 3) * 255];
                
                L_Ventral_R_Posterior = [L_Ventral_Padded; R_Posterior_Padded];
                
                L_Lateral_Medial_Padded = [L_Lateral_Medial; ones(max([0, size(L_Ventral_R_Posterior, 1) - ...
                    size(L_Lateral_Medial, 1)]), size(L_Lateral_Medial, 2), 3) * 255];
                
                L_Ventral_R_Posterior_Padded = [L_Ventral_R_Posterior; ones(max([0, size(L_Lateral_Medial, 1) - ...
                    size(L_Ventral_R_Posterior, 1)]), size(L_Ventral_R_Posterior, 2), 3) * 255];
                
                combined_views_left{frame_index} = [L_Ventral_R_Posterior_Padded, L_Lateral_Medial_Padded];
            elseif laterality == 3
                L_Posterior_Padded = L_Posterior;
                
                L_Ventral_L_Posterior = [L_Ventral_Padded;L_Posterior_Padded];
                
                L_Lateral_Medial_Padded = [L_Lateral_Medial; ones(max([0, size(L_Ventral_L_Posterior, 1) - ...
                    size(L_Lateral_Medial, 1)]), size(L_Lateral_Medial, 2), 3) * 255];
                
                L_Ventral_L_Posterior_Padded = [L_Ventral_L_Posterior; ones(max([0, size(L_Lateral_Medial_Padded, 1) - ...
                    size(L_Ventral_L_Posterior, 1)]), size(L_Ventral_L_Posterior, 2), 3) * 255];
                
                combined_views_left{frame_index} = [L_Ventral_L_Posterior_Padded, L_Lateral_Medial_Padded];
            end
            
            figure
            imshow(combined_views_left{frame_index});
            set(gca,'position',[0 0 1 1],'units','normalized')
            
            % adding time stamp text to the figure
            if laterality == 3 && ismember(overlay,[2 3])
                load(xtick_spacing_filename)
                time_txt = ['Time = ' num2str(round((T(frame_index)-1)*1000)) 'ms'];
                text(400,900,time_txt,'Fontsize',24)
            end
            fn = ['combined_viewsL' '_' num2str(frame_index)];
            print(fn, '-dtiff');
            close
            
        elseif side_index == 2
            
            R_Lateral_Padded = [R_Lateral, ones(size(R_Lateral, 1),...
                max([0, size(R_Medial, 2) - size(R_Lateral, 2)]), 3) * 255];
            
            R_Medial_Padded = [R_Medial, ones(size(R_Medial, 1),...
                max([0, size(R_Lateral, 2) - size(R_Medial, 2)]), 3) * 255];
            
            R_Lateral_Medial = [R_Medial_Padded; R_Lateral_Padded];
            
            
            R_Ventral_Padded = [R_Ventral, ones(size(R_Ventral, 1),...
                max([0, size(L_Posterior, 2) - size(R_Ventral, 2)]), 3) * 255];
            
            if laterality == 0
                L_Posterior_Padded = [L_Posterior, ones(size(L_Posterior, 1),...
                    max([0, size(R_Ventral, 2) - size(L_Posterior, 2)]), 3) * 255];
                
                R_Ventral_L_Posterior = [R_Ventral_Padded; L_Posterior_Padded];
                
                R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(max([0, size(R_Ventral_L_Posterior, 1) - ...
                    size(R_Lateral_Medial, 1)]), size(R_Lateral_Medial, 2), 3) * 255];
                
                R_Ventral_L_Posterior_Padded = [R_Ventral_L_Posterior; ones(max([0, size(R_Lateral_Medial, 1) - ...
                    size(R_Ventral_L_Posterior, 1)]), size(R_Ventral_L_Posterior, 2), 3) * 255];
                
                combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_L_Posterior_Padded];
            elseif laterality == 2
                R_Posterior_Padded = R_Posterior; %(80:495,155:395,:);
                
                R_Lateral_Medial_Padded = [R_Lateral_Medial; ones(116,545,3)*255];
                
                R_Ventral_R_Posterior = [R_Ventral_Padded;R_Posterior_Padded];
                
                combined_views_right{frame_index} = [R_Lateral_Medial_Padded, R_Ventral_R_Posterior];
                combined_views_right{frame_index} = [combined_views_right{frame_index}; 255*ones(15, 776, 3)];
            end
            
            figure
            imshow(combined_views_right{frame_index});
            set(gca,'position',[0 0 1 1],'units','normalized')
            fn = ['combined_viewsR' '_' num2str(frame_index)];
            saveas(gcf,[fn '.tif'])
            close
        end
    end
    
    if laterality == 0
        figure
        combined_views = [[combined_views_right{frame_index}; ones(max([0, size(combined_views_left{frame_index}, 1) - ...
            size(combined_views_right{frame_index}, 1)]), size(combined_views_right{frame_index}, 2), 3) * 255]...
            [combined_views_left{frame_index}; ones(max([0, size(combined_views_right{frame_index}, 1) - ...
            size(combined_views_left{frame_index}, 1)]), size(combined_views_left{frame_index}, 2), 3) * 255]];
        imshow(combined_views);
        set(gca,'position',[0 0 1 1],'units','normalized')
        axis tight
        frame(frame_index,:,:,:) = combined_views;
        if exist('j','var')
            fn = ['combined_views_full_' num2str(frame_index) '_' num2str(j)];
        else
            fn = ['combined_views_full_' num2str(frame_index)];
        end
        
        if ismember(overlay,[2 3])
            load(xtick_spacing_filename)
            time_txt = ['Time = ' num2str(round((T(frame_index)-1)*1000)) 'ms'];
            if views == 4
                uc = uicontrol('Style', 'text', 'Visible', 'off', 'FontName', 'Helvetica',...
                    'FontSize', 48, 'String', time_txt);
                text(size(combined_views, 2) - 1.5 * uc.Extent(3),...
                    size(combined_views, 1) - uc.Extent(4) - 20, time_txt, 'Fontsize', 48)
            end
        elseif overlay == 1
            T = [100 200 300 400 500 600 700 800 900 1000];
            time_txt = ['Time = ' num2str(T(frame_index)) ' ms'];
            text(txt_pos_x,txt_pos_y,time_txt,'Fontsize',20)
        end
        
        saveas(gcf,[fn '.tiff'])
    end
end
