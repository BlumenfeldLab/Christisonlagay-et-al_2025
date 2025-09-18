%==================================================================
%
%Title:         electrode_data_overlay
%Description:   This function loads a surface txt file output from Bioimage
%               Suite into Matlab, and arranges it into a readable format.
%               Repeat twice (once for cortex and once for smooth). The
%               output is a reorganized matrix of vertices.
%
%Inputs:        "vertices_num"  [int] The number of vertices in the file
%               "faces_num"     [int] The number of faces in the file
%               "data_location" [string] should point to the root folder
%                               for the seizure
%                               i.e. 'C:\IcEEGProject_ProcessedData\25\25_1'
%               "filename"      [string] Location of the vertice data
%==================================================================

function [electrode_vertex_values] = electrode_data_overlay(overlay, frame_index, all_electrodes, subtraction, file_suffix)
overlay_frame_data = [];
if overlay == 0 || overlay == 1
    %% Give every overlay value a 1
    % once aggregated across patients it gives the density map of
    % electrodes at each vertex
    %     load('labels.mat');
    try
        load('labels_all_gray_matter.mat');
    catch
        load('labels_first_surgery_all_gray_matter.mat');
    end
    overlay_frame_data = ones(length(labels),1);
    % elseif overlay == 1
    %     %% First arrival times - need to add code for 2nd and 3rd------------------
    %     load('arrival_time.mat')%, 'half_max_times');
    %     newTimes = [];
    %     for i = 1:size(half_max_times,2)
    %         arrivalOne = half_max_times{1,i};
    %         if ~isempty(half_max_times{1,i})
    %             firstTime = arrivalOne(1);
    %         else
    %             firstTime = 0;
    %         end
    %         newTimes = [newTimes, firstTime];
    %     end
    %     newTimes = newTimes';
    %     overlay_frame_data = newTimes;
elseif overlay == 2 || overlay == 3
    %% getting power zscores from mean power of each electrode
    load(['meanpower_traces_57bins_zscore_recentered_rejoutliers', file_suffix, '.mat'])
    if overlay == 2
        trialtype = 1; % confirmed perceived
    elseif overlay == 3
        trialtype = 2; % confirmed not perceived
    end
    if ~subtraction
        overlay_frame_data = squeeze(meanpower_traces(:,trialtype,3,:));
    else
        overlay_frame_data = squeeze(meanpower_traces(:,1,3,:)) - squeeze(meanpower_traces(:,2,3,:));
    end
end

% once you have the overlay data, you need to match it to an electrode from
% labels

overlay_frame_data = overlay_frame_data';
electrode_vertex_values = overlay_frame_data(frame_index, all_electrodes);

end
