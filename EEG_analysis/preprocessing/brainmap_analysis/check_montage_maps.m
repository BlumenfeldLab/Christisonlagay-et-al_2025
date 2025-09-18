

% define data location
data_location = 'E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\HNCT ERP movies\';

% load pipeline info
load('E:\HNCT\icEEG Analysis\Analysis\EEG_behavior\normal_pipeline_file_info\normal_pipeline_file_info_20_04_27.mat')

disp('checking MontageMaps...')
for i = 1 : 31
    fprintf('\n%s\t(%d / %d)', patients{i}, i, 31)

    % load montage
    load([data_location, patients{i}, '\MontageMap.mat'])
    
    % check electrode locations and ensure that none are 0, 0, 0
    for j = 1 : size(MontageMap, 1)
        if MontageMap(j, 2 : 4) == 0
            fprintf('\n\t%d no good', j)
        end
    end
end

disp('checking L_MontageMaps...')
for i = 1 : 31
    fprintf('\n%s\t(%d / %d)', patients{i}, i, 31)

    % load montage
    load([data_location, patients{i}, '\L_MontageMap.mat'])
    
    % check electrode locations and ensure that none are 0, 0, 0
    for j = 1 : size(L_MontageMap, 1)
        if L_MontageMap(j, 2 : 4) == 0
            fprintf('\n\t%d no good', j)
        end
    end
end

disp('checking R_MontageMaps...')
for i = 1 : 31
    fprintf('\n%s\t(%d / %d)', patients{i}, i, 31)

    % load montage
    load([data_location, patients{i}, '\R_MontageMap.mat'])
    
    % check electrode locations and ensure that none are 0, 0, 0
    for j = 1 : size(R_MontageMap, 1)
        if R_MontageMap(j, 2 : 4) == 0
            fprintf('\n\t%d no good', j)
        end
    end
end