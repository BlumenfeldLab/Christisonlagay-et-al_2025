function [] = electrodeOrganize(patient_folder,input_name,output_name)
cd(patient_folder)
[~, ~, raw] = xlsread([input_name '.csv']);
electrodes = raw(3:end, :);
j = 1; k = 1;
mode = 1; % 1 : looking for electrode labels. 0 : looking for coordinates
for i = 1:length(electrodes)
    if ischar(electrodes{i,1}) && mode
        newElectrodes{k,1} = electrodes{i,1};
        k = k + 1;
        mode = 0;
    end
    if isnumeric(electrodes{i,1}) == 1 && isnan(electrodes{i,1}) == 0 && ~mode
        num = num2str(electrodes{i,1});
        delim = ';';
        newElectrodes{j,1} = strcat(newElectrodes{j,1}, delim, num);
        newElectrodes{j,2} = electrodes{i,2}; % x coordinates  
        newElectrodes{j,3} = electrodes{i,3}; % y coordinates 
        newElectrodes{j,4} = electrodes{i,4}; % z coordinates 
        j = j+1;
        mode = 1;
    end
end
fn = output_name; 
xlswrite(fn, newElectrodes)
end
