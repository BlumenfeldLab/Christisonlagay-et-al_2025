% Rejects sounds with Confirmed Perceived accuracy < 60%
rootfolder = '';

folders = {} %epoch folders

b_file = ''; %behavioral folder

load(b_file)

sound_accuracy = NaN(2, 3);

for sound = 1:3
    sound_rows = allanswers_noiserun_all(allanswers_noiserun_all(:, 22) == sound, :);
    sound_accuracy(1, sound) = sum(sound_rows(:, 1) == 1 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 1);
    sound_accuracy(2, sound) = sum(sound_rows(:, 1) == 2 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 2);
end

% define trial types also in that order
trialtypes = {'BC_trials','WG_trials','blank_TN_trials'};

%% Aggregate Ordinary Epochs

BC_trials = [];
WG_trials = [];
blank_TN_trials = [];
% CG_trials = [];
% FP_trials = [];
% blank_FP_trials = [];

for i = 1:length(folders)
    disp(i);
    folder = folders{i};
    cd(folder);

    % NOTE: More sessions than recording files, using run numbers in
    % containing folders to extract appropriate rows
    folder_split = split(folders{i}, filesep);
    run_bounds = str2double(regexp(folder_split{end - 1}, '\d+', 'match'));
    if numel(run_bounds) == 1
        run_range = run_bounds;
    else
        run_range = run_bounds(1):run_bounds(2);
    end

    missing_trials = [];
    if exist('missing_trials.mat', 'file')
        load('.\missing_trials.mat')
    end
    load('Tone_playback_baseline500to0.mat');
    load('trial_identifiers.mat');

    session_rows = allanswers_noiserun_all(any(allanswers_noiserun_all(:, 13) == run_range - 1, 2), :);

    count = 0;

    delay_identifier(missing_trials) = [];
    trial_identifier(missing_trials) = [];
    type_identifier(missing_trials) = [];
    session_rows(missing_trials, :) = [];

    % go through each row in trial identifiers in order to sort epochs
    for k = 1:numel(trial_identifier)
        count = count + 1;
        sound = session_rows(count, 22);
        if sound && (sound_accuracy(1, sound) < 0.6)
            continue
        end
        % isolate that epoch
        trial_epoch =  epochs(:,:,k);
        if strcmp(trial_identifier{k,1},'Confirmed perceived') == 1
            j = 1;
        elseif strcmp(trial_identifier{k,1},'Confirmed not perceived') == 1
            j = 2;
        elseif strcmp(trial_identifier{k,1},'True Negative') == 1
            j = 3;
        else
            j = NaN;
        end

        if any(j)
            eval([trialtypes{j} ' = cat(3,trial_epoch,' trialtypes{j} ');']);
        end
    end

    clearvars epochs trial_identifiers

end

cd(rootfolder)
save('sorted_trials_sounds_restricted_CPonly_baseline500to0','BC_trials','WG_trials','blank_TN_trials','-v7.3');

%{
%% aggregate extended trials
BC_trials = [];
WG_trials = [];
blank_TN_trials = [];

for i = 1:length(folders)
    disp(i);
    folder = folders{i};
    cd(folder);

    % NOTE: More sessions than recording files, using run numbers in
    % containing folders to extract appropriate rows
    folder_split = split(folders{i}, filesep);
    run_bounds = str2double(regexp(folder_split{end - 1}, '\d+', 'match'));
    if numel(run_bounds) == 1
        run_range = run_bounds;
    else
        run_range = run_bounds(1):run_bounds(2);
    end

    missing_trials = [];
    if exist('missing_trials.mat', 'file')
        load('.\missing_trials.mat')
    end
    load('extended_tone_epochs.mat');
    load('trial_identifiers.mat');

    session_rows = allanswers_noiserun_all(any(allanswers_noiserun_all(:, 13) == run_range - 1, 2), :);

    count = 0;

    delay_identifier(missing_trials) = [];
    trial_identifier(missing_trials) = [];
    type_identifier(missing_trials) = [];
    session_rows(missing_trials, :) = [];

    % go through each row in trial identifiers in order to sort epochs
    for k = 1:numel(trial_identifier)
        count = count + 1;
        sound = session_rows(count, 22);
        if sound && (sound_accuracy(1, sound) < 0.6)
            continue
        end
        % isolate that epoch
        trial_epoch =  epochs(:,:,k);
        if strcmp(trial_identifier{k,1},'Confirmed perceived') == 1
            j = 1;
        elseif strcmp(trial_identifier{k,1},'Confirmed not perceived') == 1
            j = 2;
        elseif strcmp(trial_identifier{k,1},'True Negative') == 1
            j = 3;
        else
            j = NaN;
        end

        if any(j)
            eval([trialtypes{j} ' = cat(3,trial_epoch,' trialtypes{j} ');']);
        end
    end

    clearvars epochs trial_identifiers

end

cd(rootfolder)
save('sorted_trials_sounds_restricted_CPonly_extended','BC_trials','WG_trials','blank_TN_trials','-v7.3');


%% aggregate by hemifield

% define trials if you want to do left vs. right analysis
trialtypes_hemifield = {'BC_trials_left','WG_trials_left','blank_TN_trials_left',...
    'BC_trials_right','WG_trials_right','blank_TN_trials_right'};

trialtypenames_hemifield = {'Confirmed Perceived Left','Confirmed Not Perceived Left',...
    'True Negatives Left','Confirmed Perceived Right','Confirmed Not Perceived Right',...
    'True Negatives Right'};

BC_trials_left = [];
WG_trials_left = [];
blank_TN_trials_left = [];
BC_trials_right = [];
WG_trials_right = [];
blank_TN_trials_right = [];

for i = 1:length(folders)
    disp(i);
    folder = folders{i};
    cd(folder);
    missing_trials = [];
    if exist('missing_trials.mat', 'file')
        load('.\missing_trials.mat')
    end
    load('Tone_playback.mat');
    load('trial_identifiers.mat');
    delay_identifier(missing_trials) = [];
    trial_identifier(missing_trials) = [];
    type_identifier(missing_trials) = [];

    % go through each row in trial identifiers in order to sort epochs
    for k = 1:numel(trial_identifier)
        % isolate that epoch
        trial_epoch =  epochs(:,:,k);

        if strcmp(hemifield_identifier{k,1},'Left') == 1
            j = 0;
        elseif strcmp(hemifield_identifier{k,1},'Right') == 1
            j = 3;
        end


        if strcmp(trial_identifier{k,1},'Confirmed perceived') == 1
            j = j + 1;
        elseif strcmp(trial_identifier{k,1},'Confirmed not perceived') == 1
            j = j + 2;
        elseif strcmp(trial_identifier{k,1},'True Negative') == 1
            j = j + 3;
        else
            j = NaN;
        end

        if any(j)
            eval([trialtypes_hemifield{j} '= cat(3,trial_epoch,' trialtypes_hemifield{j} ');']);
        end
    end

    clearvars epochs trial_identifiers

end

cd(rootfolder)
save('sorted_trials_hemifield','BC_trials_left','WG_trials_left','blank_TN_trials_left',...
    'BC_trials_right','WG_trials_right','blank_TN_trials_right','-v7.3');
%}
