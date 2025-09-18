if contains(filename, '548DR')
    load ('D:\HNCT\HNCT Auditory\HNCTAuditorydata\ID\icEEG\whistle\freefield\548DR_6.mat');
elseif contains(filename, '451nh')
    load ('D:\HNCT\HNCT Auditory\HNCTAuditorydata\ID\icEEG\whistle\freefield\451nh_2.mat');
end


    subjectstrials2=size(allanswers_noiserun_all, 1);
    
    
    

for sound = 1:3
    sound_rows = allanswers_noiserun_all(allanswers_noiserun_all(:, 22) == sound, :);
    sound_accuracy2(1, sound) = sum(sound_rows(:, 1) == 1 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 1);
    sound_accuracy2(2, sound) = sum(sound_rows(:, 1) == 2 & sound_rows(:, 3) ~= sound) / sum(sound_rows(:, 1) == 2);
    sound_accuracy2(3, sound) = sum(sound_rows(:, 1) == 1 & sound_rows(:, 3) ~= sound) / sum(sound_rows(:, 1) == 1);
    sound_accuracy2(4, sound) = sum(sound_rows(:, 1) == 2 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 2);
end
sound_accuracy2(:, sound_accuracy2(1,:)<.6)=NaN;

secondaryincludedsounds=find(~isnan(sound_accuracy2(1,:)));

secondarydataforincludedsounds=allanswers_noiserun_all(ismember(allanswers_noiserun_all(:,22), secondaryincludedsounds),:);

%%reload previous data
load (filename);
allanswers_noiserun_all=[allanswers_noiserun_all;secondarydataforincludedsounds];

for sound = 1:3
    sound_rows = allanswers_noiserun_all(allanswers_noiserun_all(:, 22) == sound, :);
    sound_accuracy2(1, sound) = sum(sound_rows(:, 1) == 1 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 1);
    sound_accuracy2(2, sound) = sum(sound_rows(:, 1) == 2 & sound_rows(:, 3) ~= sound) / sum(sound_rows(:, 1) == 2);
    sound_accuracy2(3, sound) = sum(sound_rows(:, 1) == 1 & sound_rows(:, 3) ~= sound) / sum(sound_rows(:, 1) == 1);
    sound_accuracy2(4, sound) = sum(sound_rows(:, 1) == 2 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 2);
end

sound_accuracy(:, sound_accuracy(1,:)<.6)=NaN;



    subjectstrials=size(allanswers_noiserun_all, 1);
    trialcount=[trialcount;subjectstrials+subjectstrials2];
    
    
    
