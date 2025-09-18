% clear all;
clearvars -except validation_*

close all

ID_icEEGsubject_database;

%% extract accuracy data
sound_accuracy_bysound_all=[]; sound_accuracy_total_all=[];
overallsound_accuracy_all=[]; validation_all=[];
timereport=[]; stimreport=[]; trialcount=[];



%%%we need to make 2 special cases b/c 2 of our people were recalibrated
%%%part way through their behavior, and so should be considered one person,
%%%but may need to have different sounds excluded

for subj=1:length(behavior_database)
    
    filename=char(behavior_database(subj));
    load (filename);
    
    %     subjnumber1=regexp(filename,'\d*','Match')
    subjnumber=filename(max(strfind(filename, '\'))+1:max(strfind(filename, '_'))-1);
    nameforfield=['Subj' subjnumber];
    
    
    
    for sound = 1:3  %%CP, CNP, 
        sound_rows = allanswers_noiserun_all(allanswers_noiserun_all(:, 22) == sound, :);
        sound_accuracy(1, sound) = sum(sound_rows(:, 1) == 1 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 1);
        sound_accuracy(2, sound) = sum(sound_rows(:, 1) == 2 & sound_rows(:, 3) ~= sound) / sum(sound_rows(:, 1) == 2);
        sound_accuracy(3, sound) = sum(sound_rows(:, 1) == 1 & sound_rows(:, 3) ~= sound) / sum(sound_rows(:, 1) == 1);
        sound_accuracy(4, sound) = sum(sound_rows(:, 1) == 2 & sound_rows(:, 3) == sound) / sum(sound_rows(:, 1) == 2);
        
    end
    
    
    %%don't include if accuracy is less than .6
    sound_accuracy(:, sound_accuracy(1,:)<.60)=NaN;
    
    if contains(filename, '548DR')|contains(filename, '451nh')
        secondarybehavioral_file_added
    else
        subjectstrials=size(allanswers_noiserun_all, 1);
        trialcount=[trialcount;subjectstrials];
    end
    
    sound_accuracy_total=nanmean(sound_accuracy,2);
    
    %each sounds' accuracy, for each person.
    sound_accuracy_bysound_all=cat(3, sound_accuracy_bysound_all,  sound_accuracy);
    
    %average across all sounds for a given person, but this doesn't account
    %that there are slightly different # of each sound, so this is
    %essentially meaningless :)
    sound_accuracy_total_all=cat(3, sound_accuracy_total_all,  sound_accuracy_total);
    
    
    
    
    overall_accuracy=[];
    whistle1=NaN;     laser1=NaN;    water1=NaN; 


    if gt(sum(~isnan(sound_accuracy(1,:))),1)
        nameforfield
        includedsounds=find(~isnan(sound_accuracy(1,:)))
        pause
        
        
        runidentification=unique(allanswers_noiserun_all(:,13));
    firsthalf_all=[];
    secondhalf_all=[];
    for runID=1:length(runidentification)
        currentrun=runidentification(runID);
        trialsinrun=find(allanswers_noiserun_all(:,13)==currentrun);
        sizetrialsinrun=length(trialsinrun);
        firsthalf=trialsinrun(1:(ceil(sizetrialsinrun/2)));
        secondhalf=trialsinrun(((ceil(sizetrialsinrun/2)+1):end));
        firsthalf_all=[firsthalf_all; firsthalf];
        secondhalf_all=[secondhalf_all; secondhalf];
    end
        
        dataforincludedsounds=allanswers_noiserun_all(ismember(allanswers_noiserun_all(:,22), includedsounds),:);
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %%%%%%%%%%question1%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

        dataforblanks=allanswers_noiserun_all((allanswers_noiserun_all(:,22)==0),:);
        
        overallsound_accuracy = sum(dataforincludedsounds(:, 1) == 1) / sum(dataforincludedsounds(:, 1) == 1 | dataforincludedsounds(:, 1) == 2);
        overallblank_accuracy = sum(dataforblanks(:, 1) == 1) / sum(dataforblanks(:, 1) == 1 | dataforblanks(:, 1) == 2);
        overall_accuracy=[overallsound_accuracy  overallblank_accuracy];

        
        %%%%question1 timing:
        prestim=dataforincludedsounds(:,7)./dataforincludedsounds(:,8)+dataforincludedsounds(:,11);
        possibledurations=unique(round(dataforincludedsounds(:,10)-prestim));
        durations=round(dataforincludedsounds(:,10)-prestim);
        
        threesecs=dataforincludedsounds([durations==3], :);
        foursecs=dataforincludedsounds([durations==4], :);
        fivesecs=dataforincludedsounds([durations==5], :);

        timesound_percept3 = sum(threesecs(:, 1) == 1) / sum(threesecs(:, 1) == 1 | threesecs(:, 1) == 2);
        timesound_percept4 = sum(foursecs(:, 1) == 1) / sum(foursecs(:, 1) == 1 | foursecs(:, 1) == 2);
        timesound_percept5 = sum(fivesecs(:, 1) == 1) / sum(fivesecs(:, 1) == 1 | fivesecs(:, 1) == 2);

        timereport=[timereport; timesound_percept3 timesound_percept4 timesound_percept5];
        
        %%%%question1 sound:
        whistle1=dataforincludedsounds([dataforincludedsounds(:,22)==1], :);
        laser1=dataforincludedsounds([dataforincludedsounds(:,22)==2], :);
        water1=dataforincludedsounds([dataforincludedsounds(:,22)==3], :);


        whistle_percept = sum(whistle1(:, 1) == 1) / sum(whistle1(:, 1) == 1 | whistle1(:, 1) == 2);
        laser_percept = sum(laser1(:, 1) == 1) / sum(laser1(:, 1) == 1 | laser1(:, 1) == 2);
        water_percept = sum(water1(:, 1) == 1) / sum(water1(:, 1) == 1 | water1(:, 1) == 2);
        stimreport=[stimreport; whistle_percept laser_percept water_percept];

        
        
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       %%%%%%%%%%question2%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       overall_confirmedperceived = sum(dataforincludedsounds(:, 1) == 1 & dataforincludedsounds(:,3)==dataforincludedsounds(:,22)) / sum(dataforincludedsounds(:, 1) == 1);
       overall_confirmednotperceived = sum(dataforincludedsounds(:, 1) == 2 & dataforincludedsounds(:,3)~=dataforincludedsounds(:,22)) / sum(dataforincludedsounds(:, 1) == 2);

       cpcnp_accuracy=[overall_confirmedperceived 1-overall_confirmednotperceived];
        
    
        
        
        
        
    
    end
    
    validation_all=[validation_all; cpcnp_accuracy];
    overallsound_accuracy_all=[overallsound_accuracy_all; overall_accuracy];

    
    
    
end

%%summary stats
meancompletedtrials=nanmean(trialcount);
medcompletedtrials=median(trialcount);
mincompleted=min(trialcount);
maxcompleted=max(trialcount);

%%question 1
overall_accuracy=nanmean(overallsound_accuracy_all);
accuracy_SEM=nanstd(overallsound_accuracy_all)/(length(overallsound_accuracy_all)^(1/2));


timingmean=nanmean(timereport);
timingSEM=nanstd(timereport)/(length(timereport)^(1/2));
timingstats=anova1(timereport);


stimmean=nanmean(stimreport);
stimSTEM=nanstd(stimreport)/(length(stimreport)^(1/2));
stimstats=anova1(stimreport);

%%question 2

validation_mean=nanmean(validation_all);
validation_mean_1=[validation_mean; [1-validation_mean]]';
validation_SEM=nanstd(validation_all)/(length(validation_all)^(1/2));


whistlestuffCP=squeeze(sound_accuracy_bysound_all(1,1,:));
laserstuffCP=squeeze(sound_accuracy_bysound_all(1,2,:));
waterdropstuffCP=squeeze(sound_accuracy_bysound_all(1,3,:));


whistlemeanCP=nanmean(whistlestuffCP);
waterdropmeanCP=nanmean(waterdropstuffCP);
lasermeanCP=nanmean(laserstuffCP);

targetCPstats=anova1([whistlestuffCP waterdropstuffCP laserstuffCP]);

whistleCPSEM=nanstd(whistlestuffCP)/((sum(~isnan(whistlestuffCP)))^(1/2));
waterCPSEM=nanstd(waterdropstuffCP)/((sum(~isnan(waterdropstuffCP)))^(1/2));
laserCPSEM=nanstd(laserstuffCP)/((sum(~isnan(laserstuffCP)))^(1/2));


whistlestuffCNP=squeeze(sound_accuracy_bysound_all(4,1,:));
laserstuffCNP=squeeze(sound_accuracy_bysound_all(4,2,:));
waterdropstuffCNP=squeeze(sound_accuracy_bysound_all(4,3,:));


whistlemeanCNP=nanmean(whistlestuffCNP);
waterdropmeanCNP=nanmean(waterdropstuffCNP);
lasermeanCNP=nanmean(laserstuffCNP);

targetCNPstats=anova1([whistlestuffCNP waterdropstuffCNP laserstuffCNP]);


whistleCNPSEM=nanstd(whistlestuffCNP)/((sum(~isnan(whistlestuffCNP)))^(1/2));
waterCNPSEM=nanstd(waterdropstuffCNP)/((sum(~isnan(waterdropstuffCNP)))^(1/2));
laserCNPSEM=nanstd(laserstuffCNP)/((sum(~isnan(laserstuffCNP)))^(1/2));








%% plot that
figure
hold on
h=bar(validation_mean_1*100,'stacked');
Labels = {'Heard', 'Not heard'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
axis([0 3 0 100])
labels = {'Correct','Incorrect'};
legend (labels);
ylabel('Identification Accuracy (%)')
errorbar(validation_mean_1(:,1)'*100, validation_SEM*100, 'k.')

figure
hold on
h=bar(overall_accuracy*100,'stacked');
Labels = {'Heard', 'Not heard'};
set(gca, 'XTick', 1:2, 'XTickLabel', Labels);
axis([0 3 0 100])
labels = {'Correct','Incorrect'};
legend (labels);
ylabel('Identification Accuracy (%)')
errorbar(overall_accuracy'*100, validation_SEM*100, 'k.')

