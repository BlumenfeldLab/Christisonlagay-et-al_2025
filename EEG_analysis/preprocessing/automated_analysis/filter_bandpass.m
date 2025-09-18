function [ filtered_trials ] = filter_bandpass(fs, bandpass_start, bandpass_stop, BC_trials)
%filter_alpha will filter all the epochs in a data set according to the
%start and stop inputed into the filter

%   INPUTS
%   fs = sampling frequency
%   bandpass_start = frequency of the start of the bandpass filter
%   bandpass_stop = frequency of the stop of the bandpass filter
%   BC_trials = input trials (dimensions # electrode labels x # samples x trial #)

%   OUTPUTS
%   filtered_trials are the bandpass filtered trials that result

    
    %% set up the bandpass filter
    
    filterorder = 20;

    bpFilt = designfilt('bandpassiir','FilterOrder',filterorder, ...
         'HalfPowerFrequency1',bandpass_start,'HalfPowerFrequency2',bandpass_stop, ...
         'SampleRate',fs);
     
     
    %% use the filter
    filtered_trials = NaN(size(BC_trials));
           
    for i = 1:size(BC_trials,1)
        display(i)
        for k = 1:size(BC_trials,3)
            epoch = squeeze(BC_trials(i,:,k));
            % filter signal through each filter
            % place filtered results in output matrix
            filtered_trials(i,:,k) = filtfilt(bpFilt,epoch);
        end
    end
     

end

