function [ filtered_trials ] = filter_lowpass( fs, lowpass_cutoff, BC_trials)

%   fs = sampling frequency
%   lowpass_cutoff = frequency of the start of the lowpass filter
%   BC_trials = input trials (dimensions # electrode labels x # samples x trial #)
%   labels = labels for electrodes, should match first dimension of
%   BC_trials

%   output filtered_trials should be the same dimensions as original
%   BC_trials but the samples should now be filtered.


    % create low-pass filter at the frequency specified by the input
    lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',lowpass_cutoff,'PassbandRipple',0.2, ...
         'SampleRate',fs);

    filtered_trials = NaN(size(BC_trials));
           
    for i = 1:size(BC_trials,1)
        display(i)
        for k = 1:size(BC_trials,3)
            epoch = squeeze(BC_trials(i,:,k));
            % filter signal through each filter
            % place filtered results in output matrix
            filtered_trials(i,:,k) = filtfilt(lpFilt,epoch);
        end
    end

end

