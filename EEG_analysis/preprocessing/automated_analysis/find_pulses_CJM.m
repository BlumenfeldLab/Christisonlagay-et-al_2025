function [onsets, offsets] = find_pulses_CJM(ttl, fs, min_pulse_length, max_pulse_length, qp, tolerance, varargin)
%{
Created by Chris Micek on 9/10/2018

Finds the onset and offset indices of TTL pulses in the provided TTL
channel.

INPUTS:
  REQUIRED
  - ttl: The data from the EEG channel containing TTL pulses.

  - fs: The sampling rate of the TTL channel, in hertz.

  - min_pulse_length: The minimum valid pulse length, in seconds. Used to
                      exclude pulse onset and offset indices that are found
                      too close together. When paired with the 'trim'
                      option, (see 'params' below) excludes pulses shorter
                      than 'min_pulse_length' - 'tolerance'.

  - max_pulse_length: The maximum valid pulse length, in seconds. When
                      paired with the 'trim' option (see 'params' below),
                      excludes pulses with durations longer than the
                      'max_pulse_length' + 'tolerance'.

  - qp: A one- or two-element numeric vector containing the query point(s)
        used when using line intersections to detect pulse onsets/offsets.
        If 'qp' is one element, uses the same query point for detecting
        both onsets and offsets. Otherwise, the first element of 'qp'
        is the onset-detection query point, and the second is the
        offset-detection query point.

  - tolerance: Tolerance value included when determining if pulses have
               valid durations; a valid pulse has a duration >=
               'min_pulse_length' - 'tolerance' and <= 'max_pulse_length' +
               'tolerance'.

  OPTIONAL
  - params: A struct with fields containing any combination of the optional
            parameters below (passing in an empty struct is also possible,
            and can be useful if you want to clear existing fields you will
            redefine later):

            - reverse: One- or two-element boolean vector. If one element,
                       applies to both onset and offset detection; if two
                       elements, first element applies to onset detection
                       and second to offset detection. If true, detects
                       onsets/offsets from right to left instead of left
                       to right.

            - flip: One- or two-element boolean vector. If one element,
                    applies to both onset and offset detection; if two
                    elements, first element applies to onset detection and
                    second to offset detection. If true, flips detection
                    criterion for onsets/offsets (that is, assumes onsets
                    occur on a falling edge instead of a rising edge, and
                    vice versa for offsets).

            - strict_onset_detection: Removes any successive adjacent onsets
                                      from the list of discovered onsets and
                                      offsets, and assumes that of the
                                      detections that remain, the only valid
                                      offsets are adjacent to onsets. Mutually
                                      exclusive from 'strict_offset_detection'.

            - strict_offset_detection: Removes any successive adjacent offsets
                                       from the list of discovered onsets and
                                       offsets, and assumes that of the
                                       detections that remain, the only valid
                                       onsets are adjacent to offsets. Mutually
                                       exclusive from 'strict_onset_detection'.

            - trim: A single boolean flag indicating whether to trim
                    extraneous pulses from the list of pulse detections,
                    making these assumptions:
                    - The first detection should be an onset.
                    - The last detection should be an offset.
                    - Valid pulses have durations >= 'min_pulse_length' -
                      'tolerance' and <= 'max_pulse_length' +
                      'tolerance'.
                    - The minimum duration between successive pulses is
                      floor(100 * fs/4096) samples (~0.0244 seconds).

            - use_env: A single boolean flag indicating whether to use the
                       envelope of the TTL signal, instead of the signal
                       itself, for pulse detection. Pulses found in this
                       manner are then mapped to their true locations in
                       the raw signal.

            - sub_movmean: A single boolean, integer, or cell array of any
                           combination of booleans and integers. For each
                           non-zero value included, subtracts the moving
                           average of the TTL signal from the TTL signal.
                           If the nth element of sub_movmean is true
                           (logical 1), subtracts the moving average with
                           the default window size of 3 * fs samples for
                           the nth subtraction iteration. If the nth
                           element of sub_movmean is an integer, subtracts
                           the moving average with a window size of {n}
                           samples for the nth subtraction iteration.

OUTPUTS:
  - onsets: Vector of TTL pulse onset indices, with respect to the 'ttl'
            input variable.

  - offsets: Vector of TTL pulse offset indices, with respect to the 'ttl'
             input variable.
%}

% TTL Pulses for Run 2 for 415DO were rather noisy, so below is a more
% elaborate method of weeding out true pulse onsets/offsets from
% extraneous spikes due to noise (could replace the existing
% find_pulses function)

if numel(qp) < 1 || numel(qp) > 2
    throw(MException('find_pulses_CJM:IncorrectQueryPointNum', 'The number of query point values must be 1 or 2.'))
elseif numel(qp) == 1
    onset_query_point = qp;
    offset_query_point = qp;
elseif numel(qp) == 2
    onset_query_point = qp(1);
    offset_query_point = qp(2);
end

onset_reverse = false;
offset_reverse = false;

use_env = false;
strict_onset_detection = false;
strict_offset_detection = false;

trim = false;
sub_movmean = {false};

onset_flip = false;
offset_flip = false;

if numel(varargin) > 1
    throw(MException('find_pulses_CJM:IncorrectVararginNum', 'The number of variable input arguments must be 0 or 1.'))
elseif numel(varargin) == 1
    params = varargin{1};
    if isfield(params, 'reverse')
        reverse = params.reverse;
        if numel(reverse) > 2
            throw(MException('find_pulses_CJM:IncorrectReverseNum', 'The number of reverse values must be 1 or 2.'))
        elseif numel(reverse) == 1
            onset_reverse = reverse;
            offset_reverse = reverse;
        elseif numel(reverse) == 2
            onset_reverse = reverse(1);
            offset_reverse = reverse(2);
        end
    end
        
    if isfield(params, 'use_env')
        use_env = params.use_env;
        if numel(use_env) > 1 || ~islogical(use_env)
            throw(MException('find_pulses_CJM:IncorrectUse_Env', 'Use_env must be a single logical, either ''true'' or ''false''.'))
        end
    end
    
    if isfield(params, 'strict_onset_detection')
        strict_onset_detection = params.strict_onset_detection;
        if numel(strict_onset_detection) > 1 || ~islogical(strict_onset_detection)
            throw(MException('find_pulses_CJM:IncorrectStrict_Onset_Detection', 'Strict_onset_detection must be a single logical, either ''true'' or ''false''.'))
        end
    end
    
    if isfield(params, 'strict_offset_detection')
        strict_offset_detection = params.strict_offset_detection;
        if numel(strict_offset_detection) > 1 || ~islogical(strict_offset_detection)
            throw(MException('find_pulses_CJM:IncorrectStrict_Offset_Detection', 'Strict_offset_detection must be a single logical, either ''true'' or ''false''.'))
        end
    end
    
    if isfield(params, 'trim')
        trim = params.trim;
        if numel(trim) > 1 || ~islogical(trim)
            throw(MException('find_pulses_CJM:IncorrectTrim', 'Trim must be a single logical, either ''true'' or ''false''.'))
        end
    end
    
    if isfield(params, 'flip')
        flip = params.flip;
        if numel(flip) > 2
            throw(MException('find_pulses_CJM:IncorrectFlipNum', 'The number of flip values must be 1 or 2.'))
        elseif numel(flip) == 1
            onset_flip = flip;
            offset_flip = flip;
        elseif numel(flip) == 2
            onset_flip = flip(1);
            offset_flip = flip(2);
        end
    end
    
    % Added 3/4/2019 by CJM; subtracts a moving average of the TTL signal
    % to remove motion artifacts. The default number of samples to use is
    % 3 * fs samples, but an alternative can be specified in lieu of
    % 'true' or 'false'. Onsets and offsets are then calculated normally,
    % and mapped to the closest (lateral) points of the original signal at
    % the very end. This can be either a single value or a cell array of
    % any combination of logicals/integers; each element represents
    % successive iterations of averaging & subtracting.
    if isfield(params, 'sub_movmean')
        sub_movmean = params.sub_movmean;
        if ~iscell(sub_movmean)
            sub_movmean = {sub_movmean};
        end
        if (~any(cellfun(@islogical, sub_movmean)) && ~any(cellfun(@isnumeric, sub_movmean))) ...
                || any(cellfun(@numel, sub_movmean) > 1) || ...
                (any(cellfun(@isnumeric, sub_movmean)) && ...
                any(cellfun(@(x) isnumeric(x) && ~isequal(floor(x), x),...
                sub_movmean(cellfun(@isnumeric, sub_movmean)))))
            throw(MException('find_pulses_CJM:IncorrectSub_Movmean',...
                ['Sub_movmean must be contain only logicals (either ''true'' or ''false''),',...
                ' or integers indicating the number of samples to use in the moving average calculation.']))
        else
            ttl_copy = ttl;
        end
    end
end

for ii = 1:numel(sub_movmean)
    if sub_movmean{ii} || isnumeric(sub_movmean{ii})
        
        if isnumeric(sub_movmean{ii})
            k = sub_movmean{ii};
        else
            k = 3 * fs;
        end
        ttl = ttl - movmean(ttl, k);
    end
end

if onset_query_point >= 0 && ~onset_reverse && ~onset_flip
    onset_inds = [ttl, NaN] >= onset_query_point & [NaN, ttl] < onset_query_point;
else
    onset_inds = [ttl, NaN] < onset_query_point & [NaN, ttl] >= onset_query_point;
end

if onset_reverse
    onset_inds = fliplr(onset_inds);
end
onsets = find(onset_inds);
to_delete = NaN(1, numel(onsets));
for ii = 2:numel(onsets)
    if onsets(ii) - onsets(ii - 1) < min_pulse_length * fs
        to_delete(ii) = ii;
    end
end
onsets(~isnan(to_delete)) = [];
if onset_reverse
    onsets = numel(ttl) - fliplr(onsets);
end

if offset_query_point >= 0 && ~offset_reverse && ~offset_flip
    offset_inds = [ttl, NaN] < offset_query_point & [NaN, ttl] >= offset_query_point;
else
    offset_inds = [ttl, NaN] >= offset_query_point & [NaN, ttl] < offset_query_point;
end

if offset_reverse
    offset_inds = fliplr(offset_inds);
end
offsets = find(offset_inds);
to_delete = NaN(1, numel(offsets));
for ii = 2:numel(offsets)
    if offsets(ii) - offsets(ii - 1) < min_pulse_length * fs
        to_delete(ii) = ii;
    end
end
offsets(~isnan(to_delete)) = [];
if offset_reverse
    offsets = numel(ttl) - fliplr(offsets);
end

% At this point, we may not necessarily have the same number of
% onsets and offsets, and they may not necessarily appear one
% following the other like we would expect. To put things in order,
% we can associate each onset or offset sample with a label (1 for
% onset, 0 for offset), concatenate the sample numbers, and sort
% them in ascending order.

if ~use_env
    sorted_transitions = [onsets, offsets; ones(1, numel(onsets)), zeros(1, numel(offsets))];
    [~, order] = sort(sorted_transitions(1, :));
    sorted_transitions = sorted_transitions(:, order);
    if trim
        to_delete = NaN(1, size(sorted_transitions, 2));
        median_ttl = median(ttl);
        seek_samples = floor(100 * fs/4096);

        first_onset = find(sorted_transitions(2, :) == 1, 1);
        if first_onset > 1
            to_delete(1:first_onset - 1) = 1:first_onset - 1;
        end

        last_offset = find(sorted_transitions(2, :) == 0, 1, 'last');
        if last_offset < size(sorted_transitions, 2)
            to_delete(last_offset + 1:end) = last_offset + 1:size(sorted_transitions, 2);
        end

        for ii = 1:size(sorted_transitions, 2)

            if ii ~= size(sorted_transitions, 2) && sorted_transitions(2, ii) == 1
                for jj = (ii + 1):size(sorted_transitions, 2)
                    if sorted_transitions(2, jj) == 0
                        break
                    end
                end
                if (sorted_transitions(1, jj) - sorted_transitions(1, ii)) / fs > max_pulse_length + tolerance || ...
                        ~any(abs(ttl(max([1, sorted_transitions(1, ii) - seek_samples]):sorted_transitions(1, ii))) <= ...
                        abs(median_ttl) + tolerance) %|| ...
                    % (sorted_transitions(1, jj) - sorted_transitions(1, ii)) / fs < min_pulse_length - tol
                    to_delete(ii) = ii;
                end

            elseif ii ~= 1 && sorted_transitions(2, ii) == 0
                for jj = (ii - 1):-1:1
                    if sorted_transitions(2, jj) == 1
                        break
                    end
                end
                if (sorted_transitions(1, ii) - sorted_transitions(1, jj)) / fs < min_pulse_length - tolerance || ...
                        ~any(abs(ttl(sorted_transitions(1, ii):min([sorted_transitions(1, ii) + seek_samples, ...
                        size(ttl, 2)]))) <= abs(median_ttl) + tolerance) %|| ...
                    %(sorted_transitions(1, ii) - sorted_transitions(1, jj)) / fs > max_pulse_length + tol
                    to_delete(ii) = ii;
                end
            end

        end

        sorted_transitions(:, ~isnan(to_delete)) = [];
    end
    
    if strict_onset_detection
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 1 & [sorted_transitions(2, :), NaN] == 1)) = [];
        onsets = sorted_transitions(1, sorted_transitions(2, :) == 1);
        offsets = sorted_transitions(1, find(sorted_transitions(2, :) == 1) + 1);
    elseif strict_offset_detection
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 0 & [sorted_transitions(2, :), NaN] == 0)) = [];
        onsets = sorted_transitions(1, sorted_transitions(2, :) == 1);
        offsets = sorted_transitions(1, find(sorted_transitions(2, :) == 1) + 1);
    else
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 0 & [sorted_transitions(2, :), NaN] == 0)) = [];
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 1 & [sorted_transitions(2, :), NaN] == 1) - 1) = [];
        onsets = sorted_transitions(1, sorted_transitions(2, :) == 1);
        offsets = sorted_transitions(1, sorted_transitions(2, :) == 0);
    end
else
    
    env_onset_reverse = false;
    env_offset_reverse = false;
    
    if isfield(params, 'env_reverse')
        env_reverse = params.env_reverse;
        if numel(env_reverse) > 2
            throw(MException('find_pulses_CJM:IncorrectEnvReverseNum', 'The number of env_reverse values must be 1 or 2.'))
        elseif numel(env_reverse) == 1
            env_onset_reverse = env_reverse;
            env_offset_reverse = env_reverse;
        elseif numel(env_reverse) == 2
            env_onset_reverse = env_reverse(1);
            env_offset_reverse = env_reverse(2);
        end
    end
    
    env_onset_flip = false;
    env_offset_flip = false;
    
    if isfield(params, 'env_flip')
        env_flip = params.env_flip;
        if numel(env_flip) > 2
            throw(MException('find_pulses_CJM:IncorrectEnvFlipNum', 'The number of env_flip values must be 1 or 2.'))
        elseif numel(env_flip) == 1
            env_onset_flip = env_reverse;
            env_offset_flip = env_reverse;
        elseif numel(env_flip) == 2
            env_onset_flip = env_flip(1);
            env_offset_flip = env_flip(2);
        end
    end
    
    sorted_transitions = [onsets, offsets; ones(1, numel(onsets)), zeros(1, numel(offsets))];
    [~, order] = sort(sorted_transitions(1, :));
    sorted_transitions = sorted_transitions(:, order);

    to_delete = NaN(1, size(sorted_transitions, 2));

    first_onset = find(sorted_transitions(2, :) == 1, 1);
    if first_onset > 1
        to_delete(1:first_onset - 1) = 1:first_onset - 1;
    end

    last_offset = find(sorted_transitions(2, :) == 0, 1, 'last');
    if last_offset < size(sorted_transitions, 2)
        to_delete(last_offset + 1:end) = last_offset + 1:size(sorted_transitions, 2);
    end


    sorted_transitions(:, ~isnan(to_delete)) = [];

    onsets = sorted_transitions(1, sorted_transitions(2, :) == 1);
    offsets = sorted_transitions(1, sorted_transitions(2, :) == 0);

    % For extremely noisy signals, can use the signal envelope to eliminate
    % fast oscillations between pulses, check for transitions in envelope
    % function that are in the region of "true" transitions, and map back to
    % their actual locations on the TTL signal.

    zero_locs = ttl == 0;
    
    if sum(zero_locs) >= 0.5 * numel(ttl)
        ttl_2 = ttl + rand(1, numel(ttl)) .* zero_locs * max(ttl)/1e4;
        [upper, lower] = envelope(ttl_2, floor(200 * fs/4096), 'peak');
    else
        [upper, lower] = envelope(ttl, floor(200 * fs/4096), 'peak');
    end
    
    if onset_query_point >= 0
        onset_env = upper;
    else
        onset_env = lower;
    end

    if onset_query_point >= 0 && ~env_onset_reverse && ~env_onset_flip
        env_onset_inds = [onset_env, NaN] >= onset_query_point & [NaN, onset_env] < onset_query_point;
    else
        env_onset_inds = [onset_env, NaN] < onset_query_point & [NaN, onset_env] >= onset_query_point;
    end

    if env_onset_reverse
        env_onset_inds = fliplr(env_onset_inds);
    end
    env_onsets = find(env_onset_inds);

    % to_delete = NaN(1, numel(env_onsets));
    % for ii = 2:numel(env_onsets)
    %     if env_onsets(ii) - env_onsets(ii - 1) < min_pulse_length * fs
    %         to_delete(ii) = ii;
    %     end
    % end
    % env_onsets(~isnan(to_delete)) = [];

    if env_onset_reverse
        env_onsets = numel(onset_env) - fliplr(env_onsets);
    end

    if offset_query_point >= 0
        offset_env = upper;
    else
        offset_env = lower;
    end

    if offset_query_point >= 0 && ~env_offset_reverse && ~env_offset_flip
        env_offset_inds = [offset_env, NaN] < offset_query_point & [NaN, offset_env] >= offset_query_point;
    else
        env_offset_inds = [offset_env, NaN] >= offset_query_point & [NaN, offset_env] < offset_query_point;
    end

    if env_offset_reverse
        env_offset_inds = fliplr(env_offset_inds);
    end
    env_offsets = find(env_offset_inds);

%     to_delete = NaN(1, numel(offsets));
%     for ii = 2:numel(env_offsets)
%         if env_offsets(ii) - env_offsets(ii - 1) < min_pulse_length * fs
%             to_delete(ii) = ii;
%         end
%     end
%     env_offsets(~isnan(to_delete)) = [];

    if env_offset_reverse
        env_offsets = numel(offset_env) - fliplr(env_offsets);
    end

    % env_onsets = find([upper, 0] >= onset_query_point & [0, upper] < onset_query_point);
    % env_offsets = find([upper, 0] < offset_query_point & [0, upper] >= offset_query_point);

    true_onsets = NaN(1, numel(env_onsets));
    true_offsets = NaN(1, numel(env_offsets));

    for ii = 1:numel(true_onsets)
        [~, onset_ind] = min(abs(onsets - env_onsets(ii)));
        true_onsets(ii) = onsets(onset_ind);
    end
    for ii = 1:numel(true_offsets)
        [~, offset_ind] = min(abs(offsets - env_offsets(ii)));
        true_offsets(ii) = offsets(offset_ind);
    end

    onsets = true_onsets;
    offsets = true_offsets;

    % The envelope method above can lead to some false positives, so we should
    % clean up our pulse detections so they have durations we expect.
    sorted_transitions = [onsets, offsets; ones(1, numel(onsets)), zeros(1, numel(offsets))];
    [~, order] = sort(sorted_transitions(1, :));
    sorted_transitions = sorted_transitions(:, order);
    to_delete = NaN(1, size(sorted_transitions, 2));

    for ii = 1:size(sorted_transitions, 2)

        if ii ~= size(sorted_transitions, 2) && sorted_transitions(2, ii) == 1
            for jj = (ii + 1):size(sorted_transitions, 2)
                if sorted_transitions(2, jj) == 0
                    break
                end
            end
            if (sorted_transitions(1, jj) - sorted_transitions(1, ii)) / fs > max_pulse_length + tolerance
                to_delete(ii) = ii;
            end

        elseif ii ~= 1 && sorted_transitions(2, ii) == 0
            for jj = (ii - 1):-1:1
                if sorted_transitions(2, jj) == 1
                    break
                end
            end
            if (sorted_transitions(1, ii) - sorted_transitions(1, jj)) / fs < min_pulse_length - tolerance
                to_delete(ii) = ii;
            end
        end

    end
    
    first_onset = find(sorted_transitions(2, :) == 1, 1);
    if first_onset > 1
        to_delete(1:first_onset - 1) = 1:first_onset - 1;
    end

    last_offset = find(sorted_transitions(2, :) == 0, 1, 'last');
    if last_offset < size(sorted_transitions, 2)
        to_delete(last_offset + 1:end) = last_offset + 1:size(sorted_transitions, 2);
    end

    sorted_transitions(:, ~isnan(to_delete)) = [];
    if strict_onset_detection
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 1 & [sorted_transitions(2, :), NaN] == 1)) = [];
        onsets = sorted_transitions(1, sorted_transitions(2, :) == 1);
        offsets = sorted_transitions(1, find(sorted_transitions(2, :) == 1) + 1);
    elseif strict_offset_detection
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 0 & [sorted_transitions(2, :), NaN] == 0)) = [];
        onsets = sorted_transitions(1, sorted_transitions(2, :) == 1);
        offsets = sorted_transitions(1, find(sorted_transitions(2, :) == 1) + 1);
    else
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 0 & [sorted_transitions(2, :), NaN] == 0) - 1) = [];
        sorted_transitions(:, find([NaN, sorted_transitions(2, :)]  == 1 & [sorted_transitions(2, :), NaN] == 1)) = [];
        onsets = sorted_transitions(1, sorted_transitions(2, :) == 1);
        offsets = sorted_transitions(1, sorted_transitions(2, :) == 0);
    end
end

% if sub_movmean || isnumeric(sub_movmean)
%     if onset_query_point >= 0
%         true_onset_inds = find([ttl_copy, NaN] >= onset_query_point & [NaN, ttl_copy] < onset_query_point);
%     else
%         true_onset_inds = find([ttl_copy, NaN] < onset_query_point & [NaN, ttl_copy] >= onset_query_point);
%     end
%     
%     for ii = 1:numel(onsets)
%         [~, min_ind] = min(abs(true_onset_inds - onsets(ii)));
%         onsets(ii) = true_onset_inds(min_ind);
%     end
%     
%     if offset_query_point >= 0
%         true_offset_inds = find([ttl_copy, NaN] < offset_query_point & [NaN, ttl_copy] >= offset_query_point);
%     else
%         true_offset_inds = find([ttl_copy, NaN] >= offset_query_point & [NaN, ttl_copy] < offset_query_point);
%     end
%     
%     for ii = 1:numel(offsets)
%         [~, min_ind] = min(abs(true_offset_inds - offsets(ii)));
%         offsets(ii) = true_offset_inds(min_ind);
%     end
%     
% end

end