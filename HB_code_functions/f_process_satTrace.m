function [new_satSignal, Error_satSignal_defect] = f_process_satTrace(satSignal, satFreq, ...
    n_sat_val, ssec, change_time, change_range, ...
    filt_cutoff_diff, filt_low_range, filt_high_range, ...
    nan_min_dur, segment_tresh, segment_time, segments_keep, segments_all_low)
% This function cleans the saturation signal in several steps.
% This is not the optimal cleaning, however, it could be OK
% for the price/quality ratio.

%% Step 0
% Rough check of signal quality
Error_satSignal_defect = [];

cond_1 = length(unique(satSignal)) < n_sat_val;
cond_2 = (sum(satSignal >= filt_high_range) + sum(filt_low_range >= satSignal))/length(satSignal)*100 > 70;

if cond_1 | cond_2 == 1
    Error_satSignal_defect = 1;
    new_satSignal = [];
else
    %% Step 1
    % If the signal stays the same for more than ssec
    % than it is likely to be an artifact and it should be removed
    stable_interval = ssec*satFreq;
    new_satSignal = satSignal;
    Error_satSignal_defect = 0;
    for ll = 1:length(satSignal)
        % For first stable_interval points
        if stable_interval >= ll
            sigTemp_forward = satSignal(ll:(ll+stable_interval));
            sigTemp_backward = satSignal(1:ll);
            range_forward = range(sigTemp_forward);
            range_backward = range(sigTemp_backward);
            % For middle stable_interval points
        elseif ll > stable_interval & ll < (length(satSignal)-stable_interval)
            sigTemp_forward = satSignal(ll:(ll+stable_interval));
            sigTemp_backward = satSignal((ll-stable_interval):ll);
            range_forward = range(sigTemp_forward);
            range_backward = range(sigTemp_backward);
            % For last stable_interval points
        else
            sigTemp_forward = satSignal(ll:length(satSignal));
            sigTemp_backward = satSignal((ll-stable_interval):ll);
            range_forward = range(sigTemp_forward);
            range_backward = range(sigTemp_backward);
        end

        if range_forward > 0 & range_backward > 0
            new_satSignal(ll) = new_satSignal(ll);
        else
            new_satSignal(ll) = NaN;
            %plot(satSignal((ll-stable_interval):(ll+stable_interval)))
        end

    end % loop
    satSignal = new_satSignal;
    %satSignal_flat_removed = new_satSignal;

    %% Step 2
    % If the segment changes more than 4% within 1 second
    % replace with NaN
    change_interval = change_time*satFreq;
    new_satSignal = satSignal;
    for ll = 1:length(satSignal)
        % For first stable_interval points
        if change_interval >= ll
            sigTemp_forward = satSignal(ll:(ll+change_interval));
            sigTemp_backward = satSignal(1:ll);
            range_forward = range(sigTemp_forward);
            range_backward = range(sigTemp_backward);
            % For middle stable_interval points
        elseif ll > change_interval & ll < (length(satSignal)-change_interval)
            sigTemp_forward = satSignal(ll:(ll+change_interval));
            sigTemp_backward = satSignal((ll-change_interval):ll);
            range_forward = range(sigTemp_forward);
            range_backward = range(sigTemp_backward);
            % For last stable_interval points
        else
            sigTemp_forward = satSignal(ll:length(satSignal));
            sigTemp_backward = satSignal((ll-change_interval):ll);
            range_forward = range(sigTemp_forward);
            range_backward = range(sigTemp_backward);
        end

        if range_forward >= change_range | range_backward >= change_range
            new_satSignal(ll) = NaN;
        else
            new_satSignal(ll) = new_satSignal(ll);
            %plot(satSignal((ll-stable_interval):(ll+stable_interval)))
        end

    end % loop
    satSignal = new_satSignal;
    satSignal_change_removed = new_satSignal;

    %% Step 3
    % Remove outlier and replace with NA
    drop5 = abs(diff(satSignal)) >= filt_cutoff_diff;
    satSignal(drop5) = nan;
    satSignal(satSignal >= filt_high_range) = nan;
    satSignal(filt_low_range >= satSignal) = nan;
    %satSignal_highlow_removed = satSignal;

    %% Step 4
    % Replace small segments with NaNs
    % using the average of the surrounding values
    % This is only for small segments
    [roww, coll] = find(isnan(satSignal));
    find_segments = abs(diff(roww));
    find_segment_index = find(find_segments > 1);
    n_segments = length(find_segment_index)+1;

    nan_min_points = nan_min_dur*satFreq;

    for nn = 1:n_segments;

        if nn == 1
            seg_start = 1;
            seg_end = find_segment_index(nn);
        elseif nn == n_segments
            seg_start = find_segment_index(nn-1) + 1;
            seg_end = length(roww);
        else
            seg_start = find_segment_index(nn-1) + 1;
            seg_end = find_segment_index(nn);
        end

        seg_start = roww(seg_start);
        seg_end = roww(seg_end);

        seg_duration = seg_end-seg_start;
        if seg_duration > nan_min_points
            % Leave as it is
        else
            % Calculate the mean of adjacent values (+/- 1 sec) to fill NAs
            second_back_start = seg_start - satFreq;
            second_back_end = seg_start - 1;
            second_after_start = seg_end + 1;
            second_after_end = seg_end + satFreq;

            if satFreq >= seg_start % if starts less than a second from the start
                mean_fill = mean(satSignal(second_after_start:second_after_end));
            elseif (length(satSignal) - seg_end) < satFreq % if ends less than a second from the end
                mean_fill = mean(satSignal(second_back_start:second_back_end));
            else
                mean_before = mean(satSignal(second_after_start:second_after_end));
                mean_after = mean(satSignal(second_back_start:second_back_end));
                mean_fill = (mean_before+mean_after)/2;
            end

            satSignal(seg_start:seg_end) = mean_fill;
        end %if

    end
    satSignal_na_refilled = satSignal;

    %% Step 5
    % If there are segments suspicious for artifacts, remove them
    [roww, coll] = find(~isnan(satSignal));
    find_segments = abs(diff(roww));
    find_segment_index = find(find_segments > 1);
    n_segments = length(find_segment_index)+1;

    nonnan_min_points = segments_keep*satFreq;

    segment_tresh = nanmean(satSignal)-segment_tresh;

    if n_segments == 1;
        % Do nothing
    else
        for nn = 1:n_segments;
            % Define start and end
            if nn == 1
                seg_start = 1;
                seg_end = find_segment_index(nn);
            elseif nn == n_segments
                seg_start = find_segment_index(nn-1) + 1;
                seg_end = length(roww);
            else
                seg_start = find_segment_index(nn-1) + 1;
                seg_end = find_segment_index(nn);
            end

            seg_start = roww(seg_start);
            seg_end = roww(seg_end);

            % Compute duration
            seg_duration = seg_end-seg_start+1;

            % Compute properties
            satTemp = satSignal(seg_start:seg_end);
            satTemp_sel = satTemp(satTemp > segment_tresh);
            cut_off_val = length(satTemp_sel)/length(satTemp)*100;

            if nonnan_min_points > seg_duration
                satSignal(seg_start:seg_end) = NaN;
            elseif segment_time > cut_off_val
                satSignal(seg_start:seg_end) = NaN;
            else
                % Do nothing
            end %if

        end
    end

    %new_satSignal = satSignal;
    %satSignal_small_segs_removed = satSignal;

    %% Step 6
    % If there are still isolated segments below 90, they are likely to be
    % artifacts and should be removed
    [roww, coll] = find(~isnan(satSignal));
    find_segments = abs(diff(roww));
    find_segment_index = find(find_segments > 1);
    n_segments = length(find_segment_index)+1;

    % If there is more than one segment
    if n_segments > 1
        for nn = 1:n_segments;
            % Define start and end
            if nn == 1
                seg_start = 1;
                seg_end = find_segment_index(nn);
            elseif nn == n_segments
                seg_start = find_segment_index(nn-1) + 1;
                seg_end = length(roww);
            else
                seg_start = find_segment_index(nn-1) + 1;
                seg_end = find_segment_index(nn);
            end

            seg_start = roww(seg_start);
            seg_end = roww(seg_end);

            % Compute duration
            if sum((satSignal(seg_start:seg_end))<segments_all_low) == length(satSignal(seg_start:seg_end))
                satSignal(seg_start:seg_end) = NaN;
            else
                %
            end
        end
        % If there is only one segment
    elseif n_segments == 1
        if sum((satSignal(1:length(satSignal)))<segments_all_low) == length(satSignal(seg_start:seg_end))
            satSignal(1:length(satSignal)) = NaN;
        else
            %
        end
    else
        %
    end
    new_satSignal = satSignal;
end
end