function [AZB_area_sum, AZB_duration_mean, AZB_duration_sum] = f_auc_azarbarzin_calc(treshold_vals, ...
    satSignal, satFreq, respEventList_merged, startVecNew, endVecNew, ...
    AZB_window_beforeRE, AZB_window_afterRE, AZB_pp)
% This function is used to calculate the area under the curve
% according to the method by Azabarzin

satSignal_nonnan = satSignal(~isnan(satSignal));
% Calculate the duration of the recording
duration_rec = length(satSignal_nonnan)/satFreq/60/60;

% Set baseline parameters
AZB_window = AZB_window_beforeRE + AZB_window_afterRE;
AZB_area = [];
AZB_duration_mean = [];
AZB_duration_sum = [];

% Loop over the events to calculate parameters per event
for i = 1:height(respEventList_merged)

    % If the event has no treshold, it is likely to be
    % in a zone of artifacts and should not be assessde
    if isnan(treshold_vals(i)) == 1
        duration_under_the_curve = NaN;
        area_under_the_curve = NaN;
    else

        % Cut the piece of the signal
        find_start = startVecNew(i)-(AZB_window_beforeRE*satFreq);
        find_end = endVecNew(i)+(AZB_window_afterRE*satFreq);

        % Check if event is still within the signal

        if find_start < 1 | find_end > length(satSignal)
            duration_under_the_curve = NaN;
            area_under_the_curve = NaN;
        else
            satSignTmp = satSignal(find_start:find_end);

            % If the signal includes NaN, it should not be assessed
            if anynan(satSignTmp) == 1
                duration_under_the_curve = NaN;
                area_under_the_curve = NaN;
            else

                satSignTmp = [treshold_vals(i);satSignTmp;treshold_vals(i)];

                % If the signal exceeds treshold_vals, replace with treshold_vals
                satSignTmp(satSignTmp > treshold_vals(i)) = treshold_vals(i);

                % If the signal is the same, then the square and duration AUC is 0
                % Check it in the following statement
                if range(satSignTmp) == 0
                    duration_under_the_curve = 0;
                    area_under_the_curve = 0;
                else
                    % Find where the signal reaches the treshold
                    a = find(satSignTmp == treshold_vals(i));
                    b = diff(a);
                    c = find(b > 1);
                    d = [1];
                    for cc = 1:width(c)
                        d = [d;a(c+1)];
                    end
                    treshold_pos = d;

                    % Find the most probable point of the resaturation
                    % 1. Check if the resaturation to baseline is achieved
                    % from the second half of the respiratory event
                    % to pp% from the event duration - then it belongs to this very event
                    find_start = AZB_window_beforeRE*satFreq;
                    find_end = length(satSignTmp) - AZB_window_afterRE*satFreq;
                    midline_position = find_start + (find_end-find_start)/2;
                    event_dur = find_end - find_start;
                    event_after_window = event_dur*AZB_pp/100;
                    right_border = find_end + event_after_window;
                    % Here we take from 1 since we assume potential scoring issue
                    %satSignTmp_event = satSignTmp(1:right_border);

                    sel_treshold_pos = treshold_pos(treshold_pos >= midline_position & right_border >= treshold_pos);

                    if ~isempty(sel_treshold_pos)
                        % Take the first value
                        sel_treshold_pos = sel_treshold_pos(1);
                        sel_treshold_type = 1;
                    else
                        % 2. Check if the resaturation to baseline is achieved
                        % from pp% from the event duration until the end of the window
                        % - then it belongs to this very event
                        sel_treshold_pos = treshold_pos(treshold_pos > right_border);
                        if ~isempty(sel_treshold_pos)
                            % Take the first value
                            sel_treshold_pos = sel_treshold_pos(1);
                            sel_treshold_type = 2;
                        else
                            % 3. Check if the resaturation to baseline is achieved
                            % from the start to the second half of the respiratory event
                            % - then it belongs to this very event
                            sel_treshold_pos = treshold_pos(treshold_pos > find_start & midline_position > treshold_pos);
                            if ~isempty(sel_treshold_pos)
                                % Take the last value
                                sel_treshold_pos = sel_treshold_pos(length(sel_treshold_pos));
                                sel_treshold_type = 3;
                            else
                                % 4. Check if the resaturation to baseline is achieved
                                % after pp% from the event duration until the end of the window
                                % - then it belongs to this very event
                                sel_treshold_pos = treshold_pos(treshold_pos > right_border);
                                if ~isempty(sel_treshold_pos)
                                    sel_treshold_pos = sel_treshold_pos(length(sel_treshold_pos));
                                    sel_treshold_type = 4;
                                else
                                    sel_treshold_pos = length(satSignTmp);
                                    sel_treshold_type = 5;
                                end
                            end
                        end
                    end

                    % Select the signal from the start
                    % until the selected resaturation point
                    satSignTmp_calc = satSignTmp(1:sel_treshold_pos);

                    % Find the area under the curve
                    % The baseline needs to be subtracted and the signal flipped for
                    % positive peaks (in fact, the result will be the same)
                    satSignTmp_calc_AZB = -(satSignTmp_calc-treshold_vals(i));
                    area_under_the_curve = trapz(1/(satFreq*60),satSignTmp_calc_AZB);
                    duration_under_the_curve = sum(satSignTmp_calc_AZB > 0)/satFreq; % in seconds

                end % if

            end % signal

        end % treshold

    end % signal borders are not within event

    % Add to the common matrix
    AZB_area = [AZB_area; area_under_the_curve];
    AZB_duration = [AZB_area; duration_under_the_curve];

end % event loop

% Save the values with relation to the signal duration
AZB_area_sum = nansum(AZB_area)/duration_rec; %
AZB_duration_mean = nanmean(AZB_duration(AZB_duration > 0)); % seconds
AZB_duration_sum = (nansum(AZB_duration)/60/60)/duration_rec*100; % recording time

end