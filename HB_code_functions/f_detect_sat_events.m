function [desat, resat, desatResatOK] = f_detect_sat_events(satSignal, satFreq, o2Threshold, desat_duration, resat_duration, resat_drop_tresh)
% This functions detects desaturations and resaturations.

% Set parameters
o2MaxDur = desat_duration * satFreq; % this is the maximal desat search window
resatMaxDur = resat_duration * satFreq; % this is the maximal resat search window, Azarbarzin uses 100 s

% Calculate the difference between signal points
satDeriv = diff(satSignal);

% Detect desaturations
desat = zeros(1,length(satSignal)-1);
nLimit = 0;

% Loop to detect desaturations based on satDeriv
point_start = 1;
point_end = (length(satSignal)-2)-o2MaxDur*satFreq; %last seconds cannot be analyzed
for i = point_start:point_end
    dur = 0;
    if satDeriv(i) < 0
        while (satDeriv(i + dur) <= 0 && dur <= o2MaxDur)
            o2Level = satSignal(i) - satSignal(i+dur+1);
            dur = dur + 1;
            if o2Level < o2Threshold
                nLimit = dur + 1;
            end
        end
        if nLimit <= o2MaxDur-1 && o2Level >= o2Threshold
            desat(i+1:i+dur) = 1;
        end
    end
    % pause(0.00001) % the short pause might be needed, because otherwise we cannot interrupt the loop with C-c
end

% Cleanup
clear satDeriv o2MaxDur dur nLimit o2Level

% Calculate some elements from desat vector
desatStart = find(diff(desat) == 1);
desatEnd = find(diff(desat) == -1);
%desatDur = desatEnd - desatStart;
desatNumber = length(desatEnd);
desatDrop = satSignal(desatStart) - satSignal(desatEnd);
%meanDrop = mean(desatDrop);
resatThreshold = desatDrop-(resat_drop_tresh*desatDrop);

% Detect resaturations
resat = zeros(1, length(desat));
desatResatOK = zeros(1, desatNumber);

% Loop to detect resaturations based on satDeriv
for i = 1:desatNumber
    resatDur = 0;
    % don't check if the event is too close to the end of the file
    if desatEnd(i) + resatMaxDur + 2 > length(satSignal)
        break
    else

        while resatDur <= resatMaxDur && ...
                satSignal(desatEnd(i) + resatDur) < satSignal(desatStart(i)) - resatThreshold(i);
            resatDur = resatDur + 1;
            %         if desatEnd(i) + resatDur > length(satSignal)
            %             %don't search beyond the length of the file
            %             break
            %         end
            if i < desatNumber && desatEnd(i) + resatDur + 1 > desatStart(i+1) % events shouldn't overlap
                break
            end
        end

        while resatDur <= resatMaxDur && ...
                satSignal(desatEnd(i) + resatDur -1) >= satSignal(desatStart(i)) - resatThreshold(i)
            resatDur = resatDur - 1;
            if i < desatNumber && desatEnd(i) + resatDur + 1 > desatStart(i+1) % events still shouldn't overlap
                break
            end
        end
        resatFlag = 0;
        % this I don't quite understand: shouldn't we set the resatFlag when it
        % is *larger* than the threshold??
        %     if resatDur <= resatMaxDur && ...
        %             satSignal(desatEnd(i) + resatDur) < satSignal(desatStart(i)) - resatThreshold(i)
        %         resatFlag = 1;
        if resatDur <= resatMaxDur && ...
                satSignal(desatEnd(i) + resatDur) >= satSignal(desatStart(i)) - resatThreshold(i)
            resatFlag = 1;
        end

        while satSignal(desatEnd(i) + resatDur) <= satSignal(desatEnd(i) + resatDur + 1) && ...
                satSignal(desatStart(i)) >= satSignal(desatEnd(i) + resatDur + 1 ) && ...
                resatFlag == 1
            resatDur = resatDur + 1;
            if desatEnd(i) + resatDur + 1 > length(satSignal) - 1
                break
            end
            if i < desatNumber && desatEnd(i) + resatDur + 1 > desatStart(i+1)
                break
            end
        end

        while satSignal(desatEnd(i) + resatDur) <= satSignal(desatEnd(i) + resatDur - 1)
            resatDur = resatDur - 1;
        end
        if resatFlag == 1 && resatDur > 1 % this prevents zero-length events to be scored
            resat(desatEnd(i)+1:desatEnd(i) + resatDur) = 1;
            desatResatOK(i) = 1;
        end
    end
end
