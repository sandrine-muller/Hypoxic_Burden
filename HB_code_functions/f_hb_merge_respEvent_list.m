function [respEventList_merged] = f_hb_merge_respEvent_list(respEventList, FI_minimal_distance_between_events)
% This function merges the list with respiratory events
% so that it could be used for calculations according to Filchenko method

respEventList_merged = respEventList;
tresh_val = 0;

n_events = height(respEventList_merged);


if n_events > 1
    while tresh_val < FI_minimal_distance_between_events

        n_events = height(respEventList_merged);

        for rr = 1:(n_events-1)
            respEventList_merged.dif_val(rr) = respEventList_merged.ends_relative_sec(rr+1) - respEventList_merged.starts_relative_sec(rr);
        end

        for rr = 1:(n_events-1)
            if respEventList_merged.dif_val(rr) < FI_minimal_distance_between_events
                respEventList_merged.ends_relative_sec(rr) = respEventList_merged.ends_relative_sec(rr+1);
                respEventList_merged.delete(rr+1) = 1;
            else
                respEventList_merged.delete(rr+1) = 0;
            end
        end

        respEventList_merged.respEvduration_rec = respEventList_merged.ends_relative_sec - respEventList_merged.starts_relative_sec;
        respEventList_merged = respEventList_merged(respEventList_merged.delete == 0,:);
        n_events = height(respEventList_merged)-1;

        for rr = 1:(n_events-1)
            respEventList_merged.dif_val(rr) = respEventList_merged.ends_relative_sec(rr+1) - respEventList_merged.starts_relative_sec(rr);
        end

        tresh_val = min(respEventList_merged.dif_val(1:n_events));

    end % while
else
    respEventList_merged.delete = 0;
    respEventList_merged.dif_val = FI_minimal_distance_between_events-1;
end

end