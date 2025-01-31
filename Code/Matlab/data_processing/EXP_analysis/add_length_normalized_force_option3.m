function [Force_data, p1_extension_start_percent, p3_extension_start_percent, p6_extension_start_percent] = ...
    add_length_normalized_force_option3(Force_data, condition_indices, ignore_trials, ...
    Lepoch_p1_median, Lepoch_p3_median, Lepoch_p6_median)

    p1_extension_start_percent = [];
    p3_extension_start_percent = [];
    p6_extension_start_percent = [];
   
    for i = 1:length(Force_data)
        if ~ismember(i, ignore_trials)
            if ismember(i, condition_indices.P1)
                events_indxs = Force_data{1, i}.Events;
                for j = 1:length(Force_data{1, i}.F_cal_not_length_normalized)
                    force = Force_data{1, i}.F_cal_not_length_normalized{j};
                    Force_data{1, i}.F_cal_length_normalized = cat(1, ...
                        Force_data{1, i}.F_cal_length_normalized, ...
                        [interp1(events_indxs(j, 1):events_indxs(j, 3), force, ...
                            linspace(events_indxs(j, 1), events_indxs(j, 3), Lepoch_p1_median))]);
                    
                    [~, indx] = min( abs(Force_data{1, i}.F_cal_length_normalized(j, :) - ...
                        force(1, events_indxs(j, 2) - events_indxs(j, 1) + 1)) );
                    p1_extension_start_percent = cat(1, p1_extension_start_percent, indx/Lepoch_p1_median);
                end
                
            elseif ismember(i, condition_indices.P3)
                events_indxs = Force_data{1, i}.Events;
                for j = 1:length(Force_data{1, i}.F_cal_not_length_normalized)
                    force = Force_data{1, i}.F_cal_not_length_normalized{j};
                    Force_data{1, i}.F_cal_length_normalized = cat(1, ...
                        Force_data{1, i}.F_cal_length_normalized, ...
                        [interp1(events_indxs(j, 1):events_indxs(j, 3), force, ...
                            linspace(events_indxs(j, 1), events_indxs(j, 3), Lepoch_p3_median))]);
    
                    [~, indx] = min( abs(Force_data{1, i}.F_cal_length_normalized(j, :) - ...
                        force(1, events_indxs(j, 2) - events_indxs(j, 1) + 1)) );
                    p3_extension_start_percent = cat(1, p3_extension_start_percent, indx/Lepoch_p3_median);
                end
                
            elseif ismember(i, condition_indices.P6)
                events_indxs = Force_data{1, i}.Events;
                for j = 1:length(Force_data{1, i}.F_cal_not_length_normalized)
                    force = Force_data{1, i}.F_cal_not_length_normalized{j};
                    Force_data{1, i}.F_cal_length_normalized = cat(1, ...
                        Force_data{1, i}.F_cal_length_normalized, ...
                        [interp1(events_indxs(j, 1):events_indxs(j, 3), force, ...
                            linspace(events_indxs(j, 1), events_indxs(j, 3), Lepoch_p6_median))]);
    
    
                    [~, indx] = min( abs(Force_data{1, i}.F_cal_length_normalized(j, :) - ...
                        force(1, events_indxs(j, 2) - events_indxs(j, 1) + 1)) );
                    p6_extension_start_percent = cat(1, p6_extension_start_percent, indx/Lepoch_p6_median);
                end
                
            end
        end
    end

end