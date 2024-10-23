function [epoch_numbers, removed_epochs_total] = ...
    removed_epochs_count(data, P, P1, P3, P6, ...
    err_P1_flx_sorted, err_P3_flx_sorted, err_P6_flx_sorted)
    
    epoch_numbers = zeros(1, length(data));
    for i = 1:length(P1.trials)
        l = length(data{1, P1.trials(i)}.EMG_stream.Sensors_Preprocessed);
        epoch_numbers(P1.trials(i)) = l;
    end
    for i = 1:length(P3.trials)
        l = length(data{1, P3.trials(i)}.EMG_stream.Sensors_Preprocessed);
        epoch_numbers(P3.trials(i)) = l;
    end
    for i = 1:length(P6.trials)
        l = length(data{1, P6.trials(i)}.EMG_stream.Sensors_Preprocessed);
        epoch_numbers(P6.trials(i)) = l;
    end
    
    K = floor(P*size(err_P1_flx_sorted, 2));
    removed_epochs_P1 = err_P1_flx_sorted(1:4, K+1:end, 2);
    K = floor(P*size(err_P3_flx_sorted, 2));
    removed_epochs_P3 = err_P3_flx_sorted(1:4, K+1:end, 2);
    K = floor(P*size(err_P6_flx_sorted, 2));
    removed_epochs_P6 = err_P6_flx_sorted(1:4, K+1:end, 2);
    
    removed_epochs_total = zeros(4, length(data));
    for i = 1:4
    
        for j = 1:size(removed_epochs_P1,2)
            removed_epochs_total(i, removed_epochs_P1(i,j)) = ...
                removed_epochs_total(i, removed_epochs_P1(i,j)) + 1;
        end
    
        for j = 1:size(removed_epochs_P3,2)
            removed_epochs_total(i, removed_epochs_P3(i,j)) = ...
                removed_epochs_total(i, removed_epochs_P3(i,j)) + 1;
        end
    
        for j = 1:size(removed_epochs_P6,2)
            removed_epochs_total(i, removed_epochs_P6(i,j)) = ...
                removed_epochs_total(i, removed_epochs_P6(i,j)) + 1;
        end
        
    end
    
    

end