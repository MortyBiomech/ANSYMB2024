function condition_indices = condition_indices_identifier(Trials_Info, subject)

    P1 = []; P3 = []; P6 = [];
    for i = 1:length(Trials_Info)
        P = Trials_Info{1, i}.General.Pressure;
        switch P
            case 1
                if subject >= 10
                    if strcmp(Trials_Info{1, i}.General.Description, 'Experiment')
                        P1 = cat(2, P1, i);
                    end
                else
                    P1 = cat(2, P1, i);
                end
            case 3
                if subject >= 10
                    if strcmp(Trials_Info{1, i}.General.Description, 'Experiment')
                        P3 = cat(2, P3, i);
                    end
                else
                    P3 = cat(2, P3, i);
                end
            case 6
                if subject >= 10
                    if strcmp(Trials_Info{1, i}.General.Description, 'Experiment')
                        P6 = cat(2, P6, i);
                    end
                else
                    P6 = cat(2, P6, i);
                end
        end
    end
    
    condition_indices.P1 = P1;
    condition_indices.P3 = P3;
    condition_indices.P6 = P6;
end