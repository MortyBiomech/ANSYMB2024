function output = Torque_Force_timeWarping(Trials_Info, Knee_Torque_Force_Sensor, Times, data,  subject)

    output = Knee_Torque_Force_Sensor;


    %% time warping of the force sensor data 
    flex2flex_start_indx = cellfun(@(x) x.Events.EXP_stream.flextoflex_start_indx, ...
        Trials_Info, 'UniformOutput', false);

    extension_start_indx = cellfun(@(x, y) x.Events.EXP_stream.extension_start_indx(1:length(y)), ...
        Trials_Info, flex2flex_start_indx, 'UniformOutput', false);

    flex2flex_end_indx = cellfun(@(x, y) x.Events.EXP_stream.flextoflex_end_indx(1:length(y)), ...
        Trials_Info, flex2flex_start_indx, 'UniformOutput', false);

    if size(flex2flex_start_indx{1,1}, 2) == 1
        flex2flex_start_indx = cellfun(@(x) x', flex2flex_start_indx, 'UniformOutput', false);
    end
    if size(extension_start_indx{1,1}, 2) == 1
        extension_start_indx = cellfun(@(x) x', extension_start_indx, 'UniformOutput', false);
    end
    if size(flex2flex_end_indx{1,1}, 2) == 1
        flex2flex_end_indx = cellfun(@(x) x', flex2flex_end_indx, 'UniformOutput', false);
    end

    flex2flex_start_indx = cell2mat(flex2flex_start_indx)';
    extension_start_indx = cell2mat(extension_start_indx)';
    flex2flex_end_indx = cell2mat(flex2flex_end_indx)';

    flexion_parts = extension_start_indx - flex2flex_start_indx + 1;
    extension_parts = flex2flex_end_indx - extension_start_indx;

    flexion_length = floor(median(flexion_parts));
    extension_length = floor(median(extension_parts));
    upper_flexion_length = flexion_length + 3*std(flexion_parts);
    lower_flexion_length = flexion_length - 3*std(flexion_parts);
    upper_extension_length = extension_length + 3*std(extension_parts);
    lower_extension_length = extension_length - 3*std(extension_parts);

    scale = 2;

    output = cellfun(@(x) setfield(x, 'Flexion_Length', flexion_length*scale), ...
                output, 'UniformOutput', false);
    output = cellfun(@(x) setfield(x, 'Extension_Length', extension_length*scale), ...
                output, 'UniformOutput', false);
        
    % output{1, j}.Flexion_Length = flexion_length*scale;
    % output{1, j}.Extension_Length = extension_length*scale;
        
    if ismember(subject, [11, 12, 15, 16, 17, 18])
        
        Forces = cellfun(@(x) x.Force_sensor_raw, Knee_Torque_Force_Sensor, ...
            'UniformOutput', false); 
        Angles = cellfun(@(x) x.Knee_Angle_raw, Knee_Torque_Force_Sensor, ...
            'UniformOutput', false);
    
        for j = 1:length(Forces)
    
            if isempty(Knee_Torque_Force_Sensor{1, j}.Force_sensor_raw)
                continue
            end
            
            a = Trials_Info{1, j}.Events.EXP_stream.flextoflex_start_indx;
            b = Trials_Info{1, j}.Events.EXP_stream.extension_start_indx;
            b = b(1:length(a));
            c = Trials_Info{1, j}.Events.EXP_stream.flextoflex_end_indx;
    
            constraint1 = and(b-a+1 < upper_flexion_length, ...
                b-a+1 > lower_flexion_length);
            constraint2 = and(c-b < upper_extension_length, ...
                c-b > lower_extension_length);
            constraint = and(constraint1, constraint2);
    
            F = Forces{1, j}(1, constraint);
            D = Angles{1, j}(1, constraint);
            a = num2cell(a(constraint)); 
            b = num2cell(b(constraint)); 
            c = num2cell(c(constraint)); 
    
            Times_trial = Times{1, j};
            Times_trial = Times_trial(constraint);
    
            F_timeWarped = cellfun(@(x1, x2, x3, x4, x5) ...
                [interp1(x5(1:x2-x1+1), x4(1:x2-x1+1), linspace(x5(1), x5(x2-x1+1), flexion_length*scale), 'linear'), ...
                 interp1(x5(x2-x1+2:x3-x1+1), x4(x2-x1+2:x3-x1+1), linspace(x5(x2-x1+2), x5(x3-x1+1), extension_length*scale), 'linear')], ...
                a, b, c, F, Times_trial, 'UniformOutput', false);
            F_timeWarped = squeeze(cat(3, F_timeWarped{:}))';
            output{1, j}.Force_sensor_TimeWarped = F_timeWarped;

            D_timeWarped = cellfun(@(x1, x2, x3, x4, x5) ...
                [interp1(x5(1:x2-x1+1), x4(1:x2-x1+1), linspace(x5(1), x5(x2-x1+1), flexion_length*scale), 'linear'), ...
                 interp1(x5(x2-x1+2:x3-x1+1), x4(x2-x1+2:x3-x1+1), linspace(x5(x2-x1+2), x5(x3-x1+1), extension_length*scale), 'linear')], ...
                a, b, c, D, Times_trial, 'UniformOutput', false);
            D_timeWarped = squeeze(cat(3, D_timeWarped{:}))';
            output{1, j}.Knee_Angle_TimeWarped = D_timeWarped;
    
    
        end
    end



    %% Time Warping the knee torque data

    Torques0 = cellfun(@(x) x.Torque_raw, Knee_Torque_Force_Sensor, ...
        'UniformOutput', false); 
    Torques  = cell(size(Torques0));
    for j = 1:length(Torques)
        Torques{1, j}  = cellfun(@(x) x(:, 2), Torques0{1, j}, ...
            'UniformOutput', false); 
    end


    a_prim = []; b_prim = []; c_prim = [];
    for j = 1:length(Torques)

        if isempty(Knee_Torque_Force_Sensor{1, j}.Torque_raw)
            continue
        end

        a = Trials_Info{1, j}.Events.EXP_stream.flextoflex_start_indx;
        b = Trials_Info{1, j}.Events.EXP_stream.extension_start_indx;
        b = b(1:length(a));
        
        for k = 1:length(Torques{1, j})
            flex2flex_start_time = Times{1, j}{1, k}(1);
            extension_start_time = Times{1, j}{1, k}(b(k) - a(k) + 1);
            flex2flex_end_time   = Times{1, j}{1, k}(end);

            [~, a_T] = min(abs(Knee_Torque_Force_Sensor{1, j}.Torque_raw{1, k}(:, 1) - flex2flex_start_time ));
            [~, b_T] = min(abs(Knee_Torque_Force_Sensor{1, j}.Torque_raw{1, k}(:, 1) - extension_start_time ));
            [~, c_T] = min(abs(Knee_Torque_Force_Sensor{1, j}.Torque_raw{1, k}(:, 1) - flex2flex_end_time ));

            a_prim = cat(1, a_prim, a_T);
            b_prim = cat(1, b_prim, b_T);
            c_prim = cat(1, c_prim, c_T);
        end

    end

    flexion_parts_prim = b_prim - a_prim + 1;
    extension_parts_prim = c_prim - b_prim;

    flexion_length_prim = floor(median(flexion_parts_prim));
    extension_length_prim = floor(median(extension_parts_prim));
    upper_flexion_length_prim = flexion_length_prim + 3*std(flexion_parts_prim);
    lower_flexion_length_prim = flexion_length_prim - 3*std(flexion_parts_prim);
    upper_extension_length_prim = extension_length_prim + 3*std(extension_parts_prim);
    lower_extension_length_prim = extension_length_prim - 3*std(extension_parts_prim);



    for j = 1:length(Torques)

        if isempty(Knee_Torque_Force_Sensor{1, j}.Torque_raw)
            continue
        end

        % indexes of epochs which have positive knee angle
        Knee_angles = data{1, j}.EXP_stream.Encoder_angle;
        %%% remove the epochs with positive knee angles
        Knee_angles_negativity = true(1, size(Knee_angles, 2));
        % for k = 1:length(Knee_angles)
        %     if any(Knee_angles{1, k} > 0)
        %         Knee_angles_negativity(k) = false;
        %     end
        % end
        accepted_indx = 1:length(Knee_angles_negativity);
        accepted_indx = accepted_indx(Knee_angles_negativity);


        a = Trials_Info{1, j}.Events.EXP_stream.flextoflex_start_indx;
        b = Trials_Info{1, j}.Events.EXP_stream.extension_start_indx;
        b = b(1:length(a));
        
        a_prim = []; b_prim = []; c_prim = [];
        for k = 1:length(Torques{1, j})
            flex2flex_start_time = Times{1, j}{1, accepted_indx(k)}(1);
            extension_start_time = Times{1, j}{1, accepted_indx(k)}(b(accepted_indx(k)) - a(accepted_indx(k)) + 1);
            flex2flex_end_time   = Times{1, j}{1, accepted_indx(k)}(end);

            [~, a_T] = min(abs(Knee_Torque_Force_Sensor{1, j}.Torque_raw{1, k}(:, 1) - flex2flex_start_time ));
            [~, b_T] = min(abs(Knee_Torque_Force_Sensor{1, j}.Torque_raw{1, k}(:, 1) - extension_start_time ));
            [~, c_T] = min(abs(Knee_Torque_Force_Sensor{1, j}.Torque_raw{1, k}(:, 1) - flex2flex_end_time ));

            a_prim = cat(1, a_prim, a_T);
            b_prim = cat(1, b_prim, b_T);
            c_prim = cat(1, c_prim, c_T);
        end

        constraint1 = and(b_prim-a_prim+1 < upper_flexion_length_prim, ...
            b_prim-a_prim+1 > lower_flexion_length_prim);
        constraint2 = and(c_prim-b_prim < upper_extension_length_prim, ...
            c_prim-b_prim > lower_extension_length_prim);
        constraint = and(constraint1, constraint2);


        T = Torques{1, j}(1, constraint);
        a_prim = num2cell(a_prim(constraint)); 
        b_prim = num2cell(b_prim(constraint)); 
        c_prim = num2cell(c_prim(constraint)); 

        Times_trial = cellfun(@(x) x(:,1), ...
            Knee_Torque_Force_Sensor{1, j}.Torque_raw, 'UniformOutput', false);
        Times_trial = Times_trial(constraint);

        T_timeWarped = cellfun(@(x1, x2, x3, x4, x5) ...
            [interp1(x5(x1:x2), x4(x1:x2), linspace(x5(x1), x5(x2), flexion_length*scale), 'spline'), ...
             interp1(x5(x2+1:x3), x4(x2+1:x3), linspace(x5(x2+1), x5(x3), extension_length*scale), 'spline')], ...
            a_prim', b_prim', c_prim', T, Times_trial, 'UniformOutput', false);
        T_timeWarped = squeeze(cat(3, T_timeWarped{:}))';
        output{1, j}.Torque_TimeWarped = T_timeWarped;

    end



end