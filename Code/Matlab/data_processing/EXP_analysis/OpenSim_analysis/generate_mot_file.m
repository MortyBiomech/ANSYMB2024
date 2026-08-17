function generate_mot_file(labels, Time, Knee_ang, sub, trial, epoch)
    data.labels = labels;
    data.data(:,1)  = Time;
    % data.data(:, 11) = -Knee_ang;
    min_knee_angle = min(-Knee_ang);
    if min_knee_angle < 0
        data.data(:, 11) = -Knee_ang - min_knee_angle;
    else
        data.data(:, 11) = -Knee_ang;
    end
    
    for i=1:length(data.labels)
        if i==1 || i==11
        elseif i==8 || i==16
            data.data(:,i) = 90*ones(size(Time));
        else
            data.data(:,i) = zeros(size(Time));
        end
    end
    
    % disp(['Sub-',num2str(sub),'_Trial',num2str(trial),'_Epoch',num2str(epoch),'_IK.mot'])
    fname = ['Sub-',num2str(sub),'_Trial',num2str(trial),'_Epoch',num2str(epoch),'_IK.mot'];
    writeGRFsToMOT(data,fname)
end