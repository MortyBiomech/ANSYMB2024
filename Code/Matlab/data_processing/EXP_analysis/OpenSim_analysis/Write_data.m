clc; clear; close all
%%

load("Epochs_FlextoFlex_based.mat")
cond = 24;epoch = 5; sub = 18;
Time = Epochs_FlextoFlex_based{1,cond}.EXP_stream.Times{1,epoch};
Knee_ang = Epochs_FlextoFlex_based{1,cond}.EXP_stream.Encoder_angle{1,epoch};

IK = read_motionFile('IK.mot');

data.labels = IK.labels;
data.data(:,1)  = Time;
data.data(:,11) = -Knee_ang + 5;

for i=1:length(data.labels)
    if i==1 || i==11
    elseif i==8 || i==16
        data.data(:,i) = 90*ones(size(Time));
    else
        data.data(:,i) = zeros(size(Time));
    end
end

fname = ['Sub',num2str(sub),'Cond',num2str(cond),'Epoch',num2str(epoch),'IK.mot'];
writeGRFsToMOT(data,fname)


% ID4 = read_motionFile('ID10\inverse_dynamics_Results.sto');
% ID6 = read_motionFile('ID6\inverse_dynamics.sto');
% ID2 = read_motionFile('ID2\inverse_dynamics.sto');
% 
% figure
% plot(ID2.data(:,17))
% hold on
% plot(ID4.data(:,17))
% plot(ID6.data(:,17))





