function [previous_eto] = plant_dynamics_previous_eto(dataset,offset) 

%Internal parameters
samplings_per_day=1440;

% Obtain previous eto (last 24 hours)
switch dataset
        case {1}
            load('DataSet_01.mat');
        case {2}
            load('DataSet_02.mat');
        case {3}
            load('DataSet_03.mat');
        case {4}
            load('DataSet_04.mat');
        case {5}
            load('DataSet_05.mat');
        otherwise
            disp('Not valid dataset')
end

% Use filter kalman for ETo
eto_kalman=[];
eto_kalman_ant=eto(1);
p_ant=1;
Q=0.0001;
R=0.1;
for i=1:1:length(eto)
    K=(p_ant+Q)/(p_ant+Q+R);
    eto_kalman(i,1)=eto_kalman_ant+K*(eto(i,1)-eto_kalman_ant);
    p=(p_ant+Q)*(1-K);
    p_ant=p;
    eto_kalman_ant=eto_kalman(i,1);
end
previous_eto=eto_kalman(offset+1:offset+samplings_per_day);



