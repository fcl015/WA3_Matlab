clear all;

% Input parameters
Q_cost=0;  % 30 lts during 1 sample time = 1 then multiply by weight 0.1
R_cost=1;   % 3% error during 1 sample time = 1  then multiply by weight 1

low_limit=27.0;
high_limit=28.5;
sm_initial=29.0;
set_point=29.0;

ndvi_initial=0.8;
ndvi_limit=0.8;
ndvi_irr=60;

previous_eto=[];
horizon=30;
cost_J=0;

% Obtain previous eto (last 24 hours)
dataset=5;
offset=1000;
[previous_eto]=plant_dynamics_previous_eto(dataset,offset);

% Main program
verbose=1;
figure(1)
cost_sm = simulate_irrigation_sm(Q_cost,R_cost,set_point,low_limit,high_limit,previous_eto,sm_initial,ndvi_initial,horizon,verbose); 
 
figure(2)
verbose=1;
cost_ndvi = simulate_irrigation_ndvi(Q_cost,R_cost,set_point,low_limit,high_limit,previous_eto,sm_initial,ndvi_initial,ndvi_limit,ndvi_irr,horizon,verbose); 
 


