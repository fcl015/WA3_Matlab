function [cost_J] = simulate_irrigation_sm(Q_cost,R_cost,set_point,low_limit,high_limit,previous_eto,x_0,ndvi_initial,horizon,verbose) 

% Model parameters
sm_coef_above=   0.99900;
sm_coef_middle=  1.00000;
sm_coef_below=   1.00000;
ir_coef_above=   0.00500;
ir_coef_middle=  0.00550;
ir_coef_below=   0.00200;
eto_coef_above= -0.00010;
eto_coef_middle=-0.00040;
eto_coef_below= -0.00025;
time_delay_value=   20.0;
above_value=        35.0;
below_value=        23.0;
saturation_value=   50.0;
samplings_per_day=  1440;
on_irrigation=      30.0;
off_irrigation=     0.0;
c4=                 0.999;
c5=                 0.00374;
tauB=               1440;


% Run simulation for N-days
xk_ant=x_0;
sm_est=[];
ir_est=[];
eto_est=[];
ndvi_est=[];
high_limit_vector=[];
low_limit_vector=[];
set_point_vector=[];
hist_cost=[];
t=[];
ir=0;
j=1;
cost_J=0;

for i=1:1:horizon*samplings_per_day
    % Update values for calculation
    sm_est=[sm_est xk_ant];
    ir_est=[ir_est ir];
    eto_est=[eto_est previous_eto(j)];
    j=j+1;
    if j>samplings_per_day
        j=1;
    end
    % Model delay
    if (i<=time_delay_value)
        tau=i-1;
    else
        tau=time_delay_value;
    end
    % Model plant dynamics
    if xk_ant > above_value
        xk=sm_coef_above*xk_ant+ir_coef_above*ir_est(i-tau)+eto_coef_above*eto_est(i-tau);
    elseif xk_ant > below_value
        xk=sm_coef_middle*xk_ant+ir_coef_middle*ir_est(i-tau)+eto_coef_middle*eto_est(i-tau);
    else
        xk=sm_coef_below*xk_ant+ir_coef_below*ir_est(i-tau)+eto_coef_below*eto_est(i-tau);
    end
    if xk>saturation_value
        xk=saturation_value;
    end
    % NDVI dynamics
    if (i<=2*tauB)
        current_ndvi=ndvi_initial;
    else       
        current_ndvi=c4*ndvi_est(i-tauB)+c5*(sm_est(i-tauB)-sm_est(i-2*tauB));
    end
    if( current_ndvi>1.0)
        current_ndvi=1.0;
    end
    if( current_ndvi<0.0)
        current_ndvi=0.0;
    end
    ndvi_est=[ndvi_est current_ndvi];
    % Conduct control action
    if xk>high_limit
        ir=off_irrigation;
    elseif xk<low_limit
        ir=on_irrigation;
    end
    % Calculate cost
    current_cost=Q_cost*(set_point-xk)^2+R_cost*(ir)^2;
    cost_J=cost_J+current_cost;
    hist_cost=[hist_cost cost_J];
    % Update values for next iteration
    set_point_vector=[set_point_vector set_point];
    high_limit_vector=[high_limit_vector high_limit];
    low_limit_vector=[low_limit_vector low_limit];
    t=[t i];
    xk_ant=xk;
end

% Display plots
if verbose
    subplot(3,1,1)
    hold on
    plot(t/1440,sm_est,'b-','LineWidth',2);
    plot(t/1440,set_point_vector,'k:','LineWidth',2);
    %plot(t,high_limit_vector,'k:','LineWidth',2);
    %plot(t,low_limit_vector,'k:','LineWidth',2);
    axis([min(t/1440) max(t/1440) 20 40]);
    hold off

    subplot(3,1,2)
    plot(t/1440,ndvi_est,'k-','LineWidth',2);
    axis([min(t/1440) max(t/1440) 0.5 1]);

    subplot(3,1,3)
    plot(t/1440,ir_est,'r-','LineWidth',2);
    axis([min(t/1440) max(t/1440) 0 40]);

    %subplot(4,1,4)
    %plot(t/1440,eto_est,'g-','LineWidth',2);
    %axis([min(t/1440) max(t/1440) 0 20]);

    %subplot(4,1,4)
    %plot(t,hist_cost,'r-','LineWidth',2);
    %axis([min(t) max(t) 0 max(hist_cost)]);

    display('=================================');
    fprintf('Soil Moisture Irrigation Cost: %5.2f\n',cost_J);
    display('=================================');
end




