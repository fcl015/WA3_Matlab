% Data for validation for 5 dataset
clear all

sm_coef_above=   0.99900;
sm_coef_middle=  1.00000;
sm_coef_below=   1.00000;
ir_coef_above=   0.00500;
ir_coef_middle=  0.00550;
ir_coef_below=   0.00200;
eto_coef_above= -0.00010;
eto_coef_middle=-0.00040;
eto_coef_below= -0.00025;

display_date=0;
display_plot=1;
R2_matrix=[];

display('=================================');

for dataset=[1 2 3 4 5]
    switch dataset
        case {1}
            load('DataSet_01.mat');
            initial_date=datenum(2014,6,19,12,08,0);
            data_len=length(sm);
        case {2}
            load('DataSet_02.mat');
            initial_date=datenum(2014,7,10,0,0,0);
            data_len=length(sm);
        case {3}
            load('DataSet_03.mat');
            initial_date=datenum(2014,8,6,10,08,0);
            data_len=length(sm);
        case {4}
            load('DataSet_04.mat');
            initial_date=datenum(2014,9,18,0,0,0);
            data_len=length(sm);
        case {5}
            load('DataSet_05.mat');
            initial_date=datenum(2014,10,10,0,0,0);
            data_len=length(sm);
        otherwise
            disp('Not valid dataset')
    end

    data_len=length(sm);
    if display_date
        display('=================================');
        display('Initial Date');
        display(datestr(initial_date));
    end

    %Create a time vector based on sm size
    sim_date=[];
    sim_date(1,1)=initial_date;
    for i=2:1:data_len
        sim_date(i,1)=addtodate(initial_date,i-1,'minute');
    end

    % Obtain final date
    final_date=sim_date(data_len);
    if display_date
        display('Final Date');
        display(datestr(final_date));
        display('=================================');
    end

    % Use filter kalman for ETo
    eto_kalman=[];
    eto_kalman_ant=eto(1);
    p_ant=1;
    Q=0.0001;
    R=0.1;
    for i=1:1:data_len
            K=(p_ant+Q)/(p_ant+Q+R);
            eto_kalman(i,1)=eto_kalman_ant+K*(eto(i,1)-eto_kalman_ant);
            p=(p_ant+Q)*(1-K);
            p_ant=p;
            eto_kalman_ant=eto_kalman(i,1);        
    end

    % Verify mathematical model SM
    xk_ant=sm(1);
    sm_est=[];
    t=[];
    for i=1:1:length(sm)
        if (i<=20)
            tau=i-1;
        else
            tau=20;
        end
        if xk_ant > 39.0
            xk=sm_coef_above*xk_ant+ir_coef_above*ir(i-tau)+eto_coef_above*eto_kalman(i-tau);
        elseif xk_ant > 31.0
            xk=sm_coef_middle*xk_ant+ir_coef_middle*ir(i-tau)+eto_coef_middle*eto_kalman(i-tau);
        else
            xk=sm_coef_below*xk_ant+ir_coef_below*ir(i-tau)+eto_coef_below*eto_kalman(i-tau);
        end   
        if( xk>50.0)
            xk=50.0;
        end  
        sm_est=[sm_est xk];
        t=[t i];
        xk_ant=xk;
    end

    % Obtain estimated NDVI
    ndvi_est=[];
    c4=0.999;
    c5=0.00374;
    initial_ndvi=0.80;
    tauB=1440;
    t=[];
    for i=1:1:length(sm)
        if (i<=2*tauB)
            current_ndvi=initial_ndvi;
        else       
            current_ndvi=c4*ndvi_est(i-tauB)+c5*(sm(i-tauB)-sm(i-2*tauB));
        end
        if( current_ndvi>1.0)
            current_ndvi=1.0;
        end
        if( current_ndvi<0.0)
            current_ndvi=0.0;
        end
        ndvi_est=[ndvi_est current_ndvi];
        t=[t i];
    end

    % Display plots
    if display_plot
        subplot(4,1,1)
        %hold on
        plot(sim_date,sm,'k-','LineWidth',2);
        %plot(sim_date,sm_est,'b--','LineWidth',2);
        datetick('x','DD(HH)','keeplimits')
        axis([min(sim_date(:,1)) max(sim_date(:,1)) 20 50]);
        %hold off
        
        subplot(4,1,2)
        plot(sim_date,ndvi_est,'b-','LineWidth',2);
        datetick('x','DD(HH)','keeplimits')
        axis([min(sim_date(:,1)) max(sim_date(:,1)) 0.5 1.0]);

        subplot(4,1,3)
        plot(sim_date,ir,'r-','LineWidth',2);
        datetick('x','DD(HH)','keeplimits')
        axis([min(sim_date(:,1)) max(sim_date(:,1)) 0 40]);
        
         subplot(4,1,4)
         %hold on
         %plot(sim_date,eto,'g-','LineWidth',2);
         plot(sim_date,eto_kalman,'b-','LineWidth',2);
         datetick('x','DD(HH)','keeplimits')
         %hold off
         axis([min(sim_date(:,1)) max(sim_date(:,1)) 0 20]);


    
    end

    % Calculate model error
    tmp=corrcoef(sm,sm_est);
    coef_R=tmp(1,2);
    %fprintf('R^2: %5.4f\n',coef_R);
    fprintf('%5.4f\n',coef_R);
    R2_matrix(dataset)=coef_R;
    
    if display_plot
        pause;
    end

end

%display('=================================');
%fprintf('Avg R^2: %5.4f\n',mean(R2_matrix));
fprintf('%5.4f\n',mean(R2_matrix));
%display('=================================');




