function [ x_estimated ] = DKFwithForgettingFactor ( x, Time, phi, H, Z, x_hat0_min, P0_min, R, Q, lambda, Plot )
%{ 
      This function is used to estimate the the true state from noise measurements by 
      Discret Kalman filter (DKF) with constant values of variance of measurement error matrix
      and variance of Process model noise matrix and constant Forgetting Factor
      
       Process Model :
       x(k+1) = phi(k) * x(k) + w(k)
       
       Measurement noise
       z(k) = H(k) * x(k) + v(k)

        Where :
         x(k) = n x 1 state vector at t(k)
         phi(k) = n x m state-transition matrix at t(k)
         z(k) = m x 1 measurement vector at t(k)
         H(k) = m x n measurement matrix at t(k)
%} 

%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs
% x                       : matrix of Process Model states 
% Time                : Time vector of simulation
% phi                   : state-transition matrix 
% H                      : measurement matrix 
% Z                      : measurement matrix at t(k)
% x_hat0_min      : initial conditions
% P0_min             :initial error covariance matrix
% R                       : variance of measurement error matrix
% Q                      : variance of Process model noise
% lambda            : Forgetting Factor  0 < lambda < 2
% Plot          : used to get the plot of system parameter estimation
%                   1 - if Plot = 0  --> no plot needed
%                   2- if Plot = 1  -->  plot needed and Plot(2:6) == figuires number
%% outputs
% x_estimated : the estimated states by DKF

%% function body
        %% DKF
        P_min{1}=P0_min;
        x_hat_min{1}=x_hat0_min;
        for i=1:length(Z)
            %% Calculate the gain
            K{i}=P_min{i}*H'*inv(H*P_min{i}*H'+R*lambda);
            %% Update estimate
            x_hat_p{i}=x_hat_min{i}+K{i}*(Z(i)-H*x_hat_min{i});
            %% Update error
            P_p{i}=(eye(length(H(1,:)))-K{i}*H)*P_min{i}/lambda;
            %% Project ahead
            x_hat_min{i+1}=phi*x_hat_p{i};
            P_min{i+1}=phi*P_p{i}*phi'+Q;
            P_min{i+1}=(P_min{i+1}+P_min{i+1}')/2;
        end
        x_estimated=cell2mat(x_hat_p);
%% Plotting
if Plot(1) == 1
    set(0,'defaultfigureposition',[0 50 1700 630]);
    %% State plot
    figure(Plot(2));
    set(gcf,'color','w');
    for i=1:length(x(:,1))
        hold all;
        subplot(length(x(:,1)),1,i);
        plot(Time,x(i,:));
        grid on;
        ylabel(['x'  num2str(i) ' (t)'],'fontsize',18)
        legend(['x' num2str(i)]);
    end
    subplot(length(x(:,1)),1,i);
    xlabel('Time (sec)','fontsize',18)
    subplot(length(x(:,1)),1,1);
    title('States Responses','fontsize',18)
    %% State and measurment plot
    figure(Plot(3));
    set(gcf,'color','w');
    for i=1:length(Z(:,2))
        hold all;
        if length(Z(:,2)) >1
            subplot(length(Z(:,2)),1,i);
        end
        plot(Time,x(i,:));
        plot(Time(2:end),Z(i,2:end),'--');
        grid on;
        ylabel(['x'  num2str(i) ' (t), z'  num2str(i) ' (t)'],'fontsize',18)
        legend({['x' num2str(i)],['z' num2str(i)] });
    end
    if i > 1
        subplot(length(x(:,1)),1,i);
        xlabel('Time (sec)','fontsize',18)
        subplot(length(x(:,1)),1,1);
        title('States and Measurments Responses','fontsize',18)
    else
        xlabel('Time (sec)','fontsize',18)
        title('States and Measurments Responses','fontsize',18)
    end
    %% Abs. error between states and measurments
    figure(Plot(4));
    set(gcf,'color','w');
    for i=1:length(Z(:,2))
        hold all;
        if length(Z(:,2)) >1
            subplot(length(Z(:,2)),1,i);
        end
        plot(Time(2:end),100*abs((Z(i,2:end)-x(i,2:end)))./x(i,2:end));
        grid on;
        ylabel(['error'  num2str(i) '%'],'fontsize',18)
        legend({['error'  num2str(i) '% with \mu = ' num2str(mean(100*abs((Z(i,2:end)-x(i,2:end))./x(i,2:end)))) ', \sigma = ' num2str(std(100*abs((Z(i,2:end)-x(i,2:end))./x(i,2:end))))] });
    end
    if i > 1
        subplot(length(x(:,1)),1,i);
        xlabel('Time (sec)','fontsize',18,'Interpreter','latex')
        subplot(length(x(:,1)),1,1);
        title('error% between states and measurments','fontsize',18,'Interpreter','latex')
    else
        xlabel('Time (sec)','fontsize',18)
        title('error% between states and measurments','fontsize',18)
    end
    %% State, measurment and filtered data plot 
    figure(Plot(5));
    set(gcf,'color','w');
    for i=1:length(Z(:,2))
        hold all;
        if length(Z(:,2)) >1
            subplot(length(Z(:,2)),1,i);
        end
        plot(Time,x(i,:));
        plot(Time(2:end),Z(i,2:end),'--');
        plot(Time(2:end),x_estimated(i,2:end),'-.');
        grid on;
        ylabel(['x'  num2str(i) ' (t), z'  num2str(i) ' (t), x_e_s_t' num2str(i) ' (t)'],'fontsize',18)
        legend({['x' num2str(i)],['z' num2str(i)],['x_e_s_t' num2str(i)] });
    end
    if i > 1
        subplot(length(x(:,1)),1,i);
        xlabel('Time (sec)','fontsize',18)
        subplot(length(x(:,1)),1,1);
        title('States, Measurments and and Filtered data Responses','fontsize',18)
    else
        xlabel('Time (sec)','fontsize',18)
        title('States, Measurments and and Filtered data Responses','fontsize',18)
    end
    %% Abs. error between states and filtered data
    figure(Plot(6));
    set(gcf,'color','w');
    for i=1:length(x(:,2))
        hold all;
        if length(x(:,2)) >1
            subplot(length(x(:,2)),1,i);
        end
        plot(Time(2:end),100*abs((x_estimated(i,2:end)-x(i,2:end))./x(i,2:end)));
        grid on;
        ylabel(['error' num2str(i) '%'],'fontsize',18)
        legend({['error'  num2str(i) '% with \mu = ' num2str(mean(100*abs((x_estimated(i,5:end)-x(i,5:end))./x(i,5:end)))) ', \sigma = '  num2str(std(100*abs((x_estimated(i,5:end)-x(i,5:end))./x(i,5:end))))] });
    end
    if i > 1
        subplot(length(x(:,1)),1,i);
        xlabel('Time (sec)','fontsize',18)
        subplot(length(x(:,1)),1,1);
        title('error% between states and filtered data','fontsize',18)
    else
        xlabel('Time (sec)','fontsize',18)
        title('error% between states and filtered data','fontsize',18)
    end
end
end