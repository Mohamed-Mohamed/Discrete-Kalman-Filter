% This code  is used to test DKF function
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
%% Plotting Control
% Plot = 0 --> No Plotting data    Plot = 1 --> Plotting data
Plot=1;
%% Save figures control
% Save = 0 --> No Saving figures    Save = 1 --> Saving figures
% The figures will saved in running folder directory
Save=1;
%% Time
Time=linspace(0,40,40);
dt=Time(2)-Time(1);
%% initial conditions
x(1,1)=[1];
%% state-transition matrix 
phi=1;
%% variance of measurement error matrix
R=100;
%% variance of Process model noise
Q=1;
%% Process Model states valuses
m=3;
M(1)=normrnd(0, sqrt(Q));
for i=2:length(Time)
    M(i)=normrnd(0, sqrt(Q));
    x(1,i)=phi*x(1,i-1)+dt+M(i);
end
%% measurements
z(1)=0;
H=1;
L(1)=normrnd(0, sqrt(R));
for i=2:length(Time)
    L(i)=normrnd(0, sqrt(R));
    z(i)=H*x(1,i)+L(i);
end
%% DKF parameter
x_hat0_min=x(1,1);
P0_min=zeros(1,1);
%% DKF function
[ x_estimated ] = DKF ( x, Time, phi, H, z, x_hat0_min, P0_min, R, Q, [Plot 1:5] );
%% DKFwithForgettingFactor
lambda=0.9; %   0 < lambda < 2
[ x_estimated1 ] = DKFwithForgettingFactor ( x, Time, phi, H, z, x_hat0_min, P0_min, R, Q, lambda, [Plot 6:10] );
%% DKFwithVaryingForgettingFactor
lambda_s=0.9; %   0.3 < lambda < 0.95
[ x_estimated2 ] = DKFwithVaryingForgettingFactor ( x, Time, phi, H, z, x_hat0_min, P0_min, R, Q, lambda_s, [Plot 11:16] );
%% Save figures
if Save ==1 && Plot ==1
    for S=1:16
        figure(S);
        saveas(gcf, [num2str(S) '.png']);
    end
end