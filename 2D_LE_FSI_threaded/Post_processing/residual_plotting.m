clear all; close all; clc;

dt = 0.001;
fontSize = 48;

%% Load residual data
fileID1 = fopen('Disp_Residual.txt','r');
formatSpec = '%f';
data1 = fscanf(fileID1,formatSpec);

t = data1(1:3:end-2); 
rx = data1(2:3:end-1);
rt = data1(3:3:end);

%% Residual Displacement Plots
figure(1)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32);
plot(t, rx,'-ob','LineWidth',4);
xlabel('FP Iterations','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$||r_{dx}||$','Fontsize',fontSize,'Interpreter','Latex');
%legend({'$\rho* = 1$','$\rho* = 10$','$\rho* = 100$','$\rho* = 500$','$\rho* = 1000$'},'interpreter','latex','Location','northeast','FontSize',30);

%% Residual Traction Plots
figure(2)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32);
plot(t, rt,'-dr','LineWidth',4);
xlabel('FP Iterations','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$||r_{t}||$','Fontsize',fontSize,'Interpreter','Latex');
%legend({'$\rho* = 1$','$\rho* = 10$','$\rho* = 100$','$\rho* = 500$','$\rho* = 1000$'},'interpreter','latex','Location','northeast','FontSize',30);
