clear all; close all; clc; fclose('all'); 

fontSize = 48;
formatSpec = '%f';

%% Load data
addpath('../Results/')
%fileID = fopen('HC_2DL1T1Z1R2_forces.txt','r');
%fileID = fopen('HC_2DL1T1Z1R2_vel.othd','r');
fileID = fopen('HC_2DL1T1Z1R2_disp.othd','r');
data = fscanf(fileID,formatSpec);
t = data(1:3:end-2); 
dx = data(2:3:end-1);
dy = data(3:3:end);
mag = sqrt(dx.^2+dy.^2);

figure(1)
plot(t, dx,'-r','LineWidth',3);
hold on
grid off
plot(t(1:100:end), dx(1:100:end),'+r','MarkerSize',10,'LineWidth', 1);
box off
set(gca,'TickDir','out','FontSize',32);
xlabel('$t$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$\textbf{u}_{x}^{s}$','Fontsize',fontSize,'Interpreter','Latex');

figure(2)
plot(t, dy,'-b','LineWidth',3);
hold on
grid off
plot(t(1:100:end), dy(1:100:end),'*b','MarkerSize',10,'LineWidth', 1);
box off
set(gca,'TickDir','out','FontSize',32);
xlabel('$t$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$\textbf{u}_{y}^{s}$','Fontsize',fontSize,'Interpreter','Latex');

figure(3)
plot(t, mag,'-k','LineWidth',3);
hold on
grid off
plot(t(1:100:end), mag(1:100:end),'dk','MarkerSize',10,'LineWidth', 1);
box off
set(gca,'TickDir','out','FontSize',32);
xlabel('$t$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$||\textbf{u}^{s}||_2$','Fontsize',fontSize,'Interpreter','Latex');

