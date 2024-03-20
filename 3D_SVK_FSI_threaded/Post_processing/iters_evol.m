clear all; close all; clc; fclose('all'); 

%% Format plot
fontsize = 48;
formats1pec = '%f';

%% Load data
addpath('../Results/')
fileID1 = fopen('TB_3DL1T1Z1_iters.oisd','r');
data2 = fscanf(fileID1,formats1pec);
ts2 = data2(1:2:end-3); 
iter2 = data2(2:2:end-2);

%%
figure(2)
px2 = plot(ts2(1:1:end), iter2(1:1:end),'+r','LineWidth',1,'MarkerSize',8);
grid on
hold on
box off
set(gca,'TickDir','out','fontsize',32)
xlabel('Time step','fontsize',fontsize,'Interpreter','Latex');
ylabel('\# Iterations','fontsize',fontsize,'Interpreter','Latex');