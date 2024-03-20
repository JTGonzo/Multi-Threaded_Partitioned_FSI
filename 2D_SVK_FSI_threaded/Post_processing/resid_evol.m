clear all; close all; clc; fclose('all');

fontsize = 48;
formatSpec = '%f';

%% Load residual data
addpath('../Cases/Case_1/Results/')
fileID1 = fopen('HC_2DL1T1Z1R2_residuals.txt','r');
data1 = fscanf(fileID1,formatSpec);
ts1 = data1(1:7:end-6); 
ite1 = data1(2:7:end-5);
rs1 = data1(3:7:end-4);
rd1 = data1(4:7:end-3);
rv1 = data1(5:7:end-2);
rp1 = data1(6:7:end-1);
rt1 = data1(7:7:end);

ind1 = find(ite1(:) == 1);
ind2 = find(ite1(:) == 2);
ind3 = find(ite1(:) == 3);
ind4 = find(ite1(:) == 4);
ind5 = find(ite1(:) == 5);

%% Residual Norm of Interface Displacement
figure(1)
px1 = semilogy(ts1(ind1),rd1(ind1),'-b','LineWidth',2);
grid on
hold on
px2 = semilogy(ts1(ind2),rd1(ind2),'-r','LineWidth',2);
px3 = semilogy(ts1(ind3),rd1(ind3),'-g','LineWidth',2);
px4 = semilogy(ts1(ind4),rd1(ind4),'-m','LineWidth',2);
px5 = semilogy(ts1(ind5),rd1(ind5),'-k','LineWidth',2);
legend([px1 px2 px3 px4 px5],{'$iter_1$','$iter_2$','$iter_3$','$iter_4$','$iter_5$'},'Interpreter','Latex','Location','southwest','Fontsize',24)
box off
set(gca,'TickDir','out','fontsize',32);
xlabel('Time step','fontsize',fontsize,'Interpreter','Latex');
ylabel('$||r_{u^s}||_2$','fontsize',fontsize,'Interpreter','Latex');


%%  Residual Norm of Interface Tractions
figure(2)
pt1 = semilogy(ts1(ind1),rt1(ind1),'-b','LineWidth',2);
grid on
hold on
pt2 = semilogy(ts1(ind2),rt1(ind2),'-r','LineWidth',2);
pt3 = semilogy(ts1(ind3),rt1(ind3),'-g','LineWidth',2);
pt4 = semilogy(ts1(ind4),rt1(ind4),'-m','LineWidth',2);
pt5 = semilogy(ts1(ind5),rt1(ind5),'-k','LineWidth',2);
legend([pt1 pt2 pt3 pt4 pt5],{'$iter_1$','$iter_2$','$iter_3$','$iter_4$','$iter_5$'},'Interpreter','Latex','Location','southwest','Fontsize',24)
box off
set(gca,'TickDir','out','fontsize',32);
xlabel('Time step','fontsize',fontsize,'Interpreter','Latex');
ylabel('$||r_{t^f}||$','Fontsize',fontsize,'Interpreter','Latex');