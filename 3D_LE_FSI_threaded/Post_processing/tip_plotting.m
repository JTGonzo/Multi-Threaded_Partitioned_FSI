clear all; close all; clc;

dt = 0.001;
fontSize = 48;

%% Load displacement data
fileID1 = fopen('../Cases/1/Results/Flap_Disp.othd','r');
formatSpec = '%f';
data1 = fscanf(fileID1,formatSpec);

ts = data1(1:4:end-3); 
dx = data1(2:4:end-2);
dy = data(3:4:end-1);


dz = data(4:4:end);
magD = sqrt(dx.^2+dy.^2);

%% Load velocity data
fileID2 = fopen('../Cases/1/Results/Flap_Vel.othd','r');
formatSpec = '%f';
data2 = fscanf(fileID2,formatSpec);

tf = data2(1:3:end-2); 
vx = data2(2:3:end-1);
vy = data2(3:3:end);
magV = sqrt(vx.^2+vy.^2);

%% Displacement Plots
figure(1)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(ts, dx,'-b','LineWidth',4);
plot(ts(1:40:end), dx(1:40:end),'ob','MarkerSize',10,'LineWidth', 1);
xlabel('$tU_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$A_{x}/L$','Fontsize',fontSize,'Interpreter','Latex');

figure(2)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(ts, dy,'-r','LineWidth',4);
plot(ts(1:40:end), dy(1:40:end),'*r','MarkerSize',10,'LineWidth', 1);
xlabel('$tU_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$A_{y}/L$','Fontsize',fontSize,'Interpreter','Latex');

figure(3)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(ts, magD,'-g','LineWidth',4);
plot(ts(1:40:end), magD(1:40:end),'^g','MarkerSize',10,'LineWidth', 1);
xlabel('$tU_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$||\vec{D}||$','Fontsize',fontSize,'Interpreter','Latex');

%% Velocity Plots
figure(4)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(tf, vx,'-b','LineWidth',4);
plot(tf(1:50:end), vx(1:50:end),'ob','MarkerSize',10,'LineWidth', 1);
xlabel('$tU_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$V_{x}$','Fontsize',fontSize,'Interpreter','Latex');

figure(5)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(tf, vy,'-r','LineWidth',4);
plot(tf(1:50:end), vy(1:50:end),'*r','MarkerSize',10,'LineWidth', 1);
xlabel('$tU_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$V_{y}$','Fontsize',fontSize,'Interpreter','Latex');

figure(6)
hold on
grid off
box off
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(tf, magV,'-g','LineWidth',4);
plot(tf(1:50:end), magV(1:50:end),'^g','MarkerSize',10,'LineWidth', 1);
xlabel('$tU_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$||\vec{V}||$','Fontsize',fontSize,'Interpreter','Latex');
