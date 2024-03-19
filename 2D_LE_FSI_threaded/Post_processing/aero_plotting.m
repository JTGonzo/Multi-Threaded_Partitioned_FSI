clear all; close all; clc;

dt = 0.001;
fontSize = 48;

%% Load aerodynamic data
fileID = fopen('../Cases/1/Results/AeroForces.txt','r');
formatSpec = '%f';
data = fscanf(fileID,formatSpec);

t = data(1:4:end-3); 
Fx = data(2:4:end-2);
Fy = data(3:4:end-1);
Fz = data(4:4:end);
mag = sqrt(Fx.^2+Fy.^2+Fz.^2);

%% Lift Plots
figure(1)
hold on
grid off
box off
%set(gca,'TickDir','out','FontSize',32,'xlim',[0 12],'xtick',0:4:12,'ylim',[-0.03 0.005],'ytick',-0.03:0.01:0);
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(t, Fx,'-b','LineWidth',2);
xlabel('$t U_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$C_{D}$','Fontsize',fontSize,'Interpreter','Latex');

ax=axes;
set(ax,'units','normalized','position',[.35 .47 .41 .41])
box(ax,'on')
plot(t, Fx,':b','LineWidth',1,'parent',ax);
hold on
%grid on
grid minor
plot(t(1:15:end), Fx(1:15:end),'ob','MarkerSize',4,'LineWidth', 2,'parent',ax);
%set(ax,'FontSize',20,'xlim',[7.88,7.96],'ylim',[-0.029,-0.0265],'xtick',7.88:0.04:7.96,'ytick',-0.028:0.001:-0.0265);
set(ax,'FontSize',20,'xlim',[0,2],'ylim',[-50,200],'ytick',-50:100:200,'xtick',0:1:2);

%% Drag Plots
figure(2)
hold on
grid off
box off
%set(gca,'TickDir','out','FontSize',32,'xlim',[0 12],'xtick',0:4:12,'ylim',[-0.03 0.005],'ytick',-0.03:0.01:0);
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(t, Fy,'-r','LineWidth',2);
xlabel('$t U_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$C_{L}$','Fontsize',fontSize,'Interpreter','Latex');

ax=axes;
set(ax,'units','normalized','position',[.35 .3 .41 .41])
box(ax,'on')
plot(t, Fy,':r','LineWidth',1,'parent',ax);
hold on
%grid on
grid minor
plot(t(1:15:end), Fy(1:15:end),'or','MarkerSize',4,'LineWidth', 2,'parent',ax);
%set(ax,'FontSize',20,'xlim',[7.88,7.96],'ylim',[-0.029,-0.0265],'xtick',7.88:0.04:7.96,'ytick',-0.028:0.001:-0.0265);
set(ax,'FontSize',20,'xlim',[0,2],'ylim',[-50,50],'ytick',-50:50:50,'xtick',0:1:2);

%% Twist Plots
figure(3)
hold on
grid off
box off
%set(gca,'TickDir','out','FontSize',32,'xlim',[0 12],'xtick',0:4:12,'ylim',[-0.03 0.005],'ytick',-0.03:0.01:0);
set(gca,'TickDir','out','FontSize',32,'xlim',[0 5]);
plot(t, Fz,'-m','LineWidth',2);
xlabel('$t U_0/L$','Fontsize',fontSize,'Interpreter','Latex');
ylabel('$C_{Z}$','Fontsize',fontSize,'Interpreter','Latex');


