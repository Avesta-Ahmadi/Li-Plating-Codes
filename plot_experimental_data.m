%%% -----------------------------------------------------------------
% Plot experimental data 
clear; clc; close all;
data = load('data.mat');
load('sdata.mat','sdata');


% C20 plots
figure(1);
hold on;
plot(data.R1.NMR.time, data.R1.NMR.total, 'DisplayName','Total Li', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.total_solid, 'DisplayName','Total Li - Solid Phase', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.phase2, 'DisplayName','Dilute Li', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.phase1, 'DisplayName','Disordered Li', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.electrolyte, 'DisplayName','Electrolyte Li', 'LineWidth',2)
legend('interpreter','latex','FontSize',15);
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Count','interpreter','latex','FontSize',15)
title('C20 Cycle')
hold off;
saveas(gcf,'C20.png')
fig = figure(10);
left_color = [0 0 0];
right_color = [0 0.4470 0.7410];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(data.R1.echem.time, data.R1.echem.V,'LineWidth',2)
ylabel('Voltage [V]','interpreter','latex','FontSize',15)
yyaxis right 
plot(data.R1.echem.time, data.R1.echem.I, '--','LineWidth',2)
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Current [A]','interpreter','latex','FontSize',15)
title('C20 Cycle')
saveas(gcf,'C20_DC.png')





% C3 plots
figure(2);
hold on;
plot(data.R2.NMR.time, data.R2.NMR.total, 'DisplayName','Total Li', 'LineWidth',2)
plot(data.R2.NMR.time, data.R2.NMR.total_solid, 'DisplayName','Total Li - Solid Phase', 'LineWidth',2)
plot(data.R2.NMR.time, data.R2.NMR.plated, 'DisplayName','Plated Li', 'LineWidth',2)
plot(data.R2.NMR.time, data.R2.NMR.phase3, 'DisplayName','Concentrated Li', 'LineWidth',2)
plot(data.R2.NMR.time, data.R2.NMR.phase2, 'DisplayName','Dilute Li', 'LineWidth',2)
plot(data.R2.NMR.time, data.R2.NMR.phase1, 'DisplayName','Disordered Li', 'LineWidth',2)
plot(data.R2.NMR.time, data.R2.NMR.electrolyte, 'DisplayName','Electrolyte Li', 'LineWidth',2)
legend('interpreter','latex','FontSize',15);
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Count','interpreter','latex','FontSize',15)
title('C3 Cycle')
hold off;
saveas(gcf,'C3.png')
fig = figure(20);
left_color = [0 0 0];
right_color = [0 0.4470 0.7410];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(data.R2.echem.time, data.R2.echem.V,'LineWidth',2)
ylabel('Voltage [V]','interpreter','latex','FontSize',15)
yyaxis right 
plot(data.R2.echem.time, data.R2.echem.I, '--','LineWidth',2)
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Current [A]','interpreter','latex','FontSize',15)
title('C3 Cycle')
saveas(gcf,'C3_DC.png')


% C2 plots
figure(3);
hold on;
plot(data.R3.NMR.time, data.R3.NMR.total, 'DisplayName','Total Li', 'LineWidth',2)
plot(data.R3.NMR.time, data.R3.NMR.total_solid, 'DisplayName','Total Li - Solid Phase', 'LineWidth',2)
plot(data.R3.NMR.time, data.R3.NMR.plated, 'DisplayName','Plated Li', 'LineWidth',2)
plot(data.R3.NMR.time, data.R3.NMR.phase3, 'DisplayName','Concentrated Li', 'LineWidth',2)
plot(data.R3.NMR.time, data.R3.NMR.phase2, 'DisplayName','Dilute Li', 'LineWidth',2)
plot(data.R3.NMR.time, data.R3.NMR.phase1, 'DisplayName','Disordered Li', 'LineWidth',2)
plot(data.R3.NMR.time, data.R3.NMR.electrolyte, 'DisplayName','Electrolyte Li', 'LineWidth',2)
legend('interpreter','latex','FontSize',15);
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Count','interpreter','latex','FontSize',15)
title('C2 Cycle')
hold off;
saveas(gcf,'C2.png')
fig = figure(30);
left_color = [0 0 0];
right_color = [0 0.4470 0.7410];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(data.R3.echem.time, data.R3.echem.V,'LineWidth',2)
ylabel('Voltage [V]','interpreter','latex','FontSize',15)
yyaxis right
plot(data.R3.echem.time, data.R3.echem.I,'--','LineWidth',2)
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Current [A]','interpreter','latex','FontSize',15)
title('C2 Cycle')
saveas(gcf,'C2_DC.png')


% 1C plots
figure(4);
hold on;
plot(data.R4.NMR.time, data.R4.NMR.total, 'DisplayName','Total Li', 'LineWidth',2)
plot(data.R4.NMR.time, data.R4.NMR.total_solid, 'DisplayName','Total Li - Solid Phase', 'LineWidth',2)
plot(data.R4.NMR.time, data.R4.NMR.dendrite, 'DisplayName','Dendrite Li', 'LineWidth',2)
plot(data.R4.NMR.time, data.R4.NMR.plated, 'DisplayName','Plated Li', 'LineWidth',2)
plot(data.R4.NMR.time, data.R4.NMR.phase3, 'DisplayName','Concentrated Li', 'LineWidth',2)
plot(data.R4.NMR.time, data.R4.NMR.phase2, 'DisplayName','Dilute Li', 'LineWidth',2)
plot(data.R4.NMR.time, data.R4.NMR.phase1, 'DisplayName','Disordered Li', 'LineWidth',2)
plot(data.R4.NMR.time, data.R4.NMR.electrolyte, 'DisplayName','Electrolyte Li', 'LineWidth',2)
legend('interpreter','latex','FontSize',15);
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Count','interpreter','latex','FontSize',15)
title('1C Cycle')
hold off;
saveas(gcf,'1C.png')
fig = figure(40);
left_color = [0 0 0];
right_color = [0 0.4470 0.7410];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(data.R4.echem.time, data.R4.echem.V,'LineWidth',2)
ylabel('Voltage [V]','interpreter','latex','FontSize',15)
yyaxis right
plot(data.R4.echem.time, data.R4.echem.I,'--','LineWidth',2)
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Current [A]','interpreter','latex','FontSize',15)
title('1C Cycle')
saveas(gcf,'1C_DC.png')




% 2C plots
figure(5);
hold on;
plot(data.R5.NMR.time, data.R5.NMR.total, 'DisplayName','Total Li', 'LineWidth',2)
plot(data.R5.NMR.time, data.R5.NMR.total_solid, 'DisplayName','Total Li - Solid Phase', 'LineWidth',2)
plot(data.R5.NMR.time, data.R5.NMR.dendrite, 'DisplayName','Dendrite Li', 'LineWidth',2)
plot(data.R5.NMR.time, data.R5.NMR.plated, 'DisplayName','Plated Li', 'LineWidth',2)
plot(data.R5.NMR.time, data.R5.NMR.phase3, 'DisplayName','Concentrated Li', 'LineWidth',2)
plot(data.R5.NMR.time, data.R5.NMR.phase2, 'DisplayName','Dilute Li', 'LineWidth',2)
plot(data.R5.NMR.time, data.R5.NMR.phase1, 'DisplayName','Disordered Li', 'LineWidth',2)
plot(data.R5.NMR.time, data.R5.NMR.electrolyte, 'DisplayName','Electrolyte Li', 'LineWidth',2)
legend('interpreter','latex','FontSize',15);
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Count','interpreter','latex','FontSize',15)
title('2C Cycle')
hold off;
saveas(gcf,'2C.png')
fig = figure(50);
left_color = [0 0 0];
right_color = [0 0.4470 0.7410];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(data.R5.echem.time, data.R5.echem.V,'LineWidth',2)
ylabel('Voltage [V]','interpreter','latex','FontSize',15)
yyaxis right
plot(data.R5.echem.time, data.R5.echem.I, '--','LineWidth',2)
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Current [A]','interpreter','latex','FontSize',15)
title('2C Cycle')
saveas(gcf,'2C_DC.png')



% 3C plots
figure(6);
hold on;
plot(data.R6.NMR.time, data.R6.NMR.total, 'DisplayName','Total Li', 'LineWidth',2)
plot(data.R6.NMR.time, data.R6.NMR.total_solid, 'DisplayName','Total Li - Solid Phase', 'LineWidth',2)
plot(data.R6.NMR.time, data.R6.NMR.dendrite, 'DisplayName','Dendrite Li', 'LineWidth',2)
plot(data.R6.NMR.time, data.R6.NMR.plated, 'DisplayName','Plated Li', 'LineWidth',2)
plot(data.R6.NMR.time, data.R6.NMR.phase3, 'DisplayName','Concentrated Li', 'LineWidth',2)
plot(data.R6.NMR.time, data.R6.NMR.phase2, 'DisplayName','Dilute Li', 'LineWidth',2)
plot(data.R6.NMR.time, data.R6.NMR.phase1, 'DisplayName','Disordered Li', 'LineWidth',2)
plot(data.R6.NMR.time, data.R6.NMR.electrolyte, 'DisplayName','Electrolyte Li', 'LineWidth',2)
legend('interpreter','latex','FontSize',15);
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Count','interpreter','latex','FontSize',15)
title('3C Cycle')
hold off;
saveas(gcf,'3C.png')
fig = figure(60);
left_color = [0 0 0];
right_color = [0 0.4470 0.7410];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plot(data.R6.echem.time, data.R6.echem.V,'LineWidth',2)
ylabel('Voltage [V]','interpreter','latex','FontSize',15)
yyaxis right
plot(data.R6.echem.time, data.R6.echem.I,'--','LineWidth',2)
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Current [A]','interpreter','latex','FontSize',15)
title('3C Cycle')
saveas(gcf,'3C_DC.png')












%% extract c1 and c2 data and plot
% C3 cycle
figure(11);
c1tilde = data.R2.NMR.phase1 + data.R2.NMR.phase2 + data.R2.NMR.phase3;
c2tilde = data.R2.NMR.plated;
cel = data.R2.NMR.electrolyte;
timeseq = data.R2.NMR.time*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R2.echem.time);
I = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),data.R2.NMR.time,'spline');
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$\widetilde{C}_e(t)$')
legend('interpreter','latex','FontSize',15)
hold off;
title('C3 Cycle')
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Concentrations','interpreter','latex','FontSize',15)
saveas(gcf,'C3_conc.png')



% C2 cycle
figure(12);        
data.R3.NMR.phase3(isnan(data.R3.NMR.phase3)) = 0;   % replace NAN with zero
c1tilde = data.R3.NMR.phase1 + data.R3.NMR.phase2  + data.R3.NMR.phase3;
c2tilde = data.R3.NMR.plated;
cel = data.R3.NMR.electrolyte;
timeseq = data.R3.NMR.time*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R3.echem.time);
I = interp1(data.R3.echem.time(indices), data.R3.echem.I(indices),data.R3.NMR.time,'spline');
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$\widetilde{C}_e(t)$')
legend('interpreter','latex','FontSize',15)
hold off;
title('C2 Cycle')
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Concentrations','interpreter','latex','FontSize',15)
saveas(gcf,'C2_conc.png')



% 1C cycle
figure(13);        
data.R4.NMR.phase3(isnan(data.R4.NMR.phase3)) = 0;   % replace NAN with zero
data.R4.NMR.dendrite(isnan(data.R4.NMR.dendrite)) = 0;   % replace NAN with zero
c1tilde = data.R4.NMR.phase1 + data.R4.NMR.phase2 + data.R4.NMR.phase3;
c2tilde = data.R4.NMR.plated + data.R4.NMR.dendrite;
cel = data.R4.NMR.electrolyte;
timeseq = data.R4.NMR.time*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R4.echem.time);
I = interp1(data.R4.echem.time(indices), data.R4.echem.I(indices),data.R4.NMR.time,'spline');
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$\widetilde{C}_e(t)$')
legend('interpreter','latex','FontSize',15)
hold off;
title('1C Cycle')
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Concentrations','interpreter','latex','FontSize',15)
saveas(gcf,'1C_conc.png')


% 2C cycle
figure(14);       
data.R5.NMR.phase3(isnan(data.R5.NMR.phase3)) = 0;   % replace NAN with zero
data.R5.NMR.dendrite(isnan(data.R5.NMR.dendrite)) = 0;   % replace NAN with zero
c1tilde = data.R5.NMR.phase1 + data.R5.NMR.phase2 + data.R5.NMR.phase3;
c2tilde = data.R5.NMR.plated + data.R5.NMR.dendrite;
cel = data.R5.NMR.electrolyte;
timeseq = data.R5.NMR.time*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R5.echem.time);
I = interp1(data.R5.echem.time(indices), data.R5.echem.I(indices),data.R5.NMR.time,'spline');
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$\widetilde{C}_e(t)$')
legend('interpreter','latex','FontSize',15)
hold off;
title('2C Cycle')
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Concentrations','interpreter','latex','FontSize',15)
saveas(gcf,'2C_conc.png')




% 3C cycle
figure(15);            
data.R6.NMR.phase3(isnan(data.R6.NMR.phase3)) = 0;   % replace NAN with zero
data.R6.NMR.dendrite(isnan(data.R6.NMR.dendrite)) = 0;   % replace NAN with zero
c1tilde = data.R6.NMR.phase1 + data.R6.NMR.phase2 + data.R6.NMR.phase3;
c2tilde = data.R6.NMR.plated + data.R6.NMR.dendrite;
cel = data.R6.NMR.electrolyte;
timeseq = data.R6.NMR.time*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R6.echem.time);
I = interp1(data.R6.echem.time(indices), data.R6.echem.I(indices),data.R6.NMR.time,'spline');
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$\widetilde{C}_e(t)$')
legend('interpreter','latex','FontSize',15)
hold off;
title('3C Cycle')
xlabel('Time [hr]','interpreter','latex','FontSize',15)
ylabel('Li Concentrations','interpreter','latex','FontSize',15)
saveas(gcf,'3C_conc.png')
