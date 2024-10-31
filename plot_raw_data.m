function plot_raw_data(data)
%%% -----------------------------------------------------------------
% Plot raw data 
% Input: data struct
% Output: plots for each C-rate

% C20 plots
figure(1);
hold on;
plot(data.R1.NMR.time, data.R1.NMR.total, 'DisplayName','Total Li', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.total_solid, 'DisplayName','Total Li - Solid Phase', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.phase2, 'DisplayName','Dilute Li', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.phase1, 'DisplayName','Disordered Li', 'LineWidth',2)
plot(data.R1.NMR.time, data.R1.NMR.electrolyte, 'DisplayName','Electrolyte Li', 'LineWidth',2)
legend();
xlabel('Time [hr]')
ylabel('Li Count')
title('C20 Cycle')
hold off;
saveas(gcf,'C20.png')
figure(10);
plot(data.R1.echem.time, data.R1.echem.V,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Voltage [V]')
title('C20 Cycle')
saveas(gcf,'C20_voltage.png')
figure(100);
plot(data.R1.echem.time, data.R1.echem.I,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Current [A]')
title('C20 Cycle')
saveas(gcf,'C20_current.png')



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
legend();
xlabel('Time [hr]')
ylabel('Li Count')
title('C3 Cycle')
hold off;
saveas(gcf,'C3.png')
figure(20);
plot(data.R2.echem.time, data.R2.echem.V,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Voltage [V]')
title('C3 Cycle')
saveas(gcf,'C3_voltage.png')
figure(200);
plot(data.R2.echem.time, data.R2.echem.I,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Current [A]')
title('C3 Cycle')
saveas(gcf,'C3_current.png')


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
legend();
xlabel('Time [hr]')
ylabel('Li Count')
title('C2 Cycle')
hold off;
saveas(gcf,'C2.png')
figure(30);
plot(data.R3.echem.time, data.R3.echem.V,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Voltage [V]')
title('C2 Cycle')
saveas(gcf,'C2_voltage.png')
figure(300);
plot(data.R3.echem.time, data.R3.echem.I,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Current [A]')
title('C2 Cycle')
saveas(gcf,'C2_current.png')


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
legend();
xlabel('Time [hr]')
ylabel('Li Count')
title('1C Cycle')
hold off;
saveas(gcf,'1C.png')
figure(40);
plot(data.R4.echem.time, data.R4.echem.V,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Voltage [V]')
title('1C Cycle')
saveas(gcf,'1C_voltage.png')
figure(400);
plot(data.R4.echem.time, data.R4.echem.I,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Current [A]')
title('1C Cycle')
saveas(gcf,'1C_current.png')




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
legend();
xlabel('Time [hr]')
ylabel('Li Count')
title('2C Cycle')
hold off;
saveas(gcf,'2C.png')
figure(50);
plot(data.R5.echem.time, data.R5.echem.V,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Voltage [V]')
title('2C Cycle')
saveas(gcf,'2C_voltage.png')
figure(500);
plot(data.R5.echem.time, data.R5.echem.I,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Current [A]')
title('2C Cycle')
saveas(gcf,'2C_current.png')



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
legend();
xlabel('Time [hr]')
ylabel('Li Count')
title('3C Cycle')
hold off;
saveas(gcf,'3C.png')
figure(60);
plot(data.R6.echem.time, data.R6.echem.V,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Voltage [V]')
title('3C Cycle')
saveas(gcf,'3C_voltage.png')
figure(600);
plot(data.R6.echem.time, data.R6.echem.I,'LineWidth',2)
xlabel('Time [hr]')
ylabel('Current [A]')
title('3C Cycle')
saveas(gcf,'3C_current.png')












%% extract c1 and c2 data and plot
% C3 cycle
figure(11);
tt = data.R2.echem.time(data.R2.echem.steps == "Discharge C/3"); % time at which discharge is started
ss = size(data.R2.NMR.time);
ind = ss(1) - sum(data.R2.NMR.time>tt);              % index to start discharge
data.R2.NMR.phase3(isnan(data.R2.NMR.phase3)) = 0;   % replace NAN with zero
c1tilde = data.R2.NMR.phase1(1:ind) + data.R2.NMR.phase2(1:ind) + data.R2.NMR.phase3(1:ind);
c2tilde = data.R2.NMR.plated(1:ind);
cel = data.R2.NMR.electrolyte(1:ind);
timeseq = data.R2.NMR.time(1:ind)*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R2.echem.time);
I = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),data.R2.NMR.time(1:ind),'spline');
% I, c1(t), c2(t), and t are now ready.
% find c1alpha, c1beta, c2alpha, c2beta:
c1alpha = min(c1tilde);
c1beta  = max(c1tilde);
c2alpha = min(c2tilde);
c2beta  = max(c2tilde);
% discretize the space for state variables
Ngrid = 100;
c1grid = c1alpha:(c1beta - c1alpha)/Ngrid:c1beta;
c2grid = c2alpha:(c2beta - c2alpha)/Ngrid:c2beta;
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$C_e(t)$')
legend('Interpreter','latex')
hold off;
title('C3 Cycle')
xlabel('Time [hr]')
ylabel('Li Concentrations')
saveas(gcf,'C3_conc.png')



% C2 cycle
figure(12);
tt = data.R3.echem.time(data.R3.echem.steps == "C/3 Discharge");
ss = size(data.R3.NMR.time);
ind = ss(1) - sum(data.R3.NMR.time>tt);              
data.R3.NMR.phase3(isnan(data.R3.NMR.phase3)) = 0;   % replace NAN with zero
c1tilde = data.R3.NMR.phase1(1:ind) + data.R3.NMR.phase2(1:ind) + data.R3.NMR.phase3(1:ind);
c2tilde = data.R3.NMR.plated(1:ind);
cel = data.R3.NMR.electrolyte(1:ind);
timeseq = data.R3.NMR.time(1:ind)*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R3.echem.time);
I = interp1(data.R3.echem.time(indices), data.R3.echem.I(indices),data.R3.NMR.time(1:ind),'spline');
% I, c1(t), c2(t), and t are now ready.
% find c1alpha, c1beta, c2alpha, c2beta:
c1alpha = min(c1tilde);
c1beta  = max(c1tilde);
c2alpha = min(c2tilde);
c2beta  = max(c2tilde);
% discretize the space for state variables
Ngrid = 100;
c1grid = c1alpha:(c1beta - c1alpha)/Ngrid:c1beta;
c2grid = c2alpha:(c2beta - c2alpha)/Ngrid:c2beta;
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$C_e(t)$')
legend('Interpreter','latex')
hold off;
title('C2 Cycle')
xlabel('Time [hr]')
ylabel('Li Concentrations')
saveas(gcf,'C2_conc.png')



% 1C cycle
figure(13);
tt = data.R4.echem.time(data.R4.echem.steps == "C/3 Discharge");
ss = size(data.R4.NMR.time);
ind = ss(1) - sum(data.R4.NMR.time>tt);              
data.R4.NMR.phase3(isnan(data.R4.NMR.phase3)) = 0;   % replace NAN with zero
data.R4.NMR.dendrite(isnan(data.R4.NMR.dendrite)) = 0;   % replace NAN with zero
c1tilde = data.R4.NMR.phase1(1:ind) + data.R4.NMR.phase2(1:ind) + data.R4.NMR.phase3(1:ind);
c2tilde = data.R4.NMR.plated(1:ind) + data.R4.NMR.dendrite(1:ind);
cel = data.R4.NMR.electrolyte(1:ind);
timeseq = data.R4.NMR.time(1:ind)*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R4.echem.time);
I = interp1(data.R4.echem.time(indices), data.R4.echem.I(indices),data.R4.NMR.time(1:ind),'spline');
% I, c1(t), c2(t), and t are now ready.
% find c1alpha, c1beta, c2alpha, c2beta:
c1alpha = min(c1tilde);
c1beta  = max(c1tilde);
c2alpha = min(c2tilde);
c2beta  = max(c2tilde);
% discretize the space for state variables
Ngrid = 100;
c1grid = c1alpha:(c1beta - c1alpha)/Ngrid:c1beta;
c2grid = c2alpha:(c2beta - c2alpha)/Ngrid:c2beta;
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$C_e(t)$')
legend('Interpreter','latex')
hold off;
title('1C Cycle')
xlabel('Time [hr]')
ylabel('Li Concentrations')
saveas(gcf,'1C_conc.png')


% 2C cycle
figure(14);
tt = data.R5.echem.time(data.R5.echem.steps == "C/3 discahrge");
ss = size(data.R5.NMR.time);
ind = ss(1) - sum(data.R5.NMR.time>tt);              
data.R5.NMR.phase3(isnan(data.R5.NMR.phase3)) = 0;   % replace NAN with zero
data.R5.NMR.dendrite(isnan(data.R5.NMR.dendrite)) = 0;   % replace NAN with zero
c1tilde = data.R5.NMR.phase1(1:ind) + data.R5.NMR.phase2(1:ind) + data.R5.NMR.phase3(1:ind);
c2tilde = data.R5.NMR.plated(1:ind) + data.R5.NMR.dendrite(1:ind);
cel = data.R5.NMR.electrolyte(1:ind);
timeseq = data.R5.NMR.time(1:ind)*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R5.echem.time);
I = interp1(data.R5.echem.time(indices), data.R5.echem.I(indices),data.R5.NMR.time(1:ind),'spline');
% I, c1(t), c2(t), and t are now ready.
% find c1alpha, c1beta, c2alpha, c2beta:
c1alpha = min(c1tilde);
c1beta  = max(c1tilde);
c2alpha = min(c2tilde);
c2beta  = max(c2tilde);
% discretize the space for state variables
Ngrid = 100;
c1grid = c1alpha:(c1beta - c1alpha)/Ngrid:c1beta;
c2grid = c2alpha:(c2beta - c2alpha)/Ngrid:c2beta;
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$C_e(t)$')
legend('Interpreter','latex')
hold off;
title('2C Cycle')
xlabel('Time [hr]')
ylabel('Li Concentrations')
saveas(gcf,'2C_conc.png')




% 3C cycle
figure(15);
tt = data.R6.echem.time(data.R6.echem.steps == "C/3 Discharge");
ss = size(data.R6.NMR.time);
ind = ss(1) - sum(data.R6.NMR.time>tt);              
data.R6.NMR.phase3(isnan(data.R6.NMR.phase3)) = 0;   % replace NAN with zero
data.R6.NMR.dendrite(isnan(data.R6.NMR.dendrite)) = 0;   % replace NAN with zero
c1tilde = data.R6.NMR.phase1(1:ind) + data.R6.NMR.phase2(1:ind) + data.R6.NMR.phase3(1:ind);
c2tilde = data.R6.NMR.plated(1:ind) + data.R6.NMR.dendrite(1:ind);
cel = data.R6.NMR.electrolyte(1:ind);
timeseq = data.R6.NMR.time(1:ind)*3600;
% normalize concentrations
c1tilde = c1tilde./max(c1tilde);
c2tilde = c2tilde./max(c2tilde);
cel = cel./max(cel);
% interpolate current I to NMR data points
[~,indices]=unique(data.R6.echem.time);
I = interp1(data.R6.echem.time(indices), data.R6.echem.I(indices),data.R6.NMR.time(1:ind),'spline');
% I, c1(t), c2(t), and t are now ready.
% find c1alpha, c1beta, c2alpha, c2beta:
c1alpha = min(c1tilde);
c1beta  = max(c1tilde);
c2alpha = min(c2tilde);
c2beta  = max(c2tilde);
% discretize the space for state variables
Ngrid = 100;
c1grid = c1alpha:(c1beta - c1alpha)/Ngrid:c1beta;
c2grid = c2alpha:(c2beta - c2alpha)/Ngrid:c2beta;
hold on;
plot(timeseq./3600,c1tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_1(t)$')
plot(timeseq./3600,c2tilde,'LineWidth',2,'DisplayName','$\widetilde{C}_2(t)$')
plot(timeseq./3600,cel,'LineWidth',2,'DisplayName','$C_e(t)$')
legend('Interpreter','latex')
hold off;
title('3C Cycle')
xlabel('Time [hr]')
ylabel('Li Concentrations')
saveas(gcf,'3C_conc.png')
