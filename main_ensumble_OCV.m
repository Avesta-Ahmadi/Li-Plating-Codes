% Main function to perform analysis on OCV part of data, to get beta
% parameters (the parameters of the linear model).
% data is split into three phases: charge, OCV, discharge
% the OCV part of data is used for training.
% all cycles are used for training

clear; clc; close all;

global c1grid c2grid c1alpha c1beta c2alpha c2beta Ngrid sf Nt dC1 dC2 weight
global c1tilde c2tilde timeseq I C0 weight

weight = 50;
% load data
data = load('data.mat');
load('sdata.mat');




%% plot data
plt = false;
if plt
    figure()
    plot(sdata.R1.t, sdata.R1.c1,'LineWidth',2)
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$\widetilde{C}_1(t)$','interpreter','latex','FontSize',15)
    title('C20 Cycle');
    
    figure()
    plot(sdata.R2.t, sdata.R2.c1,'LineWidth',2,'DisplayName', '$\widetilde{C}_1(t)$')
    hold on 
    plot(sdata.R2.t, sdata.R2.c2,'LineWidth',2,'DisplayName', '$\widetilde{C}_2(t)$')
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$\widetilde{C}(t)$','interpreter','latex','FontSize',15)
    legend('interpreter','latex','FontSize',15)
    title('C3 Cycle');
    
    figure()
    plot(sdata.R3.t, sdata.R3.c1,'LineWidth',2,'DisplayName', '$\widetilde{C}_1(t)$')
    hold on 
    plot(sdata.R3.t, sdata.R3.c2,'LineWidth',2,'DisplayName', '$\widetilde{C}_2(t)$')
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$\widetilde{C}(t)$','interpreter','latex','FontSize',15)
    legend('interpreter','latex','FontSize',15)
    title('C2 Cycle');
    
    figure()
    plot(sdata.R4.t, sdata.R4.c1,'LineWidth',2,'DisplayName', '$\widetilde{C}_1(t)$')
    hold on 
    plot(sdata.R4.t, sdata.R4.c2,'LineWidth',2,'DisplayName', '$\widetilde{C}_2(t)$')
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$\widetilde{C}(t)$','interpreter','latex','FontSize',15)
    legend('interpreter','latex','FontSize',15)
    title('1C Cycle');
    
    figure()
    plot(sdata.R5.t, sdata.R5.c1,'LineWidth',2,'DisplayName', '$\widetilde{C}_1(t)$')
    hold on 
    plot(sdata.R5.t, sdata.R5.c2,'LineWidth',2,'DisplayName', '$\widetilde{C}_2(t)$')
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$\widetilde{C}(t)$','interpreter','latex','FontSize',15)
    legend('interpreter','latex','FontSize',15)
    title('2C Cycle');
    
    figure()
    plot(sdata.R6.t, sdata.R6.c1,'LineWidth',2,'DisplayName', '$\widetilde{C}_1(t)$')
    hold on 
    plot(sdata.R6.t, sdata.R6.c2,'LineWidth',2,'DisplayName', '$\widetilde{C}_2(t)$')
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$\widetilde{C}(t)$','interpreter','latex','FontSize',15)
    legend('interpreter','latex','FontSize',15)
    title('3C Cycle');
end





%% running iterative algorithm for finding optimal parameters of the linear model for OCV part
c1alpha = -0.5;
c1beta  = 1.5;
c2alpha = -0.2;
c2beta  = 0.5;
Ngrid   = 5000;
dC1 = (c1beta-c1alpha)/(Ngrid-1);
dC2 = (c2beta-c2alpha)/(Ngrid-1);
c1grid = c1alpha:dC1:c1beta;
c2grid = c2alpha:dC2:c2beta;

c1tilde = struct;
c2tilde = struct;
timeseq = struct;
I = struct;
C0 = struct;
Nt = struct;

c1tilde.R2 = sdata.R2.c1(sdata.R2.label == 'OCV');
c1tilde.R3 = sdata.R3.c1(sdata.R3.label == 'OCV');
c1tilde.R4 = sdata.R4.c1(sdata.R4.label == 'OCV');
c1tilde.R5 = sdata.R5.c1(sdata.R5.label == 'OCV');
c1tilde.R6 = sdata.R6.c1(sdata.R6.label == 'OCV');
c2tilde.R2 = sdata.R2.c2(sdata.R2.label == 'OCV');
c2tilde.R3 = sdata.R3.c2(sdata.R3.label == 'OCV');
c2tilde.R4 = sdata.R4.c2(sdata.R4.label == 'OCV');
c2tilde.R5 = sdata.R5.c2(sdata.R5.label == 'OCV');
c2tilde.R6 = sdata.R6.c2(sdata.R6.label == 'OCV');
timeseq.R2 = sdata.R2.t(sdata.R2.label == 'OCV');
timeseq.R3 = sdata.R3.t(sdata.R3.label == 'OCV');
timeseq.R4 = sdata.R4.t(sdata.R4.label == 'OCV');
timeseq.R5 = sdata.R5.t(sdata.R5.label == 'OCV');
timeseq.R6 = sdata.R6.t(sdata.R6.label == 'OCV');
I.R2 = sdata.R2.J(sdata.R2.label == 'OCV');
I.R3 = sdata.R3.J(sdata.R3.label == 'OCV');
I.R4 = sdata.R4.J(sdata.R4.label == 'OCV');
I.R5 = sdata.R5.J(sdata.R5.label == 'OCV');
I.R6 = sdata.R6.J(sdata.R6.label == 'OCV');
C0.R2 = [c1tilde.R2(1), c2tilde.R2(1)];
C0.R3 = [c1tilde.R3(1), c2tilde.R3(1)];
C0.R4 = [c1tilde.R4(1), c2tilde.R4(1)];
C0.R5 = [c1tilde.R5(1), c2tilde.R5(1)];
C0.R6 = [c1tilde.R6(1), c2tilde.R6(1)];
Nt.R2 = max(size(timeseq.R2));
Nt.R3 = max(size(timeseq.R3));
Nt.R4 = max(size(timeseq.R4));
Nt.R5 = max(size(timeseq.R5));
Nt.R6 = max(size(timeseq.R6));


% Run iterative procedure
Feval = 0;
plt = false;
mkdir OCV;
beta_initial = [-0.1, 0.1, 0.1, -0.1];
foldername = '/OCV';
[beta, c1forward_history, c2forward_history, pjpb_history, beta_history, Jtau_history] = iterative_ensumble_OCV(beta_initial, plt, foldername);

% save variables
save('OCV/Results_OCV.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt", "I", "timeseq",...
    "beta_initial","beta", "c1forward_history", "c2forward_history", "pjpb_history", "beta_history", "sdata", "cm", "cm1", "cm2", "tm", "sf", "weight", "Jtau_history");


beta
A = [beta(2), 0.006*beta(3); beta(4), -beta(3)];
eig(A)

%% Plot results
% plot trajectory of concentrations
figure('visible','on')
plot(timeseq.R2, c1forward_history.R2(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\beta^{(0)})$')
hold on
plot(timeseq.R2, c1forward_history.R2(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\beta})$')
plot(timeseq.R2, c1tilde.R2,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of C3 cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c1trajectory_R2.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R3, c1forward_history.R3(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\beta^{(0)})$')
hold on
plot(timeseq.R3, c1forward_history.R3(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\beta})$')
plot(timeseq.R3, c1tilde.R3,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of C2 cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c1trajectory_R3.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R4, c1forward_history.R4(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\beta^{(0)})$')
hold on
plot(timeseq.R4, c1forward_history.R4(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\beta})$')
plot(timeseq.R4, c1tilde.R4,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of 1C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c1trajectory_R4.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R5, c1forward_history.R5(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\beta^{(0)})$')
hold on
plot(timeseq.R5, c1forward_history.R5(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\beta})$')
plot(timeseq.R5, c1tilde.R5,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of 2C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c1trajectory_R5.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R6, c1forward_history.R6(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\beta^{(0)})$')
hold on
plot(timeseq.R6, c1forward_history.R6(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\beta})$')
plot(timeseq.R6, c1tilde.R6,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of 3C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c1trajectory_R6.png');
saveas(gcf, loc);


figure('visible','on')
plot(timeseq.R2, c2forward_history.R2(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\beta^{(0)})$')
hold on
plot(timeseq.R2, c2forward_history.R2(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\beta})$')
plot(timeseq.R2, c2tilde.R2,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of C3 cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c2trajectory_R2.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R3, c2forward_history.R3(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\beta^{(0)})$')
hold on
plot(timeseq.R3, c2forward_history.R3(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\beta})$')
plot(timeseq.R3, c2tilde.R3,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of C2 cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c2trajectory_R3.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R4, c2forward_history.R4(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\beta^{(0)})$')
hold on
plot(timeseq.R4, c2forward_history.R4(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\beta})$')
plot(timeseq.R4, c2tilde.R4,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of 1C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c2trajectory_R4.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R5, c2forward_history.R5(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\beta^{(0)})$')
hold on
plot(timeseq.R5, c2forward_history.R5(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\beta})$')
plot(timeseq.R5, c2tilde.R5,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of 2C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c2trajectory_R5.png');
saveas(gcf, loc);

figure('visible','on')
plot(timeseq.R6, c2forward_history.R6(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\beta^{(0)})$')
hold on
plot(timeseq.R6, c2forward_history.R6(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\beta})$')
plot(timeseq.R6, c2tilde.R6,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of 3C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c2trajectory_R6.png');
saveas(gcf, loc);

