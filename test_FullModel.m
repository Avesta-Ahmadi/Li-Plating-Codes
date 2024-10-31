% test the results of trained models. Different models are imported, tested
% on all cycles, and plotted. This includes two types of models: models
% trained on charge states, and models trained on discharge states.


clear; clc; close all;

global c1grid c2grid c1alpha c1beta c2alpha c2beta Ngrid Nt weight beta
global c1tilde c2tilde timeseq I C0

weiht=50;
% load data
data = load('data.mat');
load('sdata.mat');

% load optimal parameters of OCV part
load('Results_OCV.mat', "beta");


fields = ['R2'; 'R3'; 'R4'; 'R5'; 'R6'];
cycles = {'C3 Cycle', 'C2 Cycle', '1C Cycle', '2C Cycle', '3C Cycle'};



%% test and plot - for models trained on CHARGE phase
cost = zeros(max(size(fields))+1,max(size(fields)));

for j=1:max(size(fields))
    data_to_load = strcat('Results_Charge_',fields(j,:),'.mat');
    load(data_to_load,"omega_1", "omega_2","lambda");
    for i = 1:max(size(fields))
        c1tilde = getfield(sdata, fields(i,:)).c1(getfield(sdata, fields(i,:)).label == 'charge');
        c2tilde = getfield(sdata, fields(i,:)).c2(getfield(sdata, fields(i,:)).label == 'charge');
        timeseq = getfield(sdata, fields(i,:)).t(getfield(sdata, fields(i,:)).label == 'charge');
        I = getfield(sdata, fields(i,:)).J(getfield(sdata, fields(i,:)).label == 'charge');
        C0 = [c1tilde(1), c2tilde(1)];
        Nt = max(size(timeseq));
        c = forward_FullModel(omega_1, omega_2, lambda);
        cost(j,i) = J(c,c1tilde,c2tilde,timeseq);
    end
end

% taking ensumble model and testing it on all cycles:
load('Results_Charge_ensumble.mat',"omega_1", "omega_2","lambda");
for i = 1:max(size(fields))
    c1tilde = getfield(sdata, fields(i,:)).c1(getfield(sdata, fields(i,:)).label == 'charge');
    c2tilde = getfield(sdata, fields(i,:)).c2(getfield(sdata, fields(i,:)).label == 'charge');
    timeseq = getfield(sdata, fields(i,:)).t(getfield(sdata, fields(i,:)).label == 'charge');
    I = getfield(sdata, fields(i,:)).J(getfield(sdata, fields(i,:)).label == 'charge');
    C0 = [c1tilde(1), c2tilde(1)];
    Nt = max(size(timeseq));
    c = forward_FullModel(omega_1, omega_2, lambda);
    cost(end,i) = J(c,c1tilde,c2tilde,timeseq);
end


% plot results
figure('visible','on')
plot(cost(1,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on C3 Cycle')
set(gca,'XTick',1:max(size(fields)),'XTickLabel',cycles)
set(gca, 'YScale', 'log')
hold on
plot(cost(2,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on C2 Cycle')
plot(cost(3,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on 1C Cycle')
plot(cost(4,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on 2C Cycle')
plot(cost(5,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on 3C Cycle')
plot(cost(6,:),'LineWidth',2, 'LineStyle','--', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on all Cycles')
ylabel('$\mathcal{J}_2$','interpreter','latex','FontSize',15)
xlim([0.7, 5.3]);
legend();
hold off
% loc = append(pwd,foldername,'/_test.png');
% saveas(gcf, loc);







%% test and plot - for models trained on DISCHARGE phase
cost = zeros(max(size(fields))+1,max(size(fields)));

for j=1:max(size(fields))
    data_to_load = strcat('Results_Discharge_',fields(j,:),'.mat');
    load(data_to_load,"omega_1","omega_2","lambda");
    for i = 1:max(size(fields))
        c1tilde = getfield(sdata, fields(i,:)).c1(getfield(sdata, fields(i,:)).label == 'discharge');
        c2tilde = getfield(sdata, fields(i,:)).c2(getfield(sdata, fields(i,:)).label == 'discharge');
        timeseq = getfield(sdata, fields(i,:)).t(getfield(sdata, fields(i,:)).label == 'discharge');
        I = getfield(sdata, fields(i,:)).J(getfield(sdata, fields(i,:)).label == 'discharge');
        C0 = [c1tilde(1), c2tilde(1)];
        Nt = max(size(timeseq));
        c = forward_FullModel(omega_1, omega_2, lambda);
        cost(j,i) = J(c,c1tilde,c2tilde,timeseq);
    end
end

% taking ensumble model and testing it on all cycles:
load('Results_Discharge_ensumble.mat',"omega_1", "omega_2","lambda");
for i = 1:max(size(fields))
    c1tilde = getfield(sdata, fields(i,:)).c1(getfield(sdata, fields(i,:)).label == 'discharge');
    c2tilde = getfield(sdata, fields(i,:)).c2(getfield(sdata, fields(i,:)).label == 'discharge');
    timeseq = getfield(sdata, fields(i,:)).t(getfield(sdata, fields(i,:)).label == 'discharge');
    I = getfield(sdata, fields(i,:)).J(getfield(sdata, fields(i,:)).label == 'discharge');
    C0 = [c1tilde(1), c2tilde(1)];
    Nt = max(size(timeseq));
    c = forward_FullModel(omega_1, omega_2, lambda);
    cost(end,i) = J(c,c1tilde,c2tilde,timeseq);
end


% plot results
figure('visible','on')
plot(cost(1,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on C3 Cycle')
set(gca,'XTick',1:max(size(fields)),'XTickLabel',cycles)
set(gca, 'YScale', 'log')
hold on
plot(cost(2,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on C2 Cycle')
plot(cost(3,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on 1C Cycle')
plot(cost(4,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on 2C Cycle')
plot(cost(5,:),'LineWidth',2, 'LineStyle','-', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on 3C Cycle')
plot(cost(6,:),'LineWidth',2, 'LineStyle','--', 'Marker', '.', 'MarkerSize', 20, 'DisplayName', 'Trained on all Cycles')
ylabel('$\mathcal{J}_2$','interpreter','latex','FontSize',15)
xlim([0.7, 5.3]);
legend();
hold off
% loc = append(pwd,foldername,'/_test.png');
% saveas(gcf, loc);











%% take ensumble model for charge and discharge, and OCV model for OCV part, compute trajectory of concentrations using these three models.

cycles = ['C3 Cycle'; 'C2 Cycle'; '1C Cycle'; '2C Cycle'; '3C Cycle'];
for i = 1:max(size(fields))
    
    % charge
    load('Results_Charge_ensumble.mat',"omega_1", "omega_2","lambda");
    c1tilde = getfield(sdata, fields(i,:)).c1(getfield(sdata, fields(i,:)).label == 'charge');
    c2tilde = getfield(sdata, fields(i,:)).c2(getfield(sdata, fields(i,:)).label == 'charge');
    timeseq = getfield(sdata, fields(i,:)).t(getfield(sdata, fields(i,:)).label == 'charge');
    I = getfield(sdata, fields(i,:)).J(getfield(sdata, fields(i,:)).label == 'charge');
    C0 = [c1tilde(1), c2tilde(1)];
    Nt = max(size(timeseq));
    c_charge = forward_FullModel(omega_1, omega_2, lambda);
    
    % OCV
    c1tilde = getfield(sdata, fields(i,:)).c1(getfield(sdata, fields(i,:)).label == 'OCV');
    c2tilde = getfield(sdata, fields(i,:)).c2(getfield(sdata, fields(i,:)).label == 'OCV');
    timeseq = getfield(sdata, fields(i,:)).t(getfield(sdata, fields(i,:)).label == 'OCV');
    I = getfield(sdata, fields(i,:)).J(getfield(sdata, fields(i,:)).label == 'OCV');
    C0 = c_charge(end,:);
    Nt = max(size(timeseq));
    c_OCV = forward_ensumble_OCV(beta, timeseq, C0);
        
    % discharge
    load('Results_Discharge_ensumble.mat',"omega_1", "omega_2","lambda");
    c1tilde = getfield(sdata, fields(i,:)).c1(getfield(sdata, fields(i,:)).label == 'discharge');
    c2tilde = getfield(sdata, fields(i,:)).c2(getfield(sdata, fields(i,:)).label == 'discharge');
    timeseq = getfield(sdata, fields(i,:)).t(getfield(sdata, fields(i,:)).label == 'discharge');
    I = getfield(sdata, fields(i,:)).J(getfield(sdata, fields(i,:)).label == 'discharge');
    C0 = c_OCV(end,:);
    Nt = max(size(timeseq));
    c_discharge = forward_FullModel(omega_1, omega_2, lambda);
    
    c = [c_charge; c_OCV; c_discharge];
    timeseq = getfield(sdata, fields(i,:)).t;
    c1tilde = getfield(sdata, fields(i,:)).c1;
    c2tilde = getfield(sdata, fields(i,:)).c2;
    
    % plot
    figure()
    plot(timeseq, c(:,1),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\beta},\bar{\alpha},\bar{\omega})$')
    hold on
    plot(timeseq, c1tilde,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', strcat('$\widetilde{C}_1(t)$'," for ",cycles(i,:)))
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
    legend('interpreter','latex','Location','best','FontSize',15);
    hold off
    
    
    figure()
    plot(timeseq, c(:,2),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\beta},\bar{\alpha},\bar{\omega})$')
    hold on
    plot(timeseq, c2tilde,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', strcat('$\widetilde{C}_2(t)$'," for ",cycles(i,:)))
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
    legend('interpreter','latex','Location','best','FontSize',15);
    hold off

end





%% ploting omega times lambda
% ploting initial guess for omega, and reconstructed omega!
load('Results_Charge_ensumble.mat',"omega_1", "omega_2","lambda");
Nn=50;
c1gridn = c1alpha:(c1beta-c1alpha)/Nn:c1beta;
c2gridn = c2alpha:(c2beta-c2alpha)/Nn:c2beta;
omega_1n = interp1(c1grid, omega_1, c1gridn,'spline');
omega_2n = interp1(c2grid, omega_2, c2gridn,'spline');
[X,Y] = meshgrid(omega_1n,omega_2n);
omegan = X.*Y;
omegan = omegan.*lambda;
figure('visible','on')
surf(c1gridn, c2gridn, omegan, 'DisplayName', 'Reconstructed, $\bar{\omega}(C_1,C_2)$', 'EdgeColor', "#0072BD", 'FaceColor', "#0072BD", 'FaceAlpha',0.5)
legend('interpreter','latex','Location','best','FontSize',15)
xlabel('$C_1$','interpreter','latex','FontSize',15)
ylabel('$C_2$','interpreter','latex','FontSize',15)
zlabel('$\omega(C_1,C_2)$','interpreter','latex','FontSize',15)



load('Results_Discharge_ensumble.mat',"omega_1", "omega_2","lambda");
Nn=50;
c1gridn = c1alpha:(c1beta-c1alpha)/Nn:c1beta;
c2gridn = c2alpha:(c2beta-c2alpha)/Nn:c2beta;
omega_1n = interp1(c1grid, omega_1, c1gridn,'spline');
omega_2n = interp1(c2grid, omega_2, c2gridn,'spline');
[X,Y] = meshgrid(omega_1n,omega_2n);
omegan = X.*Y;
omegan = omegan.*lambda;
figure('visible','on')
surf(c1gridn, c2gridn, omegan, 'DisplayName', 'Reconstructed, $\bar{\omega}(C_1,C_2)$', 'EdgeColor', "#0072BD", 'FaceColor', "#0072BD", 'FaceAlpha',0.5)
legend('interpreter','latex','Location','best','FontSize',15)
xlabel('$C_1$','interpreter','latex','FontSize',15)
ylabel('$C_2$','interpreter','latex','FontSize',15)
zlabel('$\omega(C_1,C_2)$','interpreter','latex','FontSize',15)

