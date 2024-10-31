% Function to test the results of optimal reconstruction on other cycles and generate metrics for performance! 
close all;
clear; 
clc;
global c1grid c2grid c1alpha c1beta c2alpha c2beta Ngrid A L_a F a rho_1 rho_2 sf lambda Nt dC1 dC2 weight
global c1tilde c2tilde timeseq I C0


%%% Import a trained model, test on other cycles, plot forward
%%% concentrations vs true ones for that particular model:
% load('Results_R4.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
%     "omega_1_initial","omega_2_initial","omega_1", "omega_2", "omega", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
%     "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", "cm", "tm", "sf", ...
%     "A", "L_a", "r", "F", "a", "rho_1", "rho_2", "weight", "lambda_initial", "lambda", "lambda_history");
% 
% fields = ['R2'; 'R3'; 'R4'; 'R5'; 'R6'];
% costs = zeros(1,max(size(fields)));
% for i = 1:max(size(fields))
%     c1tilde = getfield(sdata, fields(i,:)).c1;
%     c2tilde = getfield(sdata, fields(i,:)).c2;
%     timeseq = getfield(sdata, fields(i,:)).t;
%     I = getfield(sdata, fields(i,:)).J;
%     C0 = [getfield(sdata, fields(i,:)).c1(1), getfield(sdata, fields(i,:)).c2(1)];
%     Nt = max(size(timeseq));
%     c = forward(omega_1, omega_2, lambda);
%     
%     % plot
%     figure()
%     plot(timeseq, c(:,1),'LineWidth',1.5,'LineStyle','-','Color',"#77AC30",'DisplayName', '$C_1(t;\widehat{\lambda},\widehat{\omega}_1,\widehat{\omega}_2)$')
%     hold on
%     plot(timeseq, c1tilde,'LineWidth',1.5,'LineStyle',':','Color',"#A2142F",'DisplayName', '$\widetilde{C}_1(t)$')
%     xlabel('$t$','interpreter','latex','FontSize',15)
%     ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
%     legend('interpreter','latex','Location','best');
%     hold off
%     %loc = append(pwd,foldername,'/_c1trajectory.png');
%     %saveas(gcf, loc);
%     
%     figure()
%     plot(timeseq, c(:,2),'LineWidth',1.5,'LineStyle','-','Color',"#77AC30",'DisplayName', '$C_2(t;\widehat{\lambda},\widehat{\omega}_1,\widehat{\omega}_2)$')
%     hold on
%     plot(timeseq, c2tilde,'LineWidth',1.5,'LineStyle',':','Color',"#A2142F",'DisplayName', '$\widetilde{C}_2(t)$')
%     xlabel('$t$','interpreter','latex','FontSize',15)
%     ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
%     legend('interpreter','latex','Location','best');
%     hold off
%     %loc = append(pwd,foldername,'/_c1trajectory.png');
%     %saveas(gcf, loc);
%     
%     costs(i) = J(c,c1tilde,c2tilde,timeseq);
% end






%%% training on one cycle  and testing on all others, computing costs, and
%%% plotting all costs.

fields = ['R2'; 'R3'; 'R4'; 'R5'; 'R6'];
cost = zeros(max(size(fields))+1,max(size(fields)));
for j=1:max(size(fields))
    data_to_load = strcat('Results_',fields(j,:),'.mat');
    load(data_to_load,"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "omega", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", "cm", "tm", "sf", ...
    "A", "L_a", "r", "F", "a", "rho_1", "rho_2", "weight", "lambda_initial", "lambda", "lambda_history");
    for i = 1:max(size(fields))
        c1tilde = getfield(sdata, fields(i,:)).c1;
        c2tilde = getfield(sdata, fields(i,:)).c2;
        timeseq = getfield(sdata, fields(i,:)).t;
        I = getfield(sdata, fields(i,:)).J;
        C0 = [getfield(sdata, fields(i,:)).c1(1), getfield(sdata, fields(i,:)).c2(1)];
        Nt = max(size(timeseq));
        c = forward(omega_1, omega_2, lambda);
        cost(j,i) = J(c,c1tilde,c2tilde,timeseq);
    end
end
% taking ensumble model and testing it on all cycles!
load('Results_ensumble.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "omega", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", "cm", "tm", "sf", ...
    "A", "L_a", "r", "F", "a", "rho_1", "rho_2", "weight", "lambda_initial", "lambda", "lambda_history");
for i = 1:max(size(fields))
    c1tilde = getfield(sdata, fields(i,:)).c1;
    c2tilde = getfield(sdata, fields(i,:)).c2;
    timeseq = getfield(sdata, fields(i,:)).t;
    I = getfield(sdata, fields(i,:)).J;
    C0 = [getfield(sdata, fields(i,:)).c1(1), getfield(sdata, fields(i,:)).c2(1)];
    Nt = max(size(timeseq));
    c = forward(omega_1, omega_2, lambda);
    cost(end,i) = J(c,c1tilde,c2tilde,timeseq);
end
    
    
    
    
cycles = {'C3 Cycle', 'C2 Cycle', '1C Cycle', '2C Cycle', '3C Cycle'};
% plot results of train test.
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
ylabel('$\mathcal{J}$','interpreter','latex','FontSize',15)
xlim([0.7, 5.3]);
legend();
hold off
% loc = append(pwd,foldername,'/_test.png');
% saveas(gcf, loc);








