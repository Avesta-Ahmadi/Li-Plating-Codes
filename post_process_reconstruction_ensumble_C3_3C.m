function [] = post_process_reconstruction_ensumble_C3_3C(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, lambda_history, foldername)
%%% -----------------------------------------------------------------
% Post process simulation results

global c1grid c2grid timeseq c1alpha c2alpha c1beta c2beta c1tilde c2tilde

% plot omega_1 and omega_2
figure('visible','on')
plot(c1grid, omega_1_initial,'LineWidth',2, 'LineStyle','--','Color',"#A2142F",'DisplayName', 'Initial Guess, $\omega_1^{(0)}(C_1)$')
hold on
plot(c1grid, omega_1,'LineWidth',2,'Color',"#0072BD",'DisplayName', 'Reconstructed, $\bar{\omega}_1(C_1)$')
xlabel('$C_1$','interpreter','latex','FontSize',15)
ylabel('$\omega_1(C_1)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
min1 = min([min(c1forward_history.R2(end,:)), min(c1forward_history.R6(end,:))]);
max1 = max([max(c1forward_history.R2(end,:)), max(c1forward_history.R6(end,:))]);
xline(min1,'--','HandleVisibility','off');
xline(max1,'--','HandleVisibility','off');
%xlim([0,1.2]);
hold off
loc = append(pwd,foldername,'/_omega_1_optimal.png');
saveas(gcf, loc);


figure('visible','on')
plot(c2grid, omega_2_initial,'LineWidth',2, 'LineStyle','--','Color',"#A2142F",'DisplayName', 'Initial Guess, $\omega_2^{(0)}(C_2)$')
hold on
plot(c2grid, omega_2,'LineWidth', 2,'Color',"#0072BD",'DisplayName', 'Reconstructed, $\bar{\omega}_2(C_2)$')
xlabel('$C_2$','interpreter','latex','FontSize',15)
ylabel('$\omega_2(C_2)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
min2 = min([min(c2forward_history.R2(end,:)), min(c2forward_history.R6(end,:))]);
max2 = max([max(c2forward_history.R2(end,:)), max(c2forward_history.R6(end,:))]);
xline(min2,'--','HandleVisibility','off');
xline(max2,'--','HandleVisibility','off');
%xlim([0,0.2]);
hold off
loc = append(pwd,foldername,'/_omega_2_optimal.png');
saveas(gcf, loc);






% plot trajectory of concentrations
figure('visible','on')
plot(timeseq.R2, c1forward_history.R2(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
hold on
plot(timeseq.R2, c1forward_history.R2(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
plot(timeseq.R2, c1tilde.R2,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of C3 cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c1trajectory_R2.png');
saveas(gcf, loc);
% figure('visible','on')
% plot(timeseq.R3, c1forward_history.R3(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
% hold on
% plot(timeseq.R3, c1forward_history.R3(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
% plot(timeseq.R3, c1tilde.R3,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of C2 cycle')
% xlabel('$t$','interpreter','latex','FontSize',15)
% ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
% legend('interpreter','latex','Location','best','FontSize',15);
% hold off
% loc = append(pwd,foldername,'/_c1trajectory_R3.png');
% saveas(gcf, loc);
% figure('visible','on')
% plot(timeseq.R4, c1forward_history.R4(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
% hold on
% plot(timeseq.R4, c1forward_history.R4(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
% plot(timeseq.R4, c1tilde.R4,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of 1C cycle')
% xlabel('$t$','interpreter','latex','FontSize',15)
% ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
% legend('interpreter','latex','Location','best','FontSize',15);
% hold off
% loc = append(pwd,foldername,'/_c1trajectory_R4.png');
% saveas(gcf, loc);
% figure('visible','on')
% plot(timeseq.R5, c1forward_history.R5(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
% hold on
% plot(timeseq.R5, c1forward_history.R5(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
% plot(timeseq.R5, c1tilde.R5,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of 2C cycle')
% xlabel('$t$','interpreter','latex','FontSize',15)
% ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
% legend('interpreter','latex','Location','best','FontSize',15);
% hold off
% loc = append(pwd,foldername,'/_c1trajectory_R5.png');
% saveas(gcf, loc);
figure('visible','on')
plot(timeseq.R6, c1forward_history.R6(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_1(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
hold on
plot(timeseq.R6, c1forward_history.R6(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_1(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
plot(timeseq.R6, c1tilde.R6,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_1(t)$ of 3C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_1(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c1trajectory_R6.png');
saveas(gcf, loc);


figure('visible','on')
plot(timeseq.R2, c2forward_history.R2(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
hold on
plot(timeseq.R2, c2forward_history.R2(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
plot(timeseq.R2, c2tilde.R2,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of C3 cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c2trajectory_R2.png');
saveas(gcf, loc);
% figure('visible','on')
% plot(timeseq.R3, c2forward_history.R3(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
% hold on
% plot(timeseq.R3, c2forward_history.R3(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
% plot(timeseq.R3, c2tilde.R3,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of C2 cycle')
% xlabel('$t$','interpreter','latex','FontSize',15)
% ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
% legend('interpreter','latex','Location','best','FontSize',15);
% hold off
% loc = append(pwd,foldername,'/_c2trajectory_R3.png');
% saveas(gcf, loc);
% figure('visible','on')
% plot(timeseq.R4, c2forward_history.R4(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
% hold on
% plot(timeseq.R4, c2forward_history.R4(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
% plot(timeseq.R4, c2tilde.R4,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of 1C cycle')
% xlabel('$t$','interpreter','latex','FontSize',15)
% ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
% legend('interpreter','latex','Location','best','FontSize',15);
% hold off
% loc = append(pwd,foldername,'/_c2trajectory_R4.png');
% saveas(gcf, loc);
% figure('visible','on')
% plot(timeseq.R5, c2forward_history.R5(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
% hold on
% plot(timeseq.R5, c2forward_history.R5(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
% plot(timeseq.R5, c2tilde.R5,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of 2C cycle')
% xlabel('$t$','interpreter','latex','FontSize',15)
% ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
% legend('interpreter','latex','Location','best','FontSize',15);
% hold off
% loc = append(pwd,foldername,'/_c2trajectory_R5.png');
% saveas(gcf, loc);
figure('visible','on')
plot(timeseq.R6, c2forward_history.R6(1,:),'LineWidth',1.5,'LineStyle','--','Color',"#A2142F",'DisplayName', '$C_2(t;\alpha^{(0)},\omega_1^{(0)},\omega_2^{(0)})$')
hold on
plot(timeseq.R6, c2forward_history.R6(end,:),'LineWidth',1.5,'LineStyle','-','Color',"#0072BD",'DisplayName', '$C_2(t;\bar{\alpha},\bar{\omega}_1,\bar{\omega}_2)$')
plot(timeseq.R6, c2tilde.R6,'LineWidth',1.5,'LineStyle',':','Color',"#77AC30",'DisplayName', '$\widetilde{C}_2(t)$ of 3C cycle')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$C_2(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','Location','best','FontSize',15);
hold off
loc = append(pwd,foldername,'/_c2trajectory_R6.png');
saveas(gcf, loc);








% ploting initial guess for omega, and reconstructed omega!
Nn=50;
c1gridn = c1alpha:(c1beta-c1alpha)/Nn:c1beta;
c2gridn = c2alpha:(c2beta-c2alpha)/Nn:c2beta;
omega_1n = interp1(c1grid, omega_1, c1gridn,'spline');
omega_2n = interp1(c2grid, omega_2, c2gridn,'spline');
[X,Y] = meshgrid(omega_1n,omega_2n);
omegan = X.*Y;
omega_1_initialn = interp1(c1grid, omega_1_initial, c1gridn,'spline');
omega_2_initialn = interp1(c2grid, omega_2_initial, c2gridn,'spline');
[X,Y] = meshgrid(omega_1_initialn,omega_2_initialn);
omega_initial_n = X.*Y;
figure('visible','on')
surf(c1gridn, c2gridn, omega_initial_n, 'DisplayName', 'Initial Guess, $\omega^{(0)}(C_1,C_2)$', 'EdgeColor', "#A2142F", 'FaceColor', "#A2142F", 'FaceAlpha',0.5)
hold on
surf(c1gridn, c2gridn, omegan, 'DisplayName', 'Reconstructed, $\bar{\omega}(C_1,C_2)$', 'EdgeColor', "#0072BD", 'FaceColor', "#0072BD", 'FaceAlpha',0.5)
hold off
legend('interpreter','latex','Location','best','FontSize',15)
xlabel('$C_1$','interpreter','latex','FontSize',15)
ylabel('$C_2$','interpreter','latex','FontSize',15)
zlabel('$\omega(C_1,C_2)$','interpreter','latex','FontSize',15)
loc = append(pwd,foldername,'/_omegas.png');
saveas(gcf, loc);



%  plotting lambda history 
figure('visible','on')
scatter(1:max(size(lambda_history)),lambda_history,[],"blue")
xlabel('Iteration No.','interpreter','latex','FontSize',15)
ylabel('$\alpha$','interpreter','latex','FontSize',15)
loc = append(pwd,foldername,'/_lambda.png');
saveas(gcf, loc);
end