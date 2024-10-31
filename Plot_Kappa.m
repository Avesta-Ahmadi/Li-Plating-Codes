function [] = Plot_Kappa(epsilon, k1, k2, k3)
%%% -----------------------------------------------------------------
% function to plot kappa test results. 
% It recieve kappa test vectors and plots all.


figure() % for kappa_1
semilogx(epsilon, log10(abs(k1(1,:)-1)),'LineWidth',2,'LineStyle', '--', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\omega_1^\prime = (C+0.1)^3$, $N=100$')
hold on
semilogx(epsilon, log10(abs(k1(4,:)-1)),'LineWidth',2,'LineStyle', '-', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\omega_1^\prime = (C+0.1)^3$, $N=5000$')
%semilogx(epsilon, log10(abs(k1(2,:)-1)),'LineWidth',2,'Marker','o','MarkerSize',7,'DisplayName','$\omega_1 = 0.1C^2$, $N=100$')
%semilogx(epsilon, log10(abs(k1(5,:)-1)),'LineWidth',2,'Marker','o','MarkerSize',7,'DisplayName','$\omega_1 = 0.1C^2$, $N=5000$')
semilogx(epsilon, log10(abs(k1(3,:)-1)),'LineWidth',2,'LineStyle', '--', 'Color', [0.8500 0.3250 0.0980],'Marker','o','MarkerSize',7,'DisplayName','$\omega_1^\prime = C^2$, $N=100$')
semilogx(epsilon, log10(abs(k1(6,:)-1)),'LineWidth',2,'LineStyle', '-', 'Color', [0.8500 0.3250 0.0980],'Marker','o','MarkerSize',7,'DisplayName','$\omega_1^\prime = C^2$, $N=5000$')
xlabel('$\epsilon$','interpreter','latex','FontSize',15)
ylabel('$\log(|\kappa_1(\epsilon)-1|)$','interpreter','latex','FontSize',15)
xlim([1e-16,1])
legend('interpreter','latex','FontSize',15)
grid on
hold off

figure() % for kappa_2
semilogx(epsilon, log10(abs(k2(1,:)-1)),'LineWidth',2,'LineStyle', '--', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\omega_2^\prime = (C+0.1)^3$, $N=100$')
hold on
semilogx(epsilon, log10(abs(k2(4,:)-1)),'LineWidth',2,'LineStyle', '-', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\omega_2^\prime = (C+0.1)^3$, $N=5000$')
%semilogx(epsilon, log10(abs(k2(2,:)-1)),'LineWidth',2,'Marker','o','MarkerSize',7,'DisplayName','$\omega_2 = 0.1C^2$, $N=100$')
%semilogx(epsilon, log10(abs(k2(5,:)-1)),'LineWidth',2,'Marker','o','MarkerSize',7,'DisplayName','$\omega_2 = 0.1C^2$, $N=5000$')
semilogx(epsilon, log10(abs(k2(3,:)-1)),'LineWidth',2,'LineStyle', '--', 'Color', [0.8500 0.3250 0.0980],'Marker','o','MarkerSize',7,'DisplayName','$\omega_2^\prime = C^2$, $N=100$')
semilogx(epsilon, log10(abs(k2(6,:)-1)),'LineWidth',2,'LineStyle', '-', 'Color', [0.8500 0.3250 0.0980],'Marker','o','MarkerSize',7,'DisplayName','$\omega_2^\prime = C^2$, $N=5000$')
xlabel('$\epsilon$','interpreter','latex','FontSize',15)
ylabel('$\log(|\kappa_2(\epsilon)-1|)$','interpreter','latex','FontSize',15)
xlim([1e-16,1])
legend('interpreter','latex','FontSize',15)
grid on
hold off

figure() % for kappa_3
semilogx(epsilon, log10(abs(k3(1,:)-1)),'LineWidth',2,'LineStyle', '--', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\alpha^\prime = 1$, $N=100$')
hold on
semilogx(epsilon, log10(abs(k3(4,:)-1)),'LineWidth',2,'LineStyle', '-', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\alpha^\prime = 1$, $N=5000$')
xlabel('$\epsilon$','interpreter','latex','FontSize',15)
ylabel('$\log(|\kappa_3(\epsilon)-1|)$','interpreter','latex','FontSize',15)
xlim([1e-16,1])
legend('interpreter','latex','FontSize',15)
grid on
hold off

% figure() % for kappa_4
% semilogx(epsilon, log10(abs(k4(1,:)-1)),'LineWidth',2,'LineStyle', '--', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\lambda^\prime = 1$, $N=100$')
% hold on
% semilogx(epsilon, log10(abs(k4(4,:)-1)),'LineWidth',2,'LineStyle', '-', 'Color', [0 0.4470 0.7410],'Marker','.','MarkerSize',24,'DisplayName','$\lambda^\prime = 1$, $N=5000$')
% xlabel('$\epsilon$','interpreter','latex','FontSize',15)
% ylabel('$\log(|\kappa_4(\epsilon)-1|)$','interpreter','latex','FontSize',15)
% xlim([1e-16,1])
% legend('interpreter','latex','FontSize',15)
% grid on
% hold off
end
