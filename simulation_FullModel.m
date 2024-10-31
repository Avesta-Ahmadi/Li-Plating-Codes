% Script to perform simulation on manufacture experimental data
clear; clc; close all;
data = load('data.mat');
load('sdata.mat');

global timeseq I C0
global c1tilde c2tilde c1grid c2grid c1alpha c1beta c2alpha c2beta dC1 dC2 Ngrid Nt Feval beta weight

% determine constants / cell specifications
c1alpha = -0.2;                                % lower bound on C_1 in \L interval
c1beta  = 1.2;                                % upper bound on C_1 in \L interval
c2alpha = -0.2;                                % lower bound on C_2 in \L interval
c2beta  = 1.2;                                % upper bound on C_2 in \L interval
lambda = 5;                                 % parameter of the model
weight = 50;

field = 'R2';
c1tilde = getfield(sdata, field).c1(getfield(sdata, field).label == 'charge');
c2tilde = getfield(sdata, field).c2(getfield(sdata, field).label == 'charge');
timeseq = getfield(sdata, field).t(getfield(sdata, field).label == 'charge');
current = getfield(sdata, field).J(getfield(sdata, field).label == 'charge');
pol = polyfit(timeseq,current,4);
I = polyval(pol,timeseq);

C0 = [c1tilde(1), c2tilde(1)];
Nt = max(size(timeseq));

% Manufacture constitutive relations and experimental data
Ngrid = 1000;                               % discretization of concentration for plotting
Nt = 100;                                   % discretization of time for plotting

c1grid = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
c2grid = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;

% assume omega
manufactured_omega_1 = ones(1,Ngrid)*0.9 + 0.05*sin(c1grid*3) + 0.05*cos(c1grid*3);
manufactured_omega_2 = ones(1,Ngrid)*0.95 - 0.05*sin(c2grid*2) - 0.1*cos(c2grid*2);
[X,Y] = meshgrid(manufactured_omega_1,manufactured_omega_2);
manufactured_omega = X.*Y;


figure()
plot(c1grid,manufactured_omega_1,'LineWidth',2)
xlabel('$C_1$','interpreter','latex','FontSize',15)
ylabel('$\omega_1(C_1)$','interpreter','latex','FontSize',15)
figure()
plot(c2grid,manufactured_omega_2,'LineWidth',2)
xlabel('$C_2$','interpreter','latex','FontSize',15)
ylabel('$\omega_2(C_2)$','interpreter','latex','FontSize',15)
figure()
surf(c1grid, c2grid, manufactured_omega)
xlabel('$C_1$','interpreter','latex','FontSize',15)
ylabel('$C_2$','interpreter','latex','FontSize',15)
zlabel('$\omega(C_1,C_2)$','interpreter','latex','FontSize',15)

% assume beta
%load('Results_OCV.mat',"beta");
beta=[-0.1,-0.1,-0.1,-0.1];
% Construct manufactured data by solving forward problem -> going to be used as experimental data
[c_manufactured] = forward_FullModel(manufactured_omega_1, manufactured_omega_2, lambda);
figure()
plot(timeseq,c_manufactured(:,1),'LineWidth',2, 'DisplayName', '$\widetilde{C}_1(t)$')
hold on
plot(timeseq,c_manufactured(:,2),'LineWidth',2, 'DisplayName', '$\widetilde{C}_2(t)$')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$\widetilde{C}_i(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','FontSize',15)

figure()
plot(timeseq, I,'LineWidth',2, 'DisplayName', '$J_{app}(t)$')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$J_{app}(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','FontSize',15)


% Update experimental concentrations to be the manufactured data
c1tilde = c_manufactured(:,1);
c2tilde = c_manufactured(:,2);










%% perform kappa test for different parameters and initial guesses!
perform_kappa_test = true;
if perform_kappa_test
    N = [100, 5000];                              
    for i=1:2
        % Adjustable parameters
        Ngrid = N(i);

        % Construct experimental data and omega functions
        c1grid = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
        c2grid = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;

        % Construct initial guess of omega
        omega_1 = ones(1,Ngrid)*0.7;                  % initial guess of omega_1
        omega_2 = ones(1,Ngrid)*0.7;                  % initial guess of omega_2
        lambda_initial = 3;
        % perform kappa test
        plt = false;
        [epsilon, kappa_1_p, kappa_2_p, kappa_3_p, kappa_1_pp, kappa_2_pp, kappa_1_ppp, kappa_2_ppp] = kappa_FullModel(omega_1, omega_2, lambda_initial, plt);

        k1(3*(i-1)+1,:) = kappa_1_p;
        k1(3*(i-1)+2,:) = kappa_1_pp;
        k1(3*(i-1)+3,:) = kappa_1_ppp;
        k2(3*(i-1)+1,:) = kappa_2_p;
        k2(3*(i-1)+2,:) = kappa_2_pp;
        k2(3*(i-1)+3,:) = kappa_2_ppp;
        k3(3*(i-1)+1,:) = kappa_3_p;
    end
    Plot_Kappa(epsilon, k1, k2, k3)
end





%% Perform kappa test with specifying parameters and funcitons to plot gradients
if perform_kappa_test
    Ngrid = 5000;
    c1grid = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
    c2grid = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;

    % Construct initial guess of omega
    omega_1 = ones(1,Ngrid)*0.7;                  % initial guess of omega_1
    omega_2 = ones(1,Ngrid)*0.7;                  % initial guess of omega_2
    lambda_initial = 3;
    % perform kappa test
    plt = true;
    [epsilon, kappa_1_p, kappa_2_p, kappa_3_p, kappa_1_pp, kappa_2_pp, kappa_1_ppp, kappa_2_ppp] = kappa_FullModel(omega_1, omega_2, lambda_initial, plt);
end






%% Perform optimal reconstruction
% Adjust parameters
Ngrid = 5000; % discretization of state space
dC1 = (c1beta-c1alpha)/(Ngrid-1);
dC2 = (c2beta-c2alpha)/(Ngrid-1);
c1grid = c1alpha:dC1:c1beta;
c2grid = c2alpha:dC2:c2beta;
manufactured_omega_1 = ones(1,Ngrid)*0.9 + 0.05*sin(c1grid*3) + 0.05*cos(c1grid*3);
manufactured_omega_2 = ones(1,Ngrid)*0.95 - 0.05*sin(c2grid*2) - 0.1*cos(c2grid*2);
[X,Y] = meshgrid(manufactured_omega_1,manufactured_omega_2);
manufactured_omega = X.*Y;
Nt = max(size(timeseq));
[c_manufactured] = forward_FullModel(manufactured_omega_1, manufactured_omega_2, lambda);
c1tilde = c_manufactured(:,1);
c2tilde = c_manufactured(:,2);

% Construct initial guess of omega
omega_1_initial = ones(1,Ngrid)*0.4;                   % initial guess of omega_1
omega_2_initial = ones(1,Ngrid)*0.85;                  % initial guess of omega_2
lambda_initial = 0.1;

% Run iterative procedure
Feval = 0;
plt = false;
mkdir Kappa;
foldername = '/Kappa';
[omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history,...
    H1gradJ1_history, H1gradJ2_history, omega_1_history, omega_2_history, lambda_history, Jtau_history] = iterative_FullModel(omega_1_initial, omega_2_initial, lambda_initial, plt, foldername);
[X,Y] = meshgrid(omega_1,omega_2);
omega = X.*Y;

% post process simulation results
post_process_simulation(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, manufactured_omega_1, manufactured_omega_2, lambda_history, omega_1_history, omega_2_history, foldername)

% save variables 
save('Kappa/Results_kappa.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "omega", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "manufactured_omega_1", "manufactured_omega_2",...
    "manufactured_omega", "lambda_initial", "lambda", "lambda_history","Jtau_history","beta");

