% Script to perform simulation on manufacture experimental data
clear; clc; close all;
data = load('data.mat');
global timeseq I A L_a F a rho_1 rho_2 sf C0
global c1tilde c2tilde c1grid c2grid c1alpha c1beta c2alpha c2beta dC1 dC2 Ngrid Nt Feval weight





%% determine constants / cell specifications
A = 6.2e-4;                                 % surface area of electrode [m2] - based on kevin paper
L_a = 35e-6;                                % negative electrode width [m]
r = 1.5e-6;                                 % radius of negative particle [m] - Based on Kevin paper
F = 96485.3321;                             % s.A/mol
a = 3/r;                                    % specific surface area of anode particles [1/m] 
% J_app = @(I) I/A;                         % applied current density [A/m2] in time, I is the current in Amps
% J_tot = @(I) J_app(I)/(a*L_a);            % total current density averaged over the domain at each time
% rho_1 = @(I) (1/(r*F))*J_tot(I);          % constant to control intercalation, at each time
% rho_2 = @(I) (a/F)* J_tot(I);             % constant to control side reaction, at each time
rho_1 = @(I) I/(3*F*A*L_a);                 % constant to control intercalation, at each time
rho_2 = @(I) I/(F*A*L_a);                   % constant to control side reaction, at each time
c1alpha = 0;                                % lower bound on C_1 in \L interval
c1beta  = 1;                                % upper bound on C_1 in \L interval
c2alpha = 0;                                % lower bound on C_2 in \L interval
c2beta  = 0.5;                              % upper bound on C_2 in \L interval
sf=1;                                       % Scaling factor of the system of ODEs
lambda = 1;                                 % conversion factor of the system of ODEs, assumed 1 for the simulation!
weight = 50;                                % weight used in cost function for reidual of C_1 and C_2




%% Manufacture constitutive relations and experimental data
Ngrid = 1000;               % discretization of concentration for plotting
Nt = 100;                   % discretization of time for plotting
timeseq = transpose(linspace(0,1,Nt)); % construct time discretization
tt = data.R2.echem.time(string(data.R2.echem.steps) == "discharge");           % time at which discharge is started
ss = size(data.R2.NMR.time);
ind = ss(1) - sum(data.R2.NMR.time>tt);                                    % index to start discharge
% interpolate current I to NMR data points
[~,indices] = unique(data.R2.echem.time);
I = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),transpose(linspace(0,data.R2.NMR.time(ind),Nt)),'spline');
I(I==0) = 1e-15; % to prevent division by zero.
c1grid = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
c2grid = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;
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

% Construct manufactured data by solving forward problem -> going to be used as experimental data
C0 = [0.2, 0.02];
[c_manufactured] = forward(manufactured_omega_1, manufactured_omega_2, lambda);
figure()
plot(timeseq,c_manufactured(:,1),'LineWidth',2, 'DisplayName', '$\widetilde{C}_1(t)$')
hold on
plot(timeseq,c_manufactured(:,2),'LineWidth',2, 'DisplayName', '$\widetilde{C}_2(t)$')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$\widetilde{C}(t)$','interpreter','latex','FontSize',15)
legend('interpreter','latex','FontSize',15)

figure()
plot(timeseq, rho_1(I),'LineWidth',2, 'DisplayName', '$\rho_1(t)$')
hold on
plot(timeseq, rho_2(I),'LineWidth',2, 'DisplayName', '$\rho_2(t)$')
xlabel('$t$','interpreter','latex','FontSize',15)
ylabel('$\rho(t)$','interpreter','latex','FontSize',15)
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
        %Nt = 100;

        % Construct experimental data and omega functions
        timeseq = transpose(linspace(0,1,Nt)); % construct time discretization
        I = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),transpose(linspace(0,data.R2.NMR.time(ind),Nt)),'spline');
        I(I==0) = 1e-15; % to prevent division by zero.
        c1grid = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
        c2grid = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;
        manufactured_omega_1 = ones(1,Ngrid)*0.9 + 0.05*sin(c1grid*3) + 0.05*cos(c1grid*3);
        manufactured_omega_2 = ones(1,Ngrid)*0.95 - 0.05*sin(c2grid*2) - 0.1*cos(c2grid*2);
        C0 = [0.2, 0.02];
        lambda=1;
        [c_manufactured] = forward(manufactured_omega_1, manufactured_omega_2, lambda);
        c1tilde = c_manufactured(:,1);
        c2tilde = c_manufactured(:,2);

        % Construct initial guess of omega
        omega_1 = ones(1,Ngrid)*0.9;                  % initial guess of omega_1
        omega_2 = ones(1,Ngrid)*0.9;                  % initial guess of omega_2
        lambda=0.7;
        % perform kappa test
        plt = false;
        [epsilon, kappa_1_p, kappa_2_p, kappa_3_p, kappa_1_pp, kappa_2_pp, kappa_3_pp, kappa_1_ppp, kappa_2_ppp, kappa_3_ppp] = kappa(omega_1, omega_2, lambda, plt);

        k1(3*(i-1)+1,:) = kappa_1_p;
        k1(3*(i-1)+2,:) = kappa_1_pp;
        k1(3*(i-1)+3,:) = kappa_1_ppp;
        k2(3*(i-1)+1,:) = kappa_2_p;
        k2(3*(i-1)+2,:) = kappa_2_pp;
        k2(3*(i-1)+3,:) = kappa_2_ppp;
        k3(3*(i-1)+1,:) = kappa_3_p;
        k3(3*(i-1)+2,:) = kappa_3_pp;
        k3(3*(i-1)+3,:) = kappa_3_ppp;
    end
    Plot_Kappa(epsilon, k1, k2, k3)
end





%% Perform kappa test with specifying parameters and funcitons to plot gradients
if perform_kappa_test
    Ngrid = 5000;
    %Nt = 100;                                   
    timeseq = transpose(linspace(0,1,Nt)); % construct time discretization
    I = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),transpose(linspace(0,data.R2.NMR.time(ind),Nt)),'spline');
    I(I==0) = 1e-15; % to prevent division by zero.
    c1grid = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
    c2grid = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;
    manufactured_omega_1 = ones(1,Ngrid)*0.9 + 0.05*sin(c1grid*3) + 0.05*cos(c1grid*3);
    manufactured_omega_2 = ones(1,Ngrid)*0.95 - 0.05*sin(c2grid*2) - 0.1*cos(c2grid*2);
    C0 = [0.2, 0.02];
    lambda=1;
    [c_manufactured] = forward(manufactured_omega_1, manufactured_omega_2, lambda);
    c1tilde = c_manufactured(:,1);
    c2tilde = c_manufactured(:,2);

    % Construct initial guess of omega
    omega_1 = ones(1,Ngrid)*0.9;                  % initial guess of omega_1
    omega_2 = ones(1,Ngrid)*0.9;                  % initial guess of omega_2
    lambda=0.7;
    % perform kappa test
    plt = true;
    [epsilon, kappa_1_p, kappa_2_p, kappa_3_p, kappa_1_pp, kappa_2_pp, kappa_3_pp, kappa_1_ppp, kappa_2_ppp, kappa_3_ppp] = kappa(omega_1, omega_2, lambda, plt);
end






%% Perform optimal reconstruction
% Adjust parameters
Ngrid = 5000; % discretization of state space
dC1 = (c1beta-c1alpha)/(Ngrid-1);
dC2 = (c2beta-c2alpha)/(Ngrid-1);
Nt = 100;     % discretization of time                               
timeseq = transpose(linspace(0,1,Nt)); % construct time discretization
I = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),transpose(linspace(0,data.R2.NMR.time(ind),Nt)),'spline');
I(I==0) = 1e-15; % to prevent division by zero.
c1grid = c1alpha:dC1:c1beta;
c2grid = c2alpha:dC2:c2beta;
manufactured_omega_1 = ones(1,Ngrid)*0.9 + 0.05*sin(c1grid*3) + 0.05*cos(c1grid*3);
manufactured_omega_2 = ones(1,Ngrid)*0.95 - 0.05*sin(c2grid*2) - 0.1*cos(c2grid*2);
[X,Y] = meshgrid(manufactured_omega_1,manufactured_omega_2);
manufactured_omega = X.*Y;
C0 = [0.2, 0.02];
lambda=1;
[c_manufactured] = forward(manufactured_omega_1, manufactured_omega_2, lambda);
c1tilde = c_manufactured(:,1);
c2tilde = c_manufactured(:,2);

% Construct initial guess of omega
omega_1_initial = ones(1,Ngrid)*0.5;                  % initial guess of omega_1
omega_2_initial = ones(1,Ngrid)*0.85;                  % initial guess of omega_2
lambda_initial = 0.1;

% Run iterative procedure
Feval = 0;
plt = false;
mkdir figures12;
foldername = '/figures12';
[omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history,...
    H1gradJ1_history, H1gradJ2_history, omega_1_history, omega_2_history, lambda_history] = iterative(omega_1_initial, omega_2_initial, lambda_initial, plt, foldername);
[X,Y] = meshgrid(omega_1,omega_2);
omega = X.*Y;

% post process simulation results
post_process_simulation(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, manufactured_omega_1, manufactured_omega_2, lambda_history, omega_1_history, omega_2_history, foldername)

% save variables 
save('results.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "omega", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "manufactured_omega_1", "manufactured_omega_2",...
    "manufactured_omega", "lambda_initial", "lambda", "lambda_history");

