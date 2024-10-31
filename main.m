% main function to perform analysis
clear; clc; close all;

% data = preprocess();
data = load('data.mat');
% plot_raw_data(data)

global c1grid c2grid c1alpha c1beta c2alpha c2beta Ngrid A L_a F a rho_1 rho_2 sf Nt dC1 dC2 weight
global c1tilde c2tilde timeseq I C0

%% determine constants
% Infromation to be adjusted based on cell specifications
A = 6.2e-4;                                 % surface area of electrode [m2] - based on kevin paper
L_a = 35e-6;                                % negative electrode width [m]
r = 1.5e-6;                                 % radius of negative particle [m] - Based on Kevin paper
F = 96485.3321;                             % s.A/mol
a = 3/r;                                    % specific surface area of anode particles [1/m] 
rho_1 = @(I) I/(3*F*A*L_a);                 % constant to control intercalation, at each time
rho_2 = @(I) I/(F*A*L_a);                   % constant to control side reaction, at each time
weight = 50;                                % weight used in cost function for reidual of C_1 and C_2



%% prepare data and store in sdata
% compute C_1 and C_2 for each cycle, take time zero to start of discharge.
steps = ["OCV pre charge"; "charge"; "OCV post charge"; "discharge"; "OCV post discharge"];
sdata = {};


% C20 cycle
tt = data.R1.echem.time(string(data.R1.echem.steps) == "discharge");       
ind_discharge = size(data.R1.NMR.time,1) - sum(data.R1.NMR.time>tt);
c1 = data.R1.NMR.phase1(1:ind_discharge) + data.R1.NMR.phase2(1:ind_discharge);
t = data.R1.NMR.time(1:ind_discharge)*3600;
[~,indices] = unique(data.R1.echem.time);
J = interp1(data.R1.echem.time(indices), data.R1.echem.I(indices),data.R1.NMR.time(1:ind_discharge),'spline');
J(J==0) = 1e-15; 
sdata.R1 = table(t,J,c1);

% C3 cycle
tt = data.R2.echem.time(string(data.R2.echem.steps) == "discharge");       
ind_discharge = size(data.R2.NMR.time,1) - sum(data.R2.NMR.time>tt);
data.R2.NMR.phase3(isnan(data.R2.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R2.NMR.phase1(1:ind_discharge) + data.R2.NMR.phase2(1:ind_discharge) + data.R2.NMR.phase3(1:ind_discharge);
c2 = data.R2.NMR.plated(1:ind_discharge);
t = data.R2.NMR.time(1:ind_discharge)*3600;
[~,indices] = unique(data.R2.echem.time);
J = interp1(data.R2.echem.time(indices), data.R2.echem.I(indices),data.R2.NMR.time(1:ind_discharge),'spline');
J(J==0) = 1e-15; % to prevent division by zero.
sdata.R2 = table(t,J,c1,c2);

% C2 cycle
tt = data.R3.echem.time(string(data.R3.echem.steps) == "discharge");       
ind_discharge = size(data.R3.NMR.time,1) - sum(data.R3.NMR.time>tt);
data.R3.NMR.phase3(isnan(data.R3.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R3.NMR.phase1(1:ind_discharge) + data.R3.NMR.phase2(1:ind_discharge) + data.R3.NMR.phase3(1:ind_discharge);
c2 = data.R3.NMR.plated(1:ind_discharge);
t = data.R3.NMR.time(1:ind_discharge)*3600;
[~,indices] = unique(data.R3.echem.time);
J = interp1(data.R3.echem.time(indices), data.R3.echem.I(indices),data.R3.NMR.time(1:ind_discharge),'spline');
J(J==0) = 1e-15; % to prevent division by zero.
sdata.R3 = table(t,J,c1,c2);

% 1C cycle
tt = data.R4.echem.time(string(data.R4.echem.steps) == "discharge");       
ind_discharge = size(data.R4.NMR.time,1) - sum(data.R4.NMR.time>tt);
data.R4.NMR.phase3(isnan(data.R4.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R4.NMR.phase1(1:ind_discharge) + data.R4.NMR.phase2(1:ind_discharge) + data.R4.NMR.phase3(1:ind_discharge);
c2 = data.R4.NMR.plated(1:ind_discharge) + data.R4.NMR.dendrite(1:ind_discharge);
t = data.R4.NMR.time(1:ind_discharge)*3600;
[~,indices] = unique(data.R4.echem.time);
J = interp1(data.R4.echem.time(indices), data.R4.echem.I(indices),data.R4.NMR.time(1:ind_discharge),'spline');
J(J==0) = 1e-15; % to prevent division by zero.
sdata.R4 = table(t,J,c1,c2);

% 2C cycle
tt = data.R5.echem.time(string(data.R5.echem.steps) == "discharge");       
ind_discharge = size(data.R5.NMR.time,1) - sum(data.R5.NMR.time>tt);
data.R5.NMR.phase3(isnan(data.R5.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R5.NMR.phase1(1:ind_discharge) + data.R5.NMR.phase2(1:ind_discharge) + data.R5.NMR.phase3(1:ind_discharge);
c2 = data.R5.NMR.plated(1:ind_discharge) + data.R5.NMR.dendrite(1:ind_discharge);
t = data.R5.NMR.time(1:ind_discharge)*3600;
[~,indices] = unique(data.R5.echem.time);
J = interp1(data.R5.echem.time(indices), data.R5.echem.I(indices),data.R5.NMR.time(1:ind_discharge),'spline');
J(J==0) = 1e-15; % to prevent division by zero.
sdata.R5 = table(t,J,c1,c2);

% 3C cycle
tt = data.R6.echem.time(string(data.R6.echem.steps) == "discharge");       
ind_discharge = size(data.R6.NMR.time,1) - sum(data.R6.NMR.time>tt);
data.R6.NMR.phase3(isnan(data.R6.NMR.phase3)) = 0;                         % replace NAN with zero
c1 = data.R6.NMR.phase1(1:ind_discharge) + data.R6.NMR.phase2(1:ind_discharge) + data.R6.NMR.phase3(1:ind_discharge);
c2 = data.R6.NMR.plated(1:ind_discharge) + data.R6.NMR.dendrite(1:ind_discharge);
t = data.R6.NMR.time(1:ind_discharge)*3600;
[~,indices] = unique(data.R6.echem.time);
J = interp1(data.R6.echem.time(indices), data.R6.echem.I(indices),data.R6.NMR.time(1:ind_discharge),'spline');
J(J==0) = 1e-15; % to prevent division by zero.
sdata.R6 = table(t,J,c1,c2);



% Normalize variables, rescale ODE system, define tm, cm, lambda, sf
cm = max(sdata.R2.c1 + sdata.R2.c2);
tm = max(sdata.R2.t);
sf = tm/cm;
%lambda = 37992071.25;
%lambda = 29725328.737500;

%sdata.R1.t  = sdata.R1.t/tm;
%sdata.R1.c1 = sdata.R1.c1/cm;
sdata.R2.t  = sdata.R2.t/tm;
sdata.R2.c1 = sdata.R2.c1/cm;
sdata.R2.c2 = sdata.R2.c2/cm;
sdata.R3.t  = sdata.R3.t/tm;
sdata.R3.c1 = sdata.R3.c1/cm;
sdata.R3.c2 = sdata.R3.c2/cm;
sdata.R4.t  = sdata.R4.t/tm;
sdata.R4.c1 = sdata.R4.c1/cm;
sdata.R4.c2 = sdata.R4.c2/cm;
sdata.R5.t  = sdata.R5.t/tm;
sdata.R5.c1 = sdata.R5.c1/cm;
sdata.R5.c2 = sdata.R5.c2/cm;
sdata.R6.t  = sdata.R6.t/tm;
sdata.R6.c1 = sdata.R6.c1/cm;
sdata.R6.c2 = sdata.R6.c2/cm;


% plot data
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


%% State descriziation parameters
c1alpha = -0.5;
c1beta  = 1.5;
c2alpha = -0.2;
c2beta  = 0.5;
Ngrid   = 5000;
c1grid  = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
c2grid  = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;


%% run forward model for testing!
% omega_1 = ones(1,Ngrid)*0.9;
% omega_2 = ones(1,Ngrid)*0.9;
% global C0
% timeseq = sdata.R2.t;
% I = sdata.R2.J;
% C0 = [sdata.R2.c1(1), sdata.R2.c2(1)];
% c = forward(omega_1, omega_2);
% figure()
% plot(c(:,1))
% hold on 
% plot(sdata.R2.c1)
% figure()
% plot(c(:,2))
% hold on
% plot(sdata.R2.c2)



%% perform kappa test on C3 cycle
% can be done for other cycles as well!
perform_kappa_test = false;
if perform_kappa_test
    c1tilde = sdata.R2.c1;
    c2tilde = sdata.R2.c2;
    timeseq = sdata.R2.t;
    I = sdata.R2.J;
    C0 = [sdata.R2.c1(1), sdata.R2.c2(1)];
    plt = true;
    N = [100, 5000];
    for i = 1:2
        Ngrid = N(i);
        c1grid = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
        c2grid = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;
        omega_1 = ones(1,Ngrid)*0.995;                  % initial guess of omega_1
        omega_2 = ones(1,Ngrid)*0.995;                  % initial guess of omega_2
        lambda = 2.9e7;
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


%% perform optimal reconstruction on C3 cycle
% adjust parameters and data!
% set training data
traindataset = 'R6';
Ngrid = 5000;
dC1 = (c1beta-c1alpha)/(Ngrid-1);
dC2 = (c2beta-c2alpha)/(Ngrid-1);
c1grid = c1alpha:dC1:c1beta;
c2grid = c2alpha:dC2:c2beta;
c1tilde = getfield(sdata, traindataset).c1;
c2tilde = getfield(sdata, traindataset).c2;
timeseq = getfield(sdata, traindataset).t;
I = getfield(sdata, traindataset).J;
C0 = [getfield(sdata, traindataset).c1(1), getfield(sdata, traindataset).c2(1)];
Nt = max(size(timeseq));

% initial guesses
omega_1_initial = ones(1,Ngrid)*0.95;
omega_2_initial = ones(1,Ngrid)*0.99;
lambda_initial = 2.7e7;

% Run iterative procedure
Feval = 0;
plt = false;
mkdir Results_R6;
foldername = '/Results_R6';
[omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history, H1gradJ1_history, H1gradJ2_history,...
    omega_1_history, omega_2_history, lambda_history] = iterative(omega_1_initial, omega_2_initial, lambda_initial, plt, foldername);
[X,Y] = meshgrid(omega_1,omega_2);
omega = X.*Y;

% post process simulation results
post_process_reconstruction(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, lambda_history, foldername)

% save variables
save('Results_R6.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "omega", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", "cm", "tm", "sf", ...
    "A", "L_a", "r", "F", "a", "rho_1","rho_2","weight", "lambda_initial", "lambda", "lambda_history");
