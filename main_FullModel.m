% Main function to perform analysis on the whole cycle, assuming the full
% model, where we have a linear part for OCV (predefined) and a nonlinear
% part for omega
% parameters of the linear model are already determined, and will be used
% for this analysis.
% data is split into three phases: charge, OCV, discharge
% the charge, or discharge phases of data are used for training.
% note that charge phase or discharge phase are used separately for
% training, and two different models will be used for charge or discharge.
% only one cycle is used for training at a time.

clear; clc; close all;

global c1grid c2grid c1alpha c1beta c2alpha c2beta Ngrid Nt dC1 dC2 beta
global c1tilde c2tilde timeseq I C0 weight

weight = 50;
% load data
data = load('data.mat');
load('sdata.mat');

% load optimal parameters of OCV part
load('Results_OCV.mat', "beta");



%% fitting to CHARGE phase of data
% running iterative algorithm for finding optimal constitutive relations using charge portion of data
c1alpha = -0.5;
c1beta  = 1.5;
c2alpha = -0.2;
c2beta  = 0.5;
Ngrid   = 5000;
dC1 = (c1beta-c1alpha)/(Ngrid-1);
dC2 = (c2beta-c2alpha)/(Ngrid-1);
c1grid  = c1alpha:dC1:c1beta;
c2grid  = c2alpha:dC2:c2beta;


field = 'R6';
c1tilde = getfield(sdata, field).c1(getfield(sdata, field).label == 'charge');
c2tilde = getfield(sdata, field).c2(getfield(sdata, field).label == 'charge');
timeseq = getfield(sdata, field).t(getfield(sdata, field).label == 'charge');
I = getfield(sdata, field).J(getfield(sdata, field).label == 'charge');
C0 = [c1tilde(1), c2tilde(1)];
Nt = max(size(timeseq));

% Run iterative procedure
Feval = 0;
plt = false;
mkdir Charge_results/R6;
foldername = '/Charge_results/R6';

omega_1_initial = ones(1,Ngrid)* 0.25;
omega_2_initial = ones(1,Ngrid)* 0.25;
lambda_initial = 5;

[omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history, H1gradJ1_history, H1gradJ2_history, ...
    omega_1_history, omega_2_history, lambda_history, Jtau_history] = iterative_FullModel(omega_1_initial, omega_2_initial, lambda_initial, plt, foldername);

post_process_reconstruction(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, lambda_history, foldername)

% save variables
save('Charge_results/R6/Results_Charge_R6.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", ...
    "lambda_initial", "lambda", "lambda_history", "Jtau_history", "beta");



%% fitting to DISCHARGE phase of data
% running iterative algorithm for finding optimal constitutive relations using charge portion of data
c1alpha = -0.5;
c1beta  = 1.5;
c2alpha = -0.2;
c2beta  = 0.5;
Ngrid   = 5000;
c1grid  = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
c2grid  = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;
dC1 = (c1beta-c1alpha)/(Ngrid-1);
dC2 = (c2beta-c2alpha)/(Ngrid-1);

field = 'R6';
c1tilde = getfield(sdata, field).c1(getfield(sdata, field).label == 'discharge');
c2tilde = getfield(sdata, field).c2(getfield(sdata, field).label == 'discharge');
timeseq = getfield(sdata, field).t(getfield(sdata, field).label == 'discharge');
I = getfield(sdata, field).J(getfield(sdata, field).label == 'discharge');
C0 = [c1tilde(1), c2tilde(1)];
Nt = max(size(timeseq));

% Run iterative procedure
Feval = 0;
plt = false;
mkdir Discharge_results/R6;
foldername = '/Discharge_results/R6';

omega_1_initial = ones(1,Ngrid)*0.2;
omega_2_initial = ones(1,Ngrid)*0.2;
lambda_initial = 4;

[omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history, H1gradJ1_history, H1gradJ2_history, ...
    omega_1_history, omega_2_history, lambda_history, Jtau_history] = iterative_FullModel(omega_1_initial, omega_2_initial, lambda_initial, plt, foldername);

post_process_reconstruction(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, lambda_history, foldername)

% save variables
save('Discharge_results/R6/Results_Discharge_R6.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", ...
    "lambda_initial", "lambda", "lambda_history", "Jtau_history", "beta");


