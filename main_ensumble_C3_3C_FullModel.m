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

global c1grid c2grid c1alpha c1beta c2alpha c2beta Ngrid A L_a F a Nt dC1 dC2 weight beta
global c1tilde c2tilde timeseq I C0


% load data
data = load('data.mat');
load('sdata.mat');


% load optimal parameters of OCV part
load('Results_OCV.mat', "beta");
weight = 50;

%% fitting to CHARGE phase of data
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

c1tilde = struct;
c2tilde = struct;
timeseq = struct;
I = struct;
C0 = struct;
Nt = struct;

c1tilde.R2 = sdata.R2.c1(sdata.R2.label == 'charge');
% c1tilde.R3 = sdata.R3.c1(sdata.R3.label == 'charge');
% c1tilde.R4 = sdata.R4.c1(sdata.R4.label == 'charge');
% c1tilde.R5 = sdata.R5.c1(sdata.R5.label == 'charge');
c1tilde.R6 = sdata.R6.c1(sdata.R6.label == 'charge');

c2tilde.R2 = sdata.R2.c2(sdata.R2.label == 'charge');
% c2tilde.R3 = sdata.R3.c2(sdata.R3.label == 'charge');
% c2tilde.R4 = sdata.R4.c2(sdata.R4.label == 'charge');
% c2tilde.R5 = sdata.R5.c2(sdata.R5.label == 'charge');
c2tilde.R6 = sdata.R6.c2(sdata.R6.label == 'charge');

timeseq.R2 = sdata.R2.t(sdata.R2.label == 'charge');
% timeseq.R3 = sdata.R3.t(sdata.R3.label == 'charge');
% timeseq.R4 = sdata.R4.t(sdata.R4.label == 'charge');
% timeseq.R5 = sdata.R5.t(sdata.R5.label == 'charge');
timeseq.R6 = sdata.R6.t(sdata.R6.label == 'charge');

I.R2 = sdata.R2.J(sdata.R2.label == 'charge');
% I.R3 = sdata.R3.J(sdata.R3.label == 'charge');
% I.R4 = sdata.R4.J(sdata.R4.label == 'charge');
% I.R5 = sdata.R5.J(sdata.R5.label == 'charge');
I.R6 = sdata.R6.J(sdata.R6.label == 'charge');

C0.R2 = [c1tilde.R2(1), c2tilde.R2(1)];
% C0.R3 = [c1tilde.R3(1), c2tilde.R3(1)];
% C0.R4 = [c1tilde.R4(1), c2tilde.R4(1)];
% C0.R5 = [c1tilde.R5(1), c2tilde.R5(1)];
C0.R6 = [c1tilde.R6(1), c2tilde.R6(1)];

Nt.R2 = max(size(timeseq.R2));
% Nt.R3 = max(size(timeseq.R3));
% Nt.R4 = max(size(timeseq.R4));
% Nt.R5 = max(size(timeseq.R5));
Nt.R6 = max(size(timeseq.R6));


% Run iterative procedure
Feval = 0;
plt = false;
mkdir Charge_results/ensumble_C3_3C;
foldername = '/Charge_results/ensumble_C3_3C';

omega_1_initial = ones(1,Ngrid)*0.2;
omega_2_initial = ones(1,Ngrid)*0.2;
lambda_initial = 5;

[omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history, H1gradJ1_history, H1gradJ2_history, ...
    omega_1_history, omega_2_history, lambda_history, Jtau_history] = iterative_ensumble_C3_3C_FullModel(omega_1_initial, omega_2_initial, lambda_initial, plt, foldername);

post_process_reconstruction_ensumble_C3_3C(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, lambda_history, foldername)

% save variables
% save variables
save('Charge_results/ensumble_C3_3C/Results_Charge_ensumble_C3_3C.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", "cm", "tm", ...
    "A", "L_a", "F", "a","weight", "lambda_initial", "lambda", "lambda_history", "Jtau_history", "beta");




%% fitting to DISCHARGE phase of data
% running iterative algorithm for finding optimal constitutive relations using discharge portion of data
c1alpha = -0.5;
c1beta  = 1.5;
c2alpha = -0.2;
c2beta  = 0.5;
Ngrid   = 5000;
c1grid  = c1alpha:(c1beta-c1alpha)/(Ngrid-1):c1beta;
c2grid  = c2alpha:(c2beta-c2alpha)/(Ngrid-1):c2beta;
dC1 = (c1beta-c1alpha)/(Ngrid-1);
dC2 = (c2beta-c2alpha)/(Ngrid-1);

% get intial values from predictions 
cycles = ['C3 Cycle';'3C Cycle'];
fields = ['R2'; 'R6'];
initial_values = zeros(2,2);
for i = 1:max(size(fields))
    
    % charge
    load('Results_Charge_ensumble_C3_3C.mat',"omega_1", "omega_2","lambda");
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
        
    initial_values(i,:) = c_OCV(end,:);
    
    

end





c1tilde = struct;
c2tilde = struct;
timeseq = struct;
I = struct;
C0 = struct;
Nt = struct;




c1tilde.R2 = sdata.R2.c1(sdata.R2.label == 'discharge');
% c1tilde.R3 = sdata.R3.c1(sdata.R3.label == 'discharge');
% c1tilde.R4 = sdata.R4.c1(sdata.R4.label == 'discharge');
% c1tilde.R5 = sdata.R5.c1(sdata.R5.label == 'discharge');
c1tilde.R6 = sdata.R6.c1(sdata.R6.label == 'discharge');

c2tilde.R2 = sdata.R2.c2(sdata.R2.label == 'discharge');
% c2tilde.R3 = sdata.R3.c2(sdata.R3.label == 'discharge');
% c2tilde.R4 = sdata.R4.c2(sdata.R4.label == 'discharge');
% c2tilde.R5 = sdata.R5.c2(sdata.R5.label == 'discharge');
c2tilde.R6 = sdata.R6.c2(sdata.R6.label == 'discharge');

timeseq.R2 = sdata.R2.t(sdata.R2.label == 'discharge');
% timeseq.R3 = sdata.R3.t(sdata.R3.label == 'discharge');
% timeseq.R4 = sdata.R4.t(sdata.R4.label == 'discharge');
% timeseq.R5 = sdata.R5.t(sdata.R5.label == 'discharge');
timeseq.R6 = sdata.R6.t(sdata.R6.label == 'discharge');

I.R2 = sdata.R2.J(sdata.R2.label == 'discharge');
% I.R3 = sdata.R3.J(sdata.R3.label == 'discharge');
% I.R4 = sdata.R4.J(sdata.R4.label == 'discharge');
% I.R5 = sdata.R5.J(sdata.R5.label == 'discharge');
I.R6 = sdata.R6.J(sdata.R6.label == 'discharge');

% C0.R2 = [c1tilde.R2(1), c2tilde.R2(1)];
% C0.R3 = [c1tilde.R3(1), c2tilde.R3(1)];
% C0.R4 = [c1tilde.R4(1), c2tilde.R4(1)];
% C0.R5 = [c1tilde.R5(1), c2tilde.R5(1)];
% C0.R6 = [c1tilde.R6(1), c2tilde.R6(1)];
C0.R2 = initial_values(1,:);
% C0.R3 = initial_values(2,:);
% C0.R4 = initial_values(3,:);
% C0.R5 = initial_values(4,:);
C0.R6 = initial_values(2,:);

Nt.R2 = max(size(timeseq.R2));
% Nt.R3 = max(size(timeseq.R3));
% Nt.R4 = max(size(timeseq.R4));
% Nt.R5 = max(size(timeseq.R5));
Nt.R6 = max(size(timeseq.R6));


% Run iterative procedure
Feval = 0;
plt = false;
mkdir Discharge_results/ensumble_C3_3C;
foldername = '/Discharge_results/ensumble_C3_3C';

omega_1_initial = ones(1,Ngrid)*0.2;
omega_2_initial = ones(1,Ngrid)*0.2;
lambda_initial = 5;

[omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history, H1gradJ1_history, H1gradJ2_history, ...
    omega_1_history, omega_2_history, lambda_history, Jtau_history] = iterative_ensumble_C3_3C_FullModel(omega_1_initial, omega_2_initial, lambda_initial, plt, foldername);

post_process_reconstruction_ensumble_C3_3C(omega_1, omega_2, c1forward_history, c2forward_history, omega_1_initial, omega_2_initial, lambda_history, foldername)

% save variables
% save variables
save('Discharge_results/ensumble_C3_3C/Results_Discharge_ensumble_C3_3C.mat',"C0","c1alpha","c1beta","c2alpha","c2beta","c1grid","c2grid","c1tilde","c2tilde","dC1","dC2","Ngrid","Nt",...
    "omega_1_initial","omega_2_initial","omega_1", "omega_2", "c1forward_history", "c2forward_history", "L2gradJ1_history", ...
    "L2gradJ2_history","H1gradJ1_history","H1gradJ2_history", "omega_1_history", "omega_2_history", "sdata", "cm", "tm", ...
    "A", "L_a", "F", "a", "weight", "lambda_initial", "lambda", "lambda_history","Jtau_history", "beta");

