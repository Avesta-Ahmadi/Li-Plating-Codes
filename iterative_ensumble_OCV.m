function [beta, c1forward_history, c2forward_history, pjpb_history, beta_history, Jtau_history] = iterative_ensumble_OCV(beta, plt, foldername)
%%% -----------------------------------------------------------------
% the iterative scheme of minimization only for the OCV part of the model
% and data, to fit a linear model to OCV part, as dC/dt = AC, this code is
% an iterative gradient descent scheme that minimezes the cost function
% with respect to coefficients of the linear model, using an adjoint
% formulation for finding the partial derivatives.
% Input:     experimental data in time, for the OCV part of data, for all cycles
%            beta ---> intial guess of coefficients
% Output:    beta --> optimalvalues of coefficients



global timeseq Nt C0 c1tilde c2tilde

% write info in a file
fout = fopen('info.txt','w');
fprintf(fout,'Parameters of the simulation are as follows: \n');

%%% OPTIMIZATION PARAMETERS
MAX_ITER       = 500;       % Maximum number of iterations of algorithm;
MAXITER        = 100;       % Maximum number of bracketing iterations;
ALPHA0         = 1 ;        % Initial interval for bracketing;
Tol_J          = 1.0e-6;    % Tolerance of iterative algorithm.
MAX_ITER_BRENT = 50;        % Max iterations for Brent routine.
Tol_Brent      = 1.0e-4;    % Tolerance in Brent routine.
MOMENTUM       = 0.0;       % Momentum for conjugate gradient algorithm
G2_1           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for beta_1
G2_2           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for beta_2
G2_3           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for beta_3
G2_4           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for beta_4
%G2_5           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for beta_5
%G2_6           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for beta_6
GRADIENT       = 3;         % type of conjugate gradient (1 for simple gradient, 2 for Fletcher-Reeves CG, and 3 for Polak-Ribiere CG)


fprintf(fout,'     Maximum number of optimal reconstruction iterations = %6.2f \n',MAX_ITER);
fprintf(fout,'     Maximum number of bracketing iterations = %6.2f \n',MAXITER);
fprintf(fout,'     Maximum number of brent routine iterations = %6.2f \n',MAX_ITER_BRENT);
fprintf(fout,'     Initial interval for bracketing = %6.2f \n',ALPHA0);
fprintf(fout,'     Tolerance of iterative algorithm = %6.2f \n',Tol_J);
fprintf(fout,'     Tolerance of brent algorithm = %6.2f \n',Tol_Brent);
fprintf(fout,'     Type pf conjugate gradient used = %6.2f \n',GRADIENT);


% Store intermediate solutions in a table, rows correspond to iteration number
c1forward_history = struct;
c2forward_history = struct;
c1adjoint_history = struct;
c2adjoint_history = struct;
c1forward_history.R2 = zeros(1,Nt.R2);      % forward solutions for C_1
c1forward_history.R3 = zeros(1,Nt.R3);      % forward solutions for C_1
c1forward_history.R4 = zeros(1,Nt.R4);      % forward solutions for C_1
c1forward_history.R5 = zeros(1,Nt.R5);      % forward solutions for C_1
c1forward_history.R6 = zeros(1,Nt.R6);      % forward solutions for C_1
c2forward_history.R2 = zeros(1,Nt.R2);      % forward solutions for C_2
c2forward_history.R3 = zeros(1,Nt.R3);      % forward solutions for C_2
c2forward_history.R4 = zeros(1,Nt.R4);      % forward solutions for C_2
c2forward_history.R5 = zeros(1,Nt.R5);      % forward solutions for C_2
c2forward_history.R6 = zeros(1,Nt.R6);      % forward solutions for C_2
c1adjoint_history.R2 = zeros(1,Nt.R2);      % adjoint solutions for C_1
c1adjoint_history.R3 = zeros(1,Nt.R3);      % adjoint solutions for C_1
c1adjoint_history.R4 = zeros(1,Nt.R4);      % adjoint solutions for C_1
c1adjoint_history.R5 = zeros(1,Nt.R5);      % adjoint solutions for C_1
c1adjoint_history.R6 = zeros(1,Nt.R6);      % adjoint solutions for C_1
c2adjoint_history.R2 = zeros(1,Nt.R2);      % adjoint solutions for C_2
c2adjoint_history.R3 = zeros(1,Nt.R3);      % adjoint solutions for C_2
c2adjoint_history.R4 = zeros(1,Nt.R4);      % adjoint solutions for C_2
c2adjoint_history.R5 = zeros(1,Nt.R5);      % adjoint solutions for C_2
c2adjoint_history.R6 = zeros(1,Nt.R6);      % adjoint solutions for C_2
beta_history = beta;                        % partial derivaive of J with respect to lambda
Jtau_history = zeros(1,4);                  % History of cost function values
pjpb_history = zeros(1,4);                  % History of partial derivatives

Jtau_beta_1  = 1.0;
Jtau_beta_2  = 1.0;
Jtau_beta_3  = 1.0;
Jtau_beta_4  = 1.0;
%Jtau_beta_5  = 1.0;
%Jtau_beta_6  = 1.0;
Jtau_beta_40 = 0.0;


% set iterator to zero
n = 0;
% iterative scheme
while (n <= MAX_ITER) && (abs(Jtau_beta_40 - Jtau_beta_4) / abs(Jtau_beta_4) > Tol_J)
    
    % Update iterator
    n = n+1;
    fprintf(fout, '   n = %d \n',n);
    fprintf('   n = %d \n',n);
    
    % Solve forawrd problem
    cforward = struct;
    cforward.R2 = forward_ensumble_OCV(beta, timeseq.R2, C0.R2);
    cforward.R3 = forward_ensumble_OCV(beta, timeseq.R3, C0.R3);
    cforward.R4 = forward_ensumble_OCV(beta, timeseq.R4, C0.R4);
    cforward.R5 = forward_ensumble_OCV(beta, timeseq.R5, C0.R5);
    cforward.R6 = forward_ensumble_OCV(beta, timeseq.R6, C0.R6);
    c1forward_history.R2(n,:) = transpose(cforward.R2(:,1));
    c1forward_history.R3(n,:) = transpose(cforward.R3(:,1));
    c1forward_history.R4(n,:) = transpose(cforward.R4(:,1));
    c1forward_history.R5(n,:) = transpose(cforward.R5(:,1));
    c1forward_history.R6(n,:) = transpose(cforward.R6(:,1));
    c2forward_history.R2(n,:) = transpose(cforward.R2(:,2));
    c2forward_history.R3(n,:) = transpose(cforward.R3(:,2));
    c2forward_history.R4(n,:) = transpose(cforward.R4(:,2));
    c2forward_history.R5(n,:) = transpose(cforward.R5(:,2));
    c2forward_history.R6(n,:) = transpose(cforward.R6(:,2));

    fprintf(fout,'   Forward problem solved ... \n');
    
    % Solve adjoint problem
    cstar = struct;
    cstar.R2 = adjoint_ensumble_OCV(beta, cforward.R2, c1tilde.R2, c2tilde.R2, timeseq.R2);
    cstar.R3 = adjoint_ensumble_OCV(beta, cforward.R3, c1tilde.R3, c2tilde.R3, timeseq.R3);
    cstar.R4 = adjoint_ensumble_OCV(beta, cforward.R4, c1tilde.R4, c2tilde.R4, timeseq.R4);
    cstar.R5 = adjoint_ensumble_OCV(beta, cforward.R5, c1tilde.R5, c2tilde.R5, timeseq.R5);
    cstar.R6 = adjoint_ensumble_OCV(beta, cforward.R6, c1tilde.R6, c2tilde.R6, timeseq.R6);
    c1adjoint_history.R2(n,:) = transpose(cstar.R2(:,1));
    c2adjoint_history.R2(n,:) = transpose(cstar.R2(:,2));
    c1adjoint_history.R3(n,:) = transpose(cstar.R3(:,1));
    c2adjoint_history.R3(n,:) = transpose(cstar.R3(:,2));
    c1adjoint_history.R4(n,:) = transpose(cstar.R4(:,1));
    c2adjoint_history.R4(n,:) = transpose(cstar.R4(:,2));
    c1adjoint_history.R5(n,:) = transpose(cstar.R5(:,1));
    c2adjoint_history.R5(n,:) = transpose(cstar.R5(:,2));
    c1adjoint_history.R6(n,:) = transpose(cstar.R6(:,1));
    c2adjoint_history.R6(n,:) = transpose(cstar.R6(:,2));
    if plt 
        figure('visible','off')
        hold on
        plot(timeseq, cstar(:,1),'DisplayName','$C^\ast_1(t)$','LineWidth',2)
        plot(timeseq, cstar(:,2),'DisplayName','$C^\ast_2(t)$','LineWidth',2)
        xlabel('$t$','interpreter','latex','FontSize',15)
        ylabel('Concentration','interpreter','latex','FontSize',15)
        legend('interpreter','latex','FontSize',15)
        hold off
        loc = append(pwd,foldername,'/',num2str(n),'_adjoint.png');
        saveas(gcf, loc);
    end
    fprintf(fout,'   Adjoint problem solved ... \n');
    
    
    
    % computing partial derivative of J w.r.t beta
    pjpb = zeros(1,4);
    % for beta_1
    pj = zeros(1,5);
    pj(1) = trapz(timeseq.R2, -cstar.R2(:,1) );
    pj(2) = trapz(timeseq.R3, -cstar.R3(:,1) );
    pj(3) = trapz(timeseq.R4, -cstar.R4(:,1) );
    pj(4) = trapz(timeseq.R5, -cstar.R5(:,1) );
    pj(5) = trapz(timeseq.R6, -cstar.R6(:,1) );
    pjpb(1) = sum(pj);
    % for beta_2
    pj = zeros(1,5);
    pj(1) = trapz(timeseq.R2, - cforward.R2(:,1).*cstar.R2(:,1)  );
    pj(2) = trapz(timeseq.R3, - cforward.R3(:,1).*cstar.R3(:,1)  );
    pj(3) = trapz(timeseq.R4, - cforward.R4(:,1).*cstar.R4(:,1)  );
    pj(4) = trapz(timeseq.R5, - cforward.R5(:,1).*cstar.R5(:,1)  );
    pj(5) = trapz(timeseq.R6, - cforward.R6(:,1).*cstar.R6(:,1)  );
    pjpb(2) = sum(pj);
    % for beta_3
    pj = zeros(1,5);
    pj(1) = trapz(timeseq.R2, - 0.006*cforward.R2(:,2).*cstar.R2(:,1) + cforward.R2(:,2).*cstar.R2(:,2) );
    pj(2) = trapz(timeseq.R3, - 0.006*cforward.R3(:,2).*cstar.R3(:,1) + cforward.R3(:,2).*cstar.R3(:,2) );
    pj(3) = trapz(timeseq.R4, - 0.006*cforward.R4(:,2).*cstar.R4(:,1) + cforward.R4(:,2).*cstar.R4(:,2) );
    pj(4) = trapz(timeseq.R5, - 0.006*cforward.R5(:,2).*cstar.R5(:,1) + cforward.R5(:,2).*cstar.R5(:,2) );
    pj(5) = trapz(timeseq.R6, - 0.006*cforward.R6(:,2).*cstar.R6(:,1) + cforward.R6(:,2).*cstar.R6(:,2) );
    pjpb(3) = sum(pj);
    % for beta_4
    pj = zeros(1,5);
    pj(1) = trapz(timeseq.R2, - cforward.R2(:,1).*cstar.R2(:,2)  );
    pj(2) = trapz(timeseq.R3, - cforward.R3(:,1).*cstar.R3(:,2)  );
    pj(3) = trapz(timeseq.R4, - cforward.R4(:,1).*cstar.R4(:,2)  );
    pj(4) = trapz(timeseq.R5, - cforward.R5(:,1).*cstar.R5(:,2)  );
    pj(5) = trapz(timeseq.R6, - cforward.R6(:,1).*cstar.R6(:,2)  );
    pjpb(4) = sum(pj);
%     % for beta_1
%     pj = zeros(1,5);
%     pj(1) = trapz(timeseq.R2, -cstar.R2(:,1) );
%     pj(2) = trapz(timeseq.R3, -cstar.R3(:,1) );
%     pj(3) = trapz(timeseq.R4, -cstar.R4(:,1) );
%     pj(4) = trapz(timeseq.R5, -cstar.R5(:,1) );
%     pj(5) = trapz(timeseq.R6, -cstar.R6(:,1) );
%     pjpb(1) = sum(pj);
%     % for beta_2
%     pj = zeros(1,5);
%     pj(1) = trapz(timeseq.R2, -cstar.R2(:,2) );
%     pj(2) = trapz(timeseq.R3, -cstar.R3(:,2) );
%     pj(3) = trapz(timeseq.R4, -cstar.R4(:,2) );
%     pj(4) = trapz(timeseq.R5, -cstar.R5(:,2) );
%     pj(5) = trapz(timeseq.R6, -cstar.R6(:,2 ) );
%     pjpb(2) = sum(pj);
%     % for beta_3
%     pj = zeros(1,5);
%     pj(1) = trapz(timeseq.R2, - cforward.R2(:,1).*cstar.R2(:,1) );
%     pj(2) = trapz(timeseq.R3, - cforward.R3(:,1).*cstar.R3(:,1) );
%     pj(3) = trapz(timeseq.R4, - cforward.R4(:,1).*cstar.R4(:,1) );
%     pj(4) = trapz(timeseq.R5, - cforward.R5(:,1).*cstar.R5(:,1) );
%     pj(5) = trapz(timeseq.R6, - cforward.R6(:,1).*cstar.R6(:,1) );
%     pjpb(3) = sum(pj);
%     % for beta_4
%     pj = zeros(1,5);
%     pj(1) = trapz(timeseq.R2, - cforward.R2(:,2).*cstar.R2(:,1) );
%     pj(2) = trapz(timeseq.R3, - cforward.R3(:,2).*cstar.R3(:,1) );
%     pj(3) = trapz(timeseq.R4, - cforward.R4(:,2).*cstar.R4(:,1) );
%     pj(4) = trapz(timeseq.R5, - cforward.R5(:,2).*cstar.R5(:,1) );
%     pj(5) = trapz(timeseq.R6, - cforward.R6(:,2).*cstar.R6(:,1) );
%     pjpb(4) = sum(pj);
%     % for beta_5
%     pj = zeros(1,5);
%     pj(1) = trapz(timeseq.R2, - cforward.R2(:,1).*cstar.R2(:,2) );
%     pj(2) = trapz(timeseq.R3, - cforward.R3(:,1).*cstar.R3(:,2) );
%     pj(3) = trapz(timeseq.R4, - cforward.R4(:,1).*cstar.R4(:,2) );
%     pj(4) = trapz(timeseq.R5, - cforward.R5(:,1).*cstar.R5(:,2) );
%     pj(5) = trapz(timeseq.R6, - cforward.R6(:,1).*cstar.R6(:,2) );
%     pjpb(5) = sum(pj);
%     % for beta_6
%     pj = zeros(1,5);
%     pj(1) = trapz(timeseq.R2, - cforward.R2(:,2).*cstar.R2(:,2) );
%     pj(2) = trapz(timeseq.R3, - cforward.R3(:,2).*cstar.R3(:,2) );
%     pj(3) = trapz(timeseq.R4, - cforward.R4(:,2).*cstar.R4(:,2) );
%     pj(4) = trapz(timeseq.R5, - cforward.R5(:,2).*cstar.R5(:,2) );
%     pj(5) = trapz(timeseq.R6, - cforward.R6(:,2).*cstar.R6(:,2) );
%     pjpb(6) = sum(pj);
    
    pjpb_history(n,:) = pjpb;
    
    % Perform line minimazation to determine step length and compute the 
    %       updates for constitutive relations. Two step length values are 
    %       required (for omega_1 and omega_2), refer to manuscript for 
    %       coupled or decoupled formulation! The descent direction is 
    %       determined from the H1 grad of J wrt omega_1 and omega_2!
    %       First run mnbrak.m for getting a triplet of candidate step lengths, then
    %       run brent.m to get the optimized step length!
    
    % Compute the descent direction (DD computed from grad J, a vector evaluated at cgrid points, 
    %       negative of gradient is the descent direction
   
    % for beta_1
    which = 1;  % an integer (1 or 2 or 3 or 4 or 5 or 6) indicating which function or parameter is updated, beta_1 to beta_6
    G2_OLD_1 = G2_1;
    G2_1 = pjpb(1)*pjpb(1);
    if ( n == 1 )      
        DDl  = -pjpb(1);
    else
        switch GRADIENT
            case 1              % Simple gradient;
                MOMENTUM = 0.0; 
            case 2              % Fletcher-Reeves Conjugate Gradient;
                MOMENTUM = G2_1 / G2_OLD_1;
            case 3              % Polak-Ribiere Conjugate Gradient;
                G_DOT_G_OLD_1 = pjpb(1)*pjpb_history(n-1,1);
                MOMENTUM = (G2_1 - G_DOT_G_OLD_1) / G2_OLD_1;
        end
        DD0 = DDl;
        DDl  = -pjpb(1) + MOMENTUM * DD0;
    end
    tau1 = 0.0;
    if ( n == 1 )
        tau2 = ALPHA0;
    else
        tau2 = tau_beta_1;
    end
    Jtau1 = Jeval_ensumble_OCV(tau1, beta, DDl, which);
    Jtau2 = Jeval_ensumble_OCV(tau2, beta, DDl, which);
    [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak_ensumble_OCV(tau1, tau2, Jtau1, Jtau2, beta, DDl, which, MAXITER);
    [tau_beta_1, Jtau_beta_1] = Brent_ensumble_OCV(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, beta, DDl, which);
    Jtau_history(n,1) = Jtau_beta_1;
    % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
    beta_10 = beta(1);
    beta(1) =  beta_10 + tau_beta_1 * DDl;
    fprintf(fout,'   line search for beta_1 performed ... \n');   
    
    % for beta_2
    which = 2;  % an integer (1 or 2 or 3 ) indicating which function or parameter is updated, beta_1 to beta_3
    G2_OLD_2 = G2_2;
    G2_2 = pjpb(2)*pjpb(2);
    if ( n == 1 )      
        DDl  = -pjpb(2);
    else
        switch GRADIENT
            case 1              % Simple gradient;
                MOMENTUM = 0.0; 
            case 2              % Fletcher-Reeves Conjugate Gradient;
                MOMENTUM = G2_2 / G2_OLD_2;
            case 3              % Polak-Ribiere Conjugate Gradient;
                G_DOT_G_OLD_2 = pjpb(2)*pjpb_history(n-1,2);
                MOMENTUM = (G2_2 - G_DOT_G_OLD_2) / G2_OLD_2;
        end
        DD0 = DDl;
        DDl  = -pjpb(2) + MOMENTUM * DD0;
    end
    tau1 = 0.0;
    if ( n == 1 )
        tau2 = ALPHA0;
    else
        tau2 = tau_beta_2;
    end
    Jtau1 = Jeval_ensumble_OCV(tau1, beta, DDl, which);
    Jtau2 = Jeval_ensumble_OCV(tau2, beta, DDl, which);
    [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak_ensumble_OCV(tau1, tau2, Jtau1, Jtau2, beta, DDl, which, MAXITER);
    [tau_beta_2, Jtau_beta_2] = Brent_ensumble_OCV(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, beta, DDl, which);
    Jtau_history(n,2) = Jtau_beta_2;
    % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
    beta_20 = beta(2);
    beta(2) =  beta_20 + tau_beta_2 * DDl;
    fprintf(fout,'   line search for beta_3 performed ... \n');
    
    
    % for beta_3
    which = 3;  % an integer (1 or 2 or 3) indicating which function or parameter is updated, beta_1 to beta_3
    G2_OLD_3 = G2_3;
    G2_3 = pjpb(3)*pjpb(3);
    if ( n == 1 )      
        DDl  = -pjpb(3);
    else
        switch GRADIENT
            case 1              % Simple gradient;
                MOMENTUM = 0.0; 
            case 2              % Fletcher-Reeves Conjugate Gradient;
                MOMENTUM = G2_3 / G2_OLD_3;
            case 3              % Polak-Ribiere Conjugate Gradient;
                G_DOT_G_OLD_3 = pjpb(3)*pjpb_history(n-1,3);
                MOMENTUM = (G2_3 - G_DOT_G_OLD_3) / G2_OLD_3;
        end
        DD0 = DDl;
        DDl  = -pjpb(3) + MOMENTUM * DD0;
    end
    tau1 = 0.0;
    if ( n == 1 )
        tau2 = ALPHA0;
    else
        tau2 = tau_beta_3;
    end
    Jtau1 = Jeval_ensumble_OCV(tau1, beta, DDl, which);
    Jtau2 = Jeval_ensumble_OCV(tau2, beta, DDl, which);
    [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak_ensumble_OCV(tau1, tau2, Jtau1, Jtau2, beta, DDl, which, MAXITER);
    [tau_beta_3, Jtau_beta_3] = Brent_ensumble_OCV(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, beta, DDl, which);
    Jtau_history(n,3) = Jtau_beta_3;
    % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
    beta_30 = beta(3);
    beta(3) =  beta_30 + tau_beta_3 * DDl;
    fprintf(fout,'   line search for beta_3 performed ... \n');
    
    % for beta_4
    which = 4;  % an integer (1 or 2 or 3) indicating which function or parameter is updated, beta_1 to beta_4
    G2_OLD_4 = G2_4;
    G2_4 = pjpb(4)*pjpb(4);
    if ( n == 1 )      
        DDl  = -pjpb(4);
    else
        switch GRADIENT
            case 1              % Simple gradient;
                MOMENTUM = 0.0; 
            case 2              % Fletcher-Reeves Conjugate Gradient;
                MOMENTUM = G2_4 / G2_OLD_4;
            case 3              % Polak-Ribiere Conjugate Gradient;
                G_DOT_G_OLD_4 = pjpb(4)*pjpb_history(n-1,4);
                MOMENTUM = (G2_4 - G_DOT_G_OLD_4) / G2_OLD_4;
        end
        DD0 = DDl;
        DDl  = -pjpb(4) + MOMENTUM * DD0;
    end
    tau1 = 0.0;
    if ( n == 1 )
        tau2 = ALPHA0;
    else
        tau2 = tau_beta_4;
    end
    Jtau1 = Jeval_ensumble_OCV(tau1, beta, DDl, which);
    Jtau2 = Jeval_ensumble_OCV(tau2, beta, DDl, which);
    [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak_ensumble_OCV(tau1, tau2, Jtau1, Jtau2, beta, DDl, which, MAXITER);
    Jtau_beta_40 = Jtau_beta_4;
    [tau_beta_4, Jtau_beta_4] = Brent_ensumble_OCV(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, beta, DDl, which);
    Jtau_history(n,4) = Jtau_beta_4;
    % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
    beta_40 = beta(4);
    beta(4) =  beta_40 + tau_beta_4 * DDl;
    fprintf(fout,'   line search for beta_4 performed ... \n');
    
%     % for beta_5
%     which = 5;  % an integer (1 or 2 or 3) indicating which function or parameter is updated, beta_1 to beta_3
%     G2_OLD_5 = G2_5;
%     G2_5 = pjpb(5)*pjpb(5);
%     if ( n == 1 )      
%         DDl  = -pjpb(5);
%     else
%         switch GRADIENT
%             case 1              % Simple gradient;
%                 MOMENTUM = 0.0; 
%             case 2              % Fletcher-Reeves Conjugate Gradient;
%                 MOMENTUM = G2_5 / G2_OLD_5;
%             case 3              % Polak-Ribiere Conjugate Gradient;
%                 G_DOT_G_OLD_5 = pjpb(5)*pjpb_history(n-1,5);
%                 MOMENTUM = (G2_5 - G_DOT_G_OLD_5) / G2_OLD_5;
%         end
%         DD0 = DDl;
%         DDl  = -pjpb(4) + MOMENTUM * DD0;
%     end
%     tau1 = 0.0;
%     if ( n == 1 )
%         tau2 = ALPHA0;
%     else
%         tau2 = tau_beta_5;
%     end
%     Jtau1 = Jeval_ensumble_OCV(tau1, beta, DDl, which);
%     Jtau2 = Jeval_ensumble_OCV(tau2, beta, DDl, which);
%     [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak_ensumble_OCV(tau1, tau2, Jtau1, Jtau2, beta, DDl, which, MAXITER);
%     [tau_beta_5, Jtau_beta_5] = Brent_ensumble_OCV(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, beta, DDl, which);
%     Jtau_history(n,5) = Jtau_beta_5;
%     % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
%     beta_50 = beta(5);
%     beta(5) =  beta_50 + tau_beta_5 * DDl;
%     fprintf(fout,'   line search for beta_5 performed ... \n');
%     
%     % for beta_6
%     which = 6;  % an integer (1 or 2 or 3) indicating which function or parameter is updated, beta_1 to beta_3
%     G2_OLD_6 = G2_6;
%     G2_6 = pjpb(6)*pjpb(6);
%     if ( n == 1 )      
%         DDl  = -pjpb(6);
%     else
%         switch GRADIENT
%             case 1              % Simple gradient;
%                 MOMENTUM = 0.0; 
%             case 2              % Fletcher-Reeves Conjugate Gradient;
%                 MOMENTUM = G2_6 / G2_OLD_6;
%             case 3              % Polak-Ribiere Conjugate Gradient;
%                 G_DOT_G_OLD_6 = pjpb(6)*pjpb_history(n-1,6);
%                 MOMENTUM = (G2_6 - G_DOT_G_OLD_6) / G2_OLD_6;
%         end
%         DD0 = DDl;
%         DDl  = -pjpb(6) + MOMENTUM * DD0;
%     end
%     tau1 = 0.0;
%     if ( n == 1 )
%         tau2 = ALPHA0;
%     else
%         tau2 = tau_beta_6;
%     end
%     Jtau1 = Jeval_ensumble_OCV(tau1, beta, DDl, which);
%     Jtau2 = Jeval_ensumble_OCV(tau2, beta, DDl, which);
%     [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak_ensumble_OCV(tau1, tau2, Jtau1, Jtau2, beta, DDl, which, MAXITER);
%     Jtau_beta_60 = Jtau_beta_6;
%     [tau_beta_6, Jtau_beta_6] = Brent_ensumble_OCV(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, beta, DDl, which);
%     Jtau_history(n,6) = Jtau_beta_6;
%     % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
%     beta_60 = beta(6);
%     beta(6) =  beta_60 + tau_beta_6 * DDl;
%     fprintf(fout,'   line search for beta_6 performed ... \n');
    
    % save history
    beta_history(n+1,:) = beta;
    
    % print in file
    fprintf(fout,'   Constitutive relations are updated ... \n');
    fprintf(fout,'   Stopping criteria at end of this iteration:\n');
    fprintf(fout,'        iterations: %d <= %d \n', n, MAX_ITER);
    fprintf(fout,'        tolerance for beta_6: %12.8f > %12.8f \n', abs(Jtau_beta_40 - Jtau_beta_4) / abs(Jtau_beta_4), Tol_J);
    fprintf(fout,'   Proceeding to next iteration ... \n\n\n');
    
end


% plot cost function history
figure('visible','on')
semilogy(1:size(Jtau_history,1),Jtau_history(:,3)./Jtau_history(1,3),'Marker', 'o', 'MarkerSize',8)
xlabel('Iteration No.','interpreter','latex','FontSize',15)
ylabel('$\mathcal{J}(\bf\beta^{(n+1)})/ \mathcal{J}(\bf\beta^{(0)})$','interpreter','latex','FontSize',15)
loc = append(pwd,foldername,'/_J_history.png');
saveas(gcf, loc);



% printing results to output file info.text
fprintf(fout,'     Optima beta value = %18.8f \n',beta);

fprintf(fout,'   Histogram of lambda values:\n');
fprintf(fout,'%18s \r\n','lambda');
fprintf(fout,'%18.8f \r\n',beta_history);
fprintf(fout,'\n\n\n');


fprintf(fout,'   Histogram of cost function values:\n');
fprintf(fout,'%18s \r\n','Jtau');
fprintf(fout,'%18.8f \r\n',Jtau_history);
fprintf(fout,'\n\n\n');


fclose(fout);








function dydx = bvpfcn(x,y,cgrid,L2gradJ,l) 
% equation to solve BVP problem
dydx = zeros(2,1);
grad = interp1(cgrid, L2gradJ, x,'spline'); % interpolate L2 gradient for a particular concentration
dydx = [y(2)
       y(1)/(l^2) - grad/(l^2)];
end
%--------------------------------
function res = bcfcn_Dirichlet(ya,yb)
% Dirichlet boundary conditions - H1 gradient is zero at boundaries.
res = [ya(1)
    yb(1)];
end
%--------------------------------
function res = bcfcn_Newmann(ya,yb)
% Newmann boundary conditions - the gradient of H1 gradient is zero at boundaries. 
res = [ya(2)
    yb(2)];
end
end
