function [omega_1, omega_2, lambda, c1forward_history, c2forward_history, L2gradJ1_history, L2gradJ2_history, H1gradJ1_history, H1gradJ2_history, omega_1_history, omega_2_history, lambda_history] = iterative(omega_1, omega_2, lambda, plt, foldername)
%%% -----------------------------------------------------------------
% the iterative scheme of minimization
% Input:     experimental data in time
% Output:    computed omega(c1,c2)

global c1grid c2grid Ngrid timeseq Nt Feval LineMin1 LineMin2 c1alpha c1beta c2alpha c2beta I rho_1 rho_2

% write info in a file
fout = fopen('info.txt','w');
fprintf(fout,'Parameters of the simulation are as follows: \n');


%%% OPTIMIZATION PARAMETERS
MAX_ITER       = 50;        % Maximum number of iterations of algorithm;
MAXITER        = 100;       % Maximum number of bracketing iterations;
ALPHA0         = 1;         % Initial interval for bracketing;
Tol_J          = 1.0e-9;    % Tolerance of iterative algorithm.
MAX_ITER_BRENT = 50;        % Max iterations for Brent routine.
Tol_Brent      = 1.0e-4;    % Tolerance in Brent routine.
MOMENTUM       = 0.0;       % Momentum for conjugate gradient algorithm
G2_1           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for C_1
G2_2           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for C_2
G2_3           = 0.0;       % initialization for (Gradient dot gradient) in Fletcher Reeves CG method for lambda
GRADIENT       = 3;         % type of conjugate gradient (1 for simple gradient, 2 for Fletcher-Reeves CG, and 3 for Polak-Ribiere CG)


fprintf(fout,'     Maximum number of optimal reconstruction iterations = %6.2f \n',MAX_ITER);
fprintf(fout,'     Maximum number of bracketing iterations = %6.2f \n',MAXITER);
fprintf(fout,'     Maximum number of brent routine iterations = %6.2f \n',MAX_ITER_BRENT);
fprintf(fout,'     Initial interval for bracketing = %6.2f \n',ALPHA0);
fprintf(fout,'     Tolerance of iterative algorithm = %6.2f \n',Tol_J);
fprintf(fout,'     Tolerance of brent algorithm = %6.2f \n',Tol_Brent);
fprintf(fout,'     Type pf conjugate gradient used = %6.2f \n',GRADIENT);


%%% Parameters for solving BVP system for sobolev gradients:
m = 2;                      % Boundary condition type, m = 1 for Dirichlet and m = 2 for Newmann.
l = 0.1;                    % The smoothing coefficient of H1 gradients, higher l means higher smoothing.

if m==1
    fprintf(fout,'     m = %d, Dirichlet Boundary Condition for Computing H1 gradients. \n', m);
else
    fprintf(fout,'     m = %d, Newmann Boundary Condition for Computing H1 gradients. \n', m);
end
fprintf(fout,'     l = %6.2f, the penalty parameter of BVP problem, determining the smoothness of H1 gradients. \n\n\n',l);


% set iterator to zero
n = 0;

% Store intermediate solutions in a table, rows correspond to iteration number
c1forward_history = zeros(1,Nt);      % forward solutions for C_1
c2forward_history = zeros(1,Nt);      % forward solutions for C_2
c1adjoint_history = zeros(1,Nt);      % adjoint solutions for C_1
c2adjoint_history = zeros(1,Nt);      % adjoint solutions for C_2
L2gradJ1_history = zeros(1,Ngrid);    % L2 gradient of J with respect to omega_1
L2gradJ2_history = zeros(1,Ngrid);    % L2 gradient of J with respect to omega_2
H1gradJ1_history = zeros(1,Ngrid);    % H1 gradient of J with respect to omega_1
H1gradJ2_history = zeros(1,Ngrid);    % H1 gradient of J with respect to omega_2
pjpl_history = zeros(1,1);            % partial derivaive of J with respect to lambda

Jtau_omega_1  = 1.0;
Jtau_omega_10 = 0.0;
Jtau_omega_2  = 1.0;
Jtau_omega_20 = 0.0;
Jtau_lambda  = 1.0;
Jtau_lambda0 = 0.0;


Jtau_history_1 = zeros(1,2);          % History of cost function values for omega_1
Jtau_history_2 = zeros(1,2);          % History of cost function values for omega_2
Jtau_history_3 = zeros(1,2);          % History of cost function values for lambda

omega_1_history = omega_1;
omega_2_history = omega_2;
lambda_history = lambda;

% iterative scheme
while (n <= MAX_ITER) && (abs(Jtau_omega_20 - Jtau_omega_2) / abs(Jtau_omega_2) > Tol_J)
    
    % Update iterator
    n = n+1;
    fprintf(fout, '   n = %d \n',n);
    fprintf('   n = %d \n',n);
    
    % Solve forawrd problem
    [cforward] = forward(omega_1, omega_2, lambda);
    c1forward_history(n,:) = transpose(cforward(:,1));
    c2forward_history(n,:) = transpose(cforward(:,2));
    if plt
        figure('visible','off')
        hold on
        plot(timeseq, cforward(:,1),'DisplayName','$C_1(t)$','LineWidth',2)
        plot(timeseq, cforward(:,2),'DisplayName','$C_2(t)$','LineWidth',2)
        legend('interpreter','latex','FontSize',15)
        xlabel('$t$','interpreter','latex','FontSize',15)
        ylabel('Concentration','interpreter','latex','FontSize',15)
        hold off
        loc = append(pwd,foldername,'/',num2str(n),'_forward.png');
        saveas(gcf, loc);
    end
    fprintf(fout,'   Forward problem solved ... \n');
    
    
    % Solve adjoint problem
    [cstar] = adjoint(omega_1, omega_2, cforward, lambda);
    c1adjoint_history(n,:) = transpose(cstar(:,1));
    c2adjoint_history(n,:) = transpose(cstar(:,2));
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
    
    % Compute L2 gradient of J wrt constitutive relations
    L2gradJ1 = zeros(1,Ngrid);
    L2gradJ2 = zeros(1,Ngrid);
    gt1 = zeros(size(timeseq));
    gt2 = zeros(size(timeseq));
    for i = 1:size(timeseq,1)
        om1 = interp1(c1grid, omega_1, cforward(i,1),'spline');
        om2 = interp1(c2grid, omega_2, cforward(i,2),'spline');
        gt1(i) = + cstar(i,1)/om1 - 3*cstar(i,2)/om1;
        gt2(i) = + (cstar(i,1)*om1)/(3*(1-om1*om2)) - (cstar(i,2)*om1)/(1-om1*om2);
    end
    for i=1:Ngrid
        L2gradJ1(i) = non_unique_interpolation(cforward(:,1), gt1, c1grid(i));
        L2gradJ2(i) = non_unique_interpolation(cforward(:,2), gt2, c2grid(i));
    end
    L2gradJ1_history(n,:) = L2gradJ1;
    L2gradJ2_history(n,:) = L2gradJ2;
    if plt
        figure('visible','off')
        plot(c1grid, L2gradJ1,'LineWidth',2)
        xlabel('$C_1$','interpreter','latex','FontSize',15)
        ylabel('$\nabla_{\omega_1}^{L^2} \mathcal{J}$','interpreter','latex','FontSize',15)
        loc = append(pwd,foldername,'/',num2str(n),'_L2gradJ1.png');
        saveas(gcf, loc);
        figure('visible','off')
        plot(c2grid, L2gradJ2,'LineWidth',2)
        xlabel('$C_2$','interpreter','latex','FontSize',15)
        ylabel('$\nabla_{\omega_2}^{L^2} \mathcal{J}$','interpreter','latex','FontSize',15)
        loc = append(pwd,foldername,'/',num2str(n),'_L2gradJ2.png');
        saveas(gcf, loc);
    end
    fprintf(fout,'   L2 gradients computed ... \n');
    
    % computing partial derivative of J w.r.t lambda
    integrand = zeros(size(timeseq));
    for i = 1:size(timeseq,1)
        om1 = interp1(c1grid, omega_1, cforward(i,1),'spline');
        om2 = interp1(c2grid, omega_2, cforward(i,2),'spline');
        integrand(i) = rho_1(I(i))*om1*om2*cstar(i,1) + rho_2(I(i))*(1-om1*om2)*cstar(i,2);
    end
    pjpl = trapz(timeseq, integrand);
    pjpl_history(n,1) = pjpl;
    
    % Solve BVP to obtain H1 gradient of J
    %opts = bvpset('RelTol',1e-5,'Stats','off');
    %guess = [0; 0];
    %xmesh = linspace(0,1,100);
    %solinit = bvpinit(xmesh, guess);
    %switch m
    %    case 1
    %        sol1 = bvp5c(@(x,y) bvpfcn(x,y,c1grid,L2gradJ1,l), @bcfcn_Dirichlet, solinit, opts);
    %        sol2 = bvp5c(@(x,y) bvpfcn(x,y,c2grid,L2gradJ2,l), @bcfcn_Dirichlet, solinit, opts);
    %    case 2
    %        sol1 = bvp5c(@(x,y) bvpfcn(x,y,c1grid,L2gradJ1,l), @bcfcn_Newmann, solinit, opts);
    %        sol2 = bvp5c(@(x,y) bvpfcn(x,y,c2grid,L2gradJ2,l), @bcfcn_Newmann, solinit, opts);
    %    otherwise
    %        disp('     Select an option for boundary conditions!')
    %end
    %H1gradJ1 = interp1(sol1.x, sol1.y(1,:), c1grid,'spline');
    %H1gradJ2 = interp1(sol2.x, sol2.y(1,:), c2grid,'spline');
    % My BVP solver, using Newmann BCs, 
    meshsize = 500;
    H1gradJ1 = zeros(Ngrid,1);
    H1gradJ2 = zeros(Ngrid,1);
    H1gradJ1 = BVP_solve(meshsize, c1alpha, c1beta, l, L2gradJ1, c1grid, 1);
    H1gradJ2 = BVP_solve(meshsize, c2alpha, c2beta, l, L2gradJ2, c2grid, 2);
    H1gradJ1_history(n,:) = H1gradJ1;
    H1gradJ2_history(n,:) = H1gradJ2;
    if plt
        figure('visible','off')
        plot(c1grid, H1gradJ1,'LineWidth',2)
        xlabel('$C_1$','interpreter','latex','FontSize',15)
        ylabel('$\nabla_{\omega_1}^{H^1} \mathcal{J}$','interpreter','latex','FontSize',15)
        loc = append(pwd,foldername,'/',num2str(n),'_H1gradJ1.png');
        saveas(gcf, loc);
        figure('visible','off')
        plot(c2grid, H1gradJ2,'LineWidth',2)
        xlabel('$C_2$','interpreter','latex','FontSize',15)
        ylabel('$\nabla_{\omega_2}^{H^1} \mathcal{J}$','interpreter','latex','FontSize',15)
        loc = append(pwd,foldername,'/',num2str(n),'_H1gradJ2.png');
        saveas(gcf, loc);
    end
    fprintf(fout,'   H1 gradients computed ... \n');
    
    % Perform line minimazation to determine step length and compute the 
    %       updates for constitutive relations. Two step length values are 
    %       required (for omega_1 and omega_2), refer to manuscript for 
    %       coupled or decoupled formulation! The descent direction is 
    %       determined from the H1 grad of J wrt omega_1 and omega_2!
    %       First run mnbrak.m for getting a triplet of candidate step lengths, then
    %       run brent.m to get the optimized step length!
    
    % Compute the descent direction (DD computed from grad J, a vector evaluated at cgrid points, 
    %       negative of gradient is the descent direction
    
    % for lambda
    which = 3;  % an integer (1 or 2 or 3) indicating which function or parameter is updated, omega_1 or omega_2 or lambda
    G2_OLD_3 = G2_3;
    G2_3 = pjpl*pjpl;
    if ( n == 1 )      
        DDl  = -pjpl;
    else
        switch GRADIENT
            case 1              % Simple gradient;
                MOMENTUM = 0.0; 
            case 2              % Fletcher-Reeves Conjugate Gradient;
                MOMENTUM = G2_3 / G2_OLD_3;
            case 3              % Polak-Ribiere Conjugate Gradient;
                G_DOT_G_OLD_3 = pjpl*pjpl_history(n-1,1);
                MOMENTUM = (G2_3 - G_DOT_G_OLD_3) / G2_OLD_3;
        end
        DD0 = DDl;
        DDl  = -pjpl + MOMENTUM * DD0;
    end
    tau1 = 0.0;
    if ( n == 1 )
        tau2 = ALPHA0;
    else
        tau2 = tau_lambda;
    end
    lambda
    Jtau1 = Jeval(tau1, omega_1, omega_2, lambda, DDl, which);
    Jtau2 = Jeval(tau2, omega_1, omega_2, lambda, DDl, which);
    % mnbrak for omega_1
    [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak(tau1, tau2, Jtau1, Jtau2, omega_1, omega_2, lambda, DDl, which, MAXITER);
    % brent for omega_1
    Jtau_lambda0 = Jtau_lambda;
    [tau_lambda, Jtau_lambda] = Brent(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, omega_1, omega_2, lambda, DDl, which);
    Jtau_history_3(1,n) = Jtau_lambda;
    % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
    lambda_0 = lambda;
    lambda =  lambda_0 + tau_lambda * DDl;
    fprintf(fout,'   line search for lambda performed ... \n');   
    
    
    
    
    % for omega_1
    which = 1;     
    G2_OLD_1 = G2_1;
    G2_1 = prodH1(H1gradJ1, H1gradJ1, which, l);
    if ( n == 1 )      
        DD1 = -H1gradJ1;
    else
        switch GRADIENT
            case 1              % Simple gradient;
                MOMENTUM = 0.0; 
            case 2              % Fletcher-Reeves Conjugate Gradient;
                MOMENTUM = G2_1 / G2_OLD_1;
            case 3              % Polak-Ribiere Conjugate Gradient;
                G_DOT_G_OLD_1 = prodH1(H1gradJ1, H1gradJ1_history(n-1,:), which, l);
                MOMENTUM = (G2_1 - G_DOT_G_OLD_1) / G2_OLD_1;
        end
        DD0 = DD1;
        DD1  = -H1gradJ1 + MOMENTUM * DD0;
    end 
    fprintf(fout,'   Descent direction for omega_1 computed ... \n');
    % Determine intial guesses for line minimization
    tau1 = 0.0;
    if ( n == 1 )
        tau2 = ALPHA0;
    else
        tau2 = tau_omega_1;
    end
    % Evaluate cost function values for these guesses (compute the updated constitutive 
    %       relation and evaulate cost function by solving the forward problem)
    Jtau1 = Jeval(tau1, omega_1, omega_2, lambda, DD1, which);
    Jtau2 = Jeval(tau2, omega_1, omega_2, lambda, DD1, which);
    % mnbrak for omega_1
    [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak(tau1, tau2, Jtau1, Jtau2, omega_1, omega_2, lambda, DD1, which, MAXITER);
    % brent for omega_1
    Jtau_omega_10 = Jtau_omega_1;
    [tau_omega_1, Jtau_omega_1] = Brent(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, omega_1, omega_2, lambda, DD1, which);
    Jtau_history_1(1,n) = Jtau_omega_1;
    % update omega_1 based on the descent direction and step length (use conjugate gradient for this update)
    omega_1_0 = omega_1;
    omega_1 =  omega_1_0 + tau_omega_1 * DD1;
    fprintf(fout,'   line search for omega_1 performed ... \n');
    
    
    
    %for omega_2:
    which = 2;
    G2_OLD_2 = G2_2;
    G2_2 = prodH1(H1gradJ2, H1gradJ2, which, l);
    if ( n == 1 )      
        DD2 = -H1gradJ2;
    else
        switch GRADIENT
            case 1              % Simple gradient;
                MOMENTUM = 0.0; 
            case 2              % Fletcher-Reeves Conjugate Gradient;
                MOMENTUM = G2_2 / G2_OLD_2;
            case 3              % Polak-Ribiere Conjugate Gradient;
                G_DOT_G_OLD_2 = prodH1(H1gradJ2, H1gradJ2_history(n-1,:), which, l);
                MOMENTUM = (G2_2 - G_DOT_G_OLD_2) / G2_OLD_2;
        end
        DD0 = DD2;
        DD2  = -H1gradJ2 + MOMENTUM * DD0;
    end
    fprintf(fout,'   Descent direction for omega_2 computed ... \n');
    tau1 = 0.0;
    if ( n == 1 )
        tau2 = ALPHA0;
    else
        tau2 = tau_omega_2;
    end
    Jtau1 = Jeval(tau1, omega_1, omega_2, lambda, DD2, which);
    Jtau2 = Jeval(tau2, omega_1, omega_2, lambda, DD2, which);
    [tau1, Jtau1, tau2, Jtau2, tau3, Jtau3] = mnbrak(tau1, tau2, Jtau1, Jtau2, omega_1, omega_2, lambda, DD2, which, MAXITER);
    Jtau_omega_20 = Jtau_omega_2;
    [tau_omega_2, Jtau_omega_2] = Brent(tau1, tau2, tau3, Jtau1, Jtau2, Jtau3, Tol_Brent, MAX_ITER_BRENT, omega_1, omega_2, lambda, DD2, which);
    Jtau_history_2(1,n) = Jtau_omega_2;
    omega_2_0 = omega_2;
    omega_2 =  omega_2_0 + tau_omega_2 * DD2;
    fprintf(fout,'   Line search for omega_2 performed ... \n');
    
    
    
    
    omega_1_history(n+1,:) = omega_1;
    omega_2_history(n+1,:) = omega_2;
    lambda_history(n+1,1) = lambda;
    
    
    
    if plt
        figure('visible','off')
        plot(c1grid, omega_1,'LineWidth',2)
        xlabel('$C_1$','interpreter','latex','FontSize',15)
        ylabel('$\omega_1(C_1)$','interpreter','latex','FontSize',15)
        loc = append(pwd,foldername,'/',num2str(n),'_omega_1.png');
        saveas(gcf, loc);
        figure('visible','off')
        plot(c2grid, omega_2,'LineWidth',2)
        xlabel('$C_2$','interpreter','latex','FontSize',15)
        ylabel('$\omega_2(C_2)$','interpreter','latex','FontSize',15)
        loc = append(pwd,foldername,'/',num2str(n),'_omega_2.png');
        saveas(gcf, loc);
    end
    
    fprintf(fout,'   Constitutive relations are updated ... \n');
    fprintf(fout,'   Stopping criteria at end of this iteration:\n');
    fprintf(fout,'        iterations: %d <= %d \n', n, MAX_ITER);
    fprintf(fout,'        tolerance for lambda:  %12.8f > %12.8f \n', (Jtau_lambda0 - Jtau_lambda) / abs(Jtau_lambda), Tol_J);
    fprintf(fout,'        tolerance for omega_1: %12.8f > %12.8f \n', (Jtau_omega_10 - Jtau_omega_1) / abs(Jtau_omega_1), Tol_J);
    fprintf(fout,'        tolerance for omega_2: %12.8f > %12.8f \n', (Jtau_omega_20 - Jtau_omega_2) / abs(Jtau_omega_2), Tol_J);
    fprintf(fout,'   Proceeding to next iteration ... \n\n\n');
    
end


figure('visible','on')
semilogy(Jtau_history_3./Jtau_history_3(1),'Marker', 'o', 'MarkerSize',8)
ylim([0.01,1])
xlabel('Iteration No.','interpreter','latex','FontSize',15)
ylabel('$\mathcal{J}(\lambda^{(n+1)},\omega_1^{(n)},\omega_2^{(n)})/ \mathcal{J}(\lambda^{(0)},\omega_1^{(0)},\omega_2^{(0)})$','interpreter','latex','FontSize',15)
loc = append(pwd,foldername,'/_J3_history.png');
saveas(gcf, loc);
figure('visible','on')
semilogy(Jtau_history_1./Jtau_history_1(1),'Marker', 'o','MarkerSize',8)
ylim([0.01,1])
xlabel('Iteration No.','interpreter','latex','FontSize',15)
ylabel('$\mathcal{J}(\lambda^{(n+1)},\omega_1^{(n+1)},\omega_2^{(n)})/ \mathcal{J}(\lambda^{(0)},\omega_1^{(0)},\omega_2^{(0)})$','interpreter','latex','FontSize',15)
loc = append(pwd,foldername,'/_J1_history.png');
saveas(gcf, loc);
figure('visible','on')
semilogy(Jtau_history_2./Jtau_history_2(1),'Marker', 'o', 'MarkerSize',8)
ylim([0.01,1])
xlabel('Iteration No.','interpreter','latex','FontSize',15)
ylabel('$\mathcal{J}(\lambda^{(n+1)},\omega_1^{(n+1)},\omega_2^{(n+1)})/ \mathcal{J}(\lambda^{(0)},\omega_1^{(0)},\omega_2^{(0)})$','interpreter','latex','FontSize',15)
loc = append(pwd,foldername,'/_J2_history.png');
saveas(gcf, loc);



% printing results to output file info.text
fprintf(fout,'     Optima lambda value = %18.8f \n',lambda);

fprintf(fout,'   Histogram of lambda values:\n');
fprintf(fout,'%18s \r\n','lambda');
fprintf(fout,'%18.8f \r\n',lambda_history);
fprintf(fout,'\n\n\n');


fprintf(fout,'   Histogram of cost function values with respect to omega_1:\n');
fprintf(fout,'%18s \r\n','Jtau');
fprintf(fout,'%18.8f \r\n',Jtau_history_1);
fprintf(fout,'\n\n\n');

fprintf(fout,'   Histogram of cost function values with respect to omega_2:\n');
fprintf(fout,'%18s \r\n','Jtau');
fprintf(fout,'%18.8f \r\n',Jtau_history_2);
fprintf(fout,'\n\n\n');

fprintf(fout,'   Histogram of cost function values with respect to lambda:\n');
fprintf(fout,'%18s \r\n','Jtau');
fprintf(fout,'%18.8f \r\n',Jtau_history_3);
fprintf(fout,'\n\n\n');

fprintf(fout,'   Total number of function evaluations = %6.2f \n\n\n',Feval);
fprintf(fout,'   Line minimization results for omega_1:\n');
fprintf(fout,'%12s %18s \r\n','tau','Ftau');
fprintf(fout,'%12.8f %18.8f \r\n',LineMin1);
fprintf(fout,'\n\n\n');

fprintf(fout,'   Line minimization results for omega_2:\n');
fprintf(fout,'%12s %18s \r\n','tau','Ftau');
fprintf(fout,'%12.8f %18.8f \r\n',LineMin2);
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
