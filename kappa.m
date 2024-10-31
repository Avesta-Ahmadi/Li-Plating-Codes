function [epsilon, kappa_1_p, kappa_2_p, kappa_3_p, kappa_1_pp, kappa_2_pp, kappa_3_pp, kappa_1_ppp, kappa_2_ppp, kappa_3_ppp] = kappa(omega_1, omega_2, lambda, plt)
%%% -----------------------------------------------------------------
% Computes cost function J1 for c1, based on separation of variables
% formulation
% Input:     omega_1 and omega_2 vectors parametrized by concentration
% Output:    results of the kappa test as a plot
% we assume a omega_prime function as a perturbation to the original system. 
% The forward problem is solved based on the omega function
% and its perturbation and the directional differential based on finite difference
% scheme is computed. Also, the H1 gradients of cost function wrt to
% constitiutve relations are computed and the directional differential is
% computed using Riesz formulation.prodH1
% The kappa will be performed only on the first iteration of the iterative
% algorithm.



global c1grid c2grid c1tilde c2tilde Ngrid timeseq c1alpha c1beta c2alpha c2beta rho_1 rho_2 I

% Parameters for solving BVP system:
% m determines boundary condition, m = 1 for Dirichlet and m = 2 for Newmann
% l determines the smoothing coefficient of H1 gradient, higher l means higher smoothing 
m = 2;
l = 0.1;

% assume a perturbation function. Three different functions are proposed for this purpose:
omega_p   = @(c) (c-0.1).^3;
omega_pp  = @(c) 0.1*c.^2;
omega_ppp = @(c) (1/(2*pi)^0.5)*exp(-0.5*(c).^2);
lambda_prime = 1;

% compute directional derivative based on the Riesz form and L2 gradient
[cforward] = forward(omega_1, omega_2, lambda);
if plt
    figure()
    hold on
    plot(timeseq, cforward(:,1),'DisplayName','$C_1(t)$','LineWidth',2)
    plot(timeseq, cforward(:,2),'DisplayName','$C_2(t)$','LineWidth',2)
    legend('interpreter','latex','FontSize',15)
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('Concentration','interpreter','latex','FontSize',15)
    hold off
end

[cstar] = adjoint(omega_1, omega_2, cforward, lambda);
if plt
    figure()
    hold on
    plot(timeseq, cstar(:,1),'DisplayName','$C^\ast_1(t)$','LineWidth',2)
    plot(timeseq, cstar(:,2),'DisplayName','$C^\ast_2(t)$','LineWidth',2)
    xlabel('$t$','interpreter','latex','FontSize',15)
    ylabel('Concentration','interpreter','latex','FontSize',15)
    legend('interpreter','latex','FontSize',15)
    hold off
end

% computing L2 gradients for omega_1 and omega_2
L2gradJ1 = zeros(Ngrid,1);
L2gradJ2 = zeros(Ngrid,1);
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

% computing partial derivative of J w.r.t lambda
integrand = zeros(size(timeseq));
for i = 1:size(timeseq,1)
    om1 = interp1(c1grid, omega_1, cforward(i,1),'spline');
    om2 = interp1(c2grid, omega_2, cforward(i,2),'spline');
    integrand(i) = rho_1(I(i))*om1*om2*cstar(i,1) + rho_2(I(i))*(1-om1*om2)*cstar(i,2);
end
pjpl = trapz(timeseq, integrand);

if plt
    figure()
    plot(c1grid, L2gradJ1,'LineWidth',2)
    xlabel('$C_1$','interpreter','latex','FontSize',15)
    ylabel('$\nabla_{\omega_1}^{L^2} \mathcal{J}$','interpreter','latex','FontSize',15)
    figure()
    plot(c2grid, L2gradJ2,'LineWidth',2)
    xlabel('$C_2$','interpreter','latex','FontSize',15)
    ylabel('$\nabla_{\omega_2}^{L^2} \mathcal{J}$','interpreter','latex','FontSize',15)
end

% opts = bvpset('RelTol',1e-5,'Stats','off');
% guess = [0; 0];
% xmesh = linspace(0,1,100);
% solinit = bvpinit(xmesh, guess);
% switch m
%     case 1
%         sol1 = bvp5c(@(x,y) bvpfcn(x,y,c1grid,L2gradJ1,l), @bcfcn_Dirichlet, solinit, opts);
%         sol2 = bvp5c(@(x,y) bvpfcn(x,y,c2grid,L2gradJ2,l), @bcfcn_Dirichlet, solinit, opts);
%     case 2
%         sol1 = bvp5c(@(x,y) bvpfcn(x,y,c1grid,L2gradJ1,l), @bcfcn_Newmann, solinit, opts);
%         sol2 = bvp5c(@(x,y) bvpfcn(x,y,c2grid,L2gradJ2,l), @bcfcn_Newmann, solinit, opts);
%     otherwise
%         disp('Select an option for boundary conditions!')
% end
% H1gradJ1 = zeros(Ngrid,1);
% H1gradJ2 = zeros(Ngrid,1);
% H1gradJ1 = interp1(sol1.x, sol1.y(1,:), c1grid,'spline');
% H1gradJ2 = interp1(sol2.x, sol2.y(1,:), c2grid,'spline');

% My BVP solver, using Newmann BCs, 
meshsize = 500;
H1gradJ1 = zeros(Ngrid,1);
H1gradJ2 = zeros(Ngrid,1);
H1gradJ1 = BVP_solve(meshsize, c1alpha, c1beta, l, L2gradJ1, c1grid, 1);
H1gradJ2 = BVP_solve(meshsize, c2alpha, c2beta, l, L2gradJ2, c2grid, 2);

if plt
    figure()
    plot(c1grid, H1gradJ1,'LineWidth',2)
    xlabel('$C_1$','interpreter','latex','FontSize',15)
    ylabel('$\nabla_{\omega_1}^{H^1} \mathcal{J}$','interpreter','latex','FontSize',15)
    figure()
    plot(c2grid, H1gradJ2,'LineWidth',2)
    xlabel('$C_2$','interpreter','latex','FontSize',15)
    ylabel('$\nabla_{\omega_2}^{H^1} \mathcal{J}$','interpreter','latex','FontSize',15)
end





% for loop to iterate for differet epsilon values for a specific omega_prime funciton
epsilon = 10.^(0:-1:-16);
kappa_1_p = zeros(1,max(size(epsilon)));
kappa_2_p = zeros(1,max(size(epsilon)));
den1 = trapz(c1grid, transpose(L2gradJ1).*omega_p(c1grid));
den2 = trapz(c2grid, transpose(L2gradJ2).*omega_p(c2grid));
den3 = pjpl*lambda_prime;
%den1 = trapz(c1grid, H1gradJ1.*omega_pp(c1grid));
%den2 = trapz(c2grid, H1gradJ2.*omega_pp(c2grid));
for j = 1:max(size(epsilon))
    disp(epsilon(j))
    % compute directional derivative based on finite differencing
    [c_perturbed] = forward(omega_1 + epsilon(j).* omega_p(c1grid), omega_2, lambda);
    num1 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    [c_perturbed] = forward(omega_1, omega_2 + epsilon(j).* omega_p(c2grid), lambda);
    num2 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    [c_perturbed] = forward(omega_1, omega_2, lambda + epsilon(j)*lambda_prime);
    num3 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    kappa_1_p(j) = num1/den1;
    kappa_2_p(j) = num2/den2;
    kappa_3_p(j) = num3/den3;
end


% for loop to iterate for differet epsilon values for a specific omega_prime funciton
kappa_1_pp = zeros(1,max(size(epsilon)));
kappa_2_pp = zeros(1,max(size(epsilon)));
den1 = trapz(c1grid, transpose(L2gradJ1).*omega_pp(c1grid));
den2 = trapz(c2grid, transpose(L2gradJ2).*omega_pp(c2grid));
%den1 = trapz(c1grid, H1gradJ1.*omega_pp(c1grid));
%den2 = trapz(c2grid, H1gradJ2.*omega_pp(c2grid));
for j = 1:max(size(epsilon))
    disp(epsilon(j))
    % compute directional derivative based on finite differencing
    [c_perturbed] = forward(omega_1 + epsilon(j).* omega_pp(c1grid), omega_2, lambda);
    num1 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    [c_perturbed] = forward(omega_1, omega_2 + epsilon(j).* omega_pp(c2grid), lambda);
    num2 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    [c_perturbed] = forward(omega_1, omega_2, lambda + epsilon(j)*lambda_prime);
    num3 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    kappa_1_pp(j) = num1/den1;
    kappa_2_pp(j) = num2/den2;
    kappa_3_pp(j) = num3/den3;
end



% for loop to iterate for differet epsilon values for a specific omega_prime funciton
kappa_1_ppp = zeros(1,max(size(epsilon)));
kappa_2_ppp = zeros(1,max(size(epsilon)));
den1 = trapz(c1grid, transpose(L2gradJ1).*omega_ppp(c1grid));
den2 = trapz(c2grid, transpose(L2gradJ2).*omega_ppp(c2grid));
%den1 = trapz(c1grid, H1gradJ1.*omega_pp(c1grid));
%den2 = trapz(c2grid, H1gradJ2.*omega_pp(c2grid));
for j = 1:max(size(epsilon))
    disp(epsilon(j))
    % compute directional derivative based on finite differencing
    [c_perturbed] = forward(omega_1 + epsilon(j).* omega_ppp(c1grid), omega_2, lambda);
    num1 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    [c_perturbed] = forward(omega_1, omega_2 + epsilon(j).* omega_ppp(c2grid), lambda);
    num2 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    [c_perturbed] = forward(omega_1, omega_2, lambda + epsilon(j)*lambda_prime);
    num3 = (J(c_perturbed,c1tilde,c2tilde,timeseq) - J(cforward,c1tilde,c2tilde,timeseq))/epsilon(j);
    
    kappa_1_ppp(j) = num1/den1;
    kappa_2_ppp(j) = num2/den2;
    kappa_3_ppp(j) = num3/den3;
end
















function dydx = bvpfcn(x,y,cgrid,L2gradJ,l) 
% equation to solve BVP problem
dydx = zeros(2,1);
grad = interp1(cgrid, L2gradJ, x,'spline'); % interpolate L2 gradient for a particular concentration
dydx = [y(2)
       y(1)/(l^2) - grad/(l^2)];
end
%--------------------------------
function res = bcfcn_Dirichlet(ya,yb)
% Dirichlet boundary conditions - H1 gradient is zero at boundaries
res = [ya(1)
    yb(1)];
end
%--------------------------------
function res = bcfcn_Newmann(ya,yb)
% Newmann boundary conditions - the gradient of H1 gradient is zero at
% boundaries. 
res = [ya(2)
    yb(2)];
end


end
