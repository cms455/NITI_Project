function [G_opt, mu_opt, mu_opt_div_G, kfit_final, confidence_interval] = fit_data_to_curve_lsq(rho, rho_l, c_l, cp, h, G0, mu0, mu_factor, f,k, VdB, f_reduced_idx, VdB_error)
if nargin < 1
    rho = 1000;    % Density of material (kg/m^3)
end
if nargin < 2
    rho_l = 1000;  % Density of liquid (kg/m^3)
end
if nargin < 3
    c_l = 1480;    % Speed of sound in liquid (m/s)
end
if nargin < 4
    cp = 1540;     % Phase velocity (m/s)
end
if nargin < 5
    h = 800E-6;    % Thickness (m)
end
if nargin < 6
    G0 = 20E3;     % Initial guess for G (Pa)
end
if nargin < 7
    mu_factor = 110; % mu = mu_factor * G;
end
if nargin < 8
    mu0 = G0 * mu_factor;    % Initial guess for mu (Pa)
end

lambda = rho * cp^2 - 2 * mu0; % Lame's constant (Pa)

% Define the objective function
lb = [1e3]; % G must be positive (Weird stuff happens for G < 1e3)
ub = []; % No upper bound in this case
%objfun = @(p) calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, p(1), p(1)*mu_factor, rho, rho_l, c_l, cp);
objfun = @(G,f) calc_niti_amode_loss(G, f, k, VdB, h, f_reduced_idx, mu_factor, rho, rho_l, c_l, cp);

% Perform optimization + Optimization Constraints
opts = optimset('Display', 'iter', ...
    'TolX', 1e-6, ...         % Final convergence tolerance for the step size
    'TolFun', 1e-6);      % Final convergence tolerance for the function value

figure;
hold on;
plot(f(f_reduced_idx), VdB);
title('Converging G Plot')

[popt, ~, residual, ~, ~, ~, J] = lsqcurvefit(objfun, G0, f, VdB, lb, ub, opts);

G_opt = popt(1);
mu_opt = G_opt * mu_factor;
mu_opt_div_G = mu_opt / G_opt;

confidence_interval = nlparci(popt, residual, 'jacobian', J);

% Plot the optimized curve
cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);
kfit = f(:) ./ cfit(:);

% Calculate original guess curve
%cfit_0 = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
%kfit_0 = f(:) ./ cfit_0(:);

%Calculate confidence interval curves
cfit_ci_1 = compute_niti_amode(f, k, h, confidence_interval(1), confidence_interval(1)*mu_factor, lambda, rho, rho_l, c_l);
kfit_ci_1 = f(:) ./ cfit_ci_1(:);
cfit_ci_2 = compute_niti_amode(f, k, h, confidence_interval(2), confidence_interval(2)*mu_factor, lambda, rho, rho_l, c_l);
kfit_ci_2 = f(:) ./ cfit_ci_2(:);
% Create a plot comparing the fitted curve with original guess and the data
% Define the purple color as a variable
purpleColor = [0.5, 0, 0.5]; % Change this array to substitute with another color

figure;
hold on;

% Plot the data points with error bars in dark green
errorbar(f(f_reduced_idx), VdB * 1e-3, VdB_error, 'o', ...
    'MarkerSize', 8, 'MarkerFaceColor', 'g', 'Color', 'g', ...
    'LineStyle', 'none', 'DisplayName', 'Measured Data');

% Plot the fitted curve using the purpleColor variable
plot(f, kfit * 1e-3, '-', 'Color', purpleColor, 'LineWidth', 2, ...
    'DisplayName', 'Fitted Curve');

% Add light shading between the confidence intervals

kfit_ci_1(1) = [];
kfit_ci_2(1) = [];
fill([1:length(kfit_ci_1), fliplr(1:length(kfit_ci_2))], ...
     [kfit_ci_1' * 1e-3, fliplr(kfit_ci_2') * 1e-3], ...
     [0.8, 0.8, 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
     'DisplayName', 'Confidence Interval');
% Plot the confidence interval bounds using the purpleColor variable
%plot(kfit_ci_1 * 1e-3, '--', 'Color', purpleColor, 'LineWidth', 2, ...
%    'DisplayName', 'CI 1');
%plot(kfit_ci_2 * 1e-3, '--', 'Color', purpleColor, 'LineWidth', 2, ...
%    'DisplayName', 'CI 2');

% Add legend and labels
legend('Location', 'best');
xlabel('Frequency [Hz]', 'FontSize', 20, 'FontName', 'Arial');
ylabel('Wavenumber [1/mm]', 'FontSize', 20, 'FontName', 'Arial');
ylim([0,2]);
% Set font and other axis properties
set(gca, 'FontSize', 14, 'FontName', 'Arial');

hold off;

kfit_final = kfit * 1e-3;

disp(sprintf('Confidence Interval: %d', confidence_interval));
disp(sprintf('G: %.2f Pa', G_opt));
disp(sprintf('Mu: %.2f Pa', mu_opt));
disp(sprintf('Mu/G ratio: %.2f', mu_opt / G_opt));
end

function kfit_select = calc_niti_amode_loss(G, f, k, VdB, h, f_reduced_idx, mu_factor, rho, rho_l, c_l, cp)
mu = mu_factor*G;
% Calculate lambda
lambda = rho * cp^2 - 2 * mu;
% Compute fitted values
cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);

%plot(f,kfit);
kfit_select = kfit(f_reduced_idx)';

end
