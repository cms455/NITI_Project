function [G_opt, mu_opt, mu_opt_div_G, kfit_final] = fit_data_to_curve_no_constraint(rho, rho_l, c_l, cp, h, G0, mu0, mu_factor, f,k, VdB, f_reduced_idx)
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
    objfun = @(p) calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, p(1), p(1)*mu_factor, rho, rho_l, c_l, cp);

    % Perform optimization + Optimization Constraints
     opts = optimset('Display', 'iter', ...
        'TolX', 1e-0, ...         % Final convergence tolerance for the step size
        'TolFun', 1e-0);      % Final convergence tolerance for the function value

    figure; 
    hold on;
    plot(f(f_reduced_idx), VdB);
    title('Converging G Plot')
    [popt, fval_opt] = fminsearch(objfun,G0,opts);

    G_opt = popt(1);
    mu_opt = G_opt * mu_factor;
    mu_opt_div_G = mu_opt / G_opt;

    % Plot the optimized curve
    cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);
    kfit = f(:) ./ cfit(:);

    % Calculate original guess curve
    cfit_0 = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
    kfit_0 = f(:) ./ cfit_0(:);

    % Create a plot comparing the fitted curve with original guess and the data
    figure;
    hold on;
    plot(f(f_reduced_idx), VdB * 1e-3, 'b-', 'LineWidth', 2);
    scatter(f(f_reduced_idx), VdB * 1e-3, 'b');
    plot(f, kfit * 1e-3, 'r-', 'LineWidth', 2);
    plot(f, kfit_0 * 1e-3, 'g--', 'LineWidth', 2)
    legend('Measured Data', 'Data Points', 'Fitted Curve', 'Original Guess', 'Location','best');
    xlabel('Frequency (1/s)')
    ylabel('Wavenumber (1/m)')
    set(gca, 'FontSize', 20);

    hold off;

    kfit_final = kfit * 1e-3;

    disp(sprintf('G: %.2f Pa', G_opt));
    disp(sprintf('Mu: %.2f Pa', mu_opt));
    disp(sprintf('Mu/G ratio: %.2f', mu_opt / G_opt));
end

function loss = calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, G, mu, rho, rho_l, c_l, cp)
% Calculate lambda
lambda = rho * cp^2 - 2 * mu;
% Compute fitted values
cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
% Calculate loss (mean squared error)
loss = 0;
%Loop through each point of data and calculate the loss
for n = 2:length(VdB)
    f_idx = f_reduced_idx(n);
    loss = loss + (VdB(n) - kfit(f_idx))^2;
end

%plot each line to confirm that it is converging correctly. 
plot(f,kfit);
end