% Constants and initial guesses
rho = 1000;    % Density of material (kg/m^3)
rho_l = 1000;  % Density of liquid (kg/m^3)
c_l = 1480;    % Speed of sound in liquid (m/s)
cp = 1540;     % Phase velocity (m/s)
h = 710E-6;    % Thickness (m)

G0 = 30E3;     % Initial guess for G (Pa)
mu0 = 30E3;    % Initial guess for mu (Pa)
num_selected_f_points = 6;

% Frequency and wavenumber settings
fmax = 4000; % Maximum frequency (Hz)
kmax = 2000; % Maximum wavenumber (1/m)
num_f_points = 100;
num_k_points = 100;
f = linspace(0, fmax, num_f_points);
k = linspace(0, kmax, num_k_points);
lambda = rho * cp^2 - 2 * mu0; % Lame's constant (Pa)

% Compute initial VdB
cfit = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
kfit = f(:) ./ cfit(:);
VdB = kfit;

% Reduce VdB to specific frequencies
select_freq = linspace(600, 1500, num_selected_f_points);
num_freq = length(select_freq);
f_reduced = zeros(1, num_freq);

for i = 1:num_freq
    [~, idx] = min(abs(f - select_freq(i)));
    f_reduced(i) = idx;
end
f_reduced_idx = f_reduced;
VdB = VdB(f_reduced);

% Noise settings
noise_levels = 0.5 * [0, 0.1, 0.2, 0.5];
num_noise = length(noise_levels);
num_iterations = 1;
noise = randn(size(VdB));
VdB_og = VdB;

% Optimization bounds
lb = [0, 0];        % Lower bounds for G and mu
ub = [1e7, Inf];    % Upper bounds for G and no constraint for mu
options = optimoptions('lsqcurvefit');

% Storage for results
result_noise_holder = zeros(num_iterations, num_noise, 2);

diff_holder = zeros(num_noise);

figure;
hold on;

for m = 1:num_iterations
    for n = 1:num_noise
        % Add noise to VdB
        VdB = VdB_og + noise_levels(n) * max(VdB_og) * noise;

        % Define residual function for lsqcurvefit
        residual_fun = @(p, xdata) calc_residuals(p, xdata);

        % Perform optimization
        [popt, resnorm] = lsqcurvefit(residual_fun, [G0, mu0], struct('f', f, 'k', k, 'VdB', VdB, 'h', h, 'f_reduced_idx', f_reduced_idx, 'rho', rho, 'rho_l', rho_l, 'c_l', c_l, 'cp', cp), VdB, lb, ub, options);
        G_opt = popt(1);
        mu_opt = popt(2);

        % Store results
        result_noise_holder(m, n, 1) = G_opt;
        result_noise_holder(m, n, 2) = mu_opt;

        
        plot(f, kfit*1e-3, 'r-', 'LineWidth', 2);
        
        plot(f(f_reduced), VdB*1e-3, 'b-o', 'LineWidth', 2);
        plot(f, kfit_0*1e-3, 'g--', 'LineWidth', 2)

        % Display results
        fprintf('Iter %d, Noise %d\n', m, n);
        fprintf('G: %.2f, Mu: %.2f\n', G_opt, mu_opt);
    end
end

% Post-processing: calculate mean and standard deviations
G_noise_mean = mean(result_noise_holder(:, :, 1), 1);
mu_noise_mean = mean(result_noise_holder(:, :, 2), 1);
G_noise_std = std(result_noise_holder(:, :, 1), 0, 1);
mu_noise_std = std(result_noise_holder(:, :, 2), 0, 1);

% Plot results
figure;
errorbar(1:num_noise, G_noise_mean, G_noise_std, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Noise Level Index');
ylabel('G Values (Pa)');
title('G Values with Error Bars Across Noise Levels');
grid on;

figure;
errorbar(1:num_noise, mu_noise_mean, mu_noise_std, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Noise Level Index');
ylabel('\mu Values (Pa)');
title('\mu Values with Error Bars Across Noise Levels');
grid on;


% Residual function for lsqcurvefit
function residuals = calc_residuals(p, xdata) % Adjust function signature to take two inputs
    % Extract parameters and data from xdata
    f = xdata.f;
    k = xdata.k;
    VdB = xdata.VdB;
    h = xdata.h;
    f_reduced_idx = xdata.f_reduced_idx;
    rho = xdata.rho;
    rho_l = xdata.rho_l;
    c_l = xdata.c_l;
    cp = xdata.cp;

    % Extract parameters to optimize
    G = p(1);
    mu = p(2);

    % Calculate lambda
    lambda = rho * cp^2 - 2 * mu;

    % Compute fitted values
    cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
    kfit = f(:) ./ cfit(:);

    % Compute residuals
    residuals = abs(VdB(:) - kfit(f_reduced_idx(:)));
end
