%% Load and Process Data
% Import the data
file_path = '/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/550 um thickness/Cornea_NITI_water_mueqG_fwP.txt';
load_data = importdata(file_path);
data = load_data.data;

% Extract and organize data
f_duplicates = data(:, 1); % Frequencies (potential duplicates)
k_duplicates = data(:, 2); % Wavenumbers (potential duplicates)
VdB_flat = data(:, 3);     % VdB values

% Build VdB matrix
f = unique(f_duplicates) * 1e3; % Convert frequencies to Hz
k = unique(k_duplicates) * 1e3; % Convert wavenumbers to 1/m
VdB = reshape(VdB_flat, length(f), length(k));
VdB = flipud(VdB'); % Flip for correct orientation

% Apply frequency and wavenumber masks
fmax = 4000; % Maximum frequency (Hz)
kmax = 2000; % Maximum wavenumber (1/m)
fmask = (f <= fmax);
kmask = (k <= kmax);
f = f(fmask);
k = k(kmask);
VdB = VdB(kmask, fmask);

% Reduce VdB to specific frequencies
reduce_flag = true;
reduce_all_flag = true;
if reduce_flag
    %select_freq = [600, 900, 1200, 1500];
    select_freq = linspace(600, 1500, 8);
    num_freq = length(select_freq);

    % Initialize f_reduced as an array of selected frequencies
    f_reduced = zeros(1, num_freq);

    % Loop over each frequency in select_freq
    for i = 1:num_freq
        [~, idx] = min(abs(f - select_freq(i)));
        f_reduced(i) = idx;
    end

    disp(f(f_reduced));

    f_reduced_idx = f_reduced;
    
    VdB = VdB(:, f_reduced);
end

if reduce_all_flag
    reduced_k = zeros(1, num_freq);
    reduced_VdB = zeros(length(k), num_freq);

    % Arrays to store maximum peak locations and values
    max_peak_locations = zeros(1, num_freq);
    max_peak_values = zeros(1, num_freq);

    for m = 1:num_freq
        [pks, locs] = findpeaks(VdB(:, m));
        if ~isempty(pks)
            % Find the index of the maximum peak
            [max_peak_value, max_idx] = max(pks);
            % Store the location and value of the maximum peak
            max_peak_values(m) = max_peak_value;
            max_peak_locations(m) = locs(max_idx);
      
        end
    end
    VdB = (kmax/length(k))*max_peak_locations;
end


%% Optimization Setup
% Constants and initial guesses
rho = 1000;    % Density of material (kg/m^3)
rho_l = 1000;  % Density of liquid (kg/m^3)
c_l = 1480;    % Speed of sound in liquid (m/s)
cp = 1540;     % Phase velocity (m/s)
h = 710E-6;    % Thickness (m)


G0 = 30E3;     % Initial guess for G (Pa)
mu0 = 30E3;    % Initial guess for mu (Pa)


% Define noise levels and iterations
noise_levels = 0.2 * [0, 0.2, 0.5];
num_noise = length(noise_levels);
num_iterations = 3; % Number of iterations per noise level
result_noise_holder = zeros(num_noise, 2, num_iterations); % To store results

% Iterate over noise levels and random noise
for n = 1:num_noise
    for iter = 1:num_iterations
        % Add random noise to VdB
        noise = randn(size(VdB));
        VdB_noisy = VdB + noise_levels(n) * max(VdB, [], 'all') * noise;

        % Optimization
        objfun = @(p) calc_niti_amode_loss(f, k, VdB_noisy, h, f_reduced_idx, p(1), p(2), rho, rho_l, c_l, cp);
        opts = optimset('Display', 'off'); % Turn off iteration display
        [popt, ~] = fminsearch(objfun, [G0, mu0], opts);

        % Store results
        result_noise_holder(n, 1, iter) = popt(1); % G
        result_noise_holder(n, 2, iter) = popt(2); % mu
    end
end

% Calculate mean and standard deviation for each noise level
G_noise_set = mean(result_noise_holder(:, 1, :), 3);
mu_noise_set = mean(result_noise_holder(:, 2, :), 3);
G_std = std(result_noise_holder(:, 1, :), 0, 3);
mu_std = std(result_noise_holder(:, 2, :), 0, 3);

% Calculate G-coefficients and their standard deviations
G_coeffs = mu_noise_set ./ G_noise_set;
G_coeffs_std = G_coeffs .* sqrt((mu_std ./ mu_noise_set).^2 + (G_std ./ G_noise_set).^2);

% Plot bar plot with error bars for G and mu values
figure;
hold on;
bar([G_noise_set, mu_noise_set], 'grouped');
errorbar(1:num_noise, G_noise_set, G_std, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
errorbar(1:num_noise + 0.2, mu_noise_set, mu_std, 'k.', 'LineWidth', 1.5, 'CapSize', 10);
legend({'G Values', '\mu Values'}, 'Location', 'Best');
xlabel('Noise Level Index');
ylabel('Values');
title('G and \mu Values with Error Bars');
grid on;
hold off;

% Plot line plot with error bars for G-coefficients
figure;
errorbar(1:num_noise, G_coeffs, G_coeffs_std, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Noise Level Index');
ylabel('G Coefficients (\mu / G)');
title('G Coefficients with Error Bars');
grid on;



%% Functions
function loss = calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, G, mu, rho, rho_l, c_l, cp)
    % Calculate Lame's constant
    lambda = rho * cp^2 - 2 * mu;

    % Compute fitted values
    cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
    kfit = f(:)./cfit(:);
    % Calculate loss (mean squared error)
    
    loss = 0;
    for n = 1:length(VdB)
        f_idx = f_reduced_idx(n);
        loss = loss + (VdB(n) - kfit(f_idx))^2;
    end
    a = 0;
end

