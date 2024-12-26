
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
    select_freq = linspace(600, 1500, 4);
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



% Constants and initial guesses
rho = 1000;    % Density of material (kg/m^3)
rho_l = 1000;  % Density of liquid (kg/m^3)
c_l = 1480;    % Speed of sound in liquid (m/s)
cp = 1540;     % Phase velocity (m/s)
h = 710E-6;    % Thickness (m)


G0 = 30E3;     % Initial guess for G (Pa)
mu0 = 30E3;    % Initial guess for mu (Pa)


noise_levels = 0.2*[0,0.1,0.2,0.5];
num_noise = length(noise_levels);
num_iterations = 5;
noise = randn(size(VdB));
lambda = rho * cp^2 - 2 * mu0; % Lame's constant (Pa)
VdB_og = VdB;

diff_holder = zeros(num_noise);
result_noise_holder = zeros(num_iterations,num_noise,2);


for m = 1:num_iterations
for n =1:num_noise

figure;
hold on;
% Define the objective function
noise = randn(size(VdB));
VdB = VdB_og;
VdB = VdB + (noise_levels(n)*max(VdB,[],'all')*noise);

c = polyfit(f_reduced,VdB,1);

x = f_reduced;
y_est = polyval(c,x);
% Add trend line to plot
%hold on
plot(f(f_reduced),y_est*1e-3,'r--','LineWidth',2);
scatter(f(f_reduced), VdB*1e-3, 'bx', 'LineWidth', 4);
%hold off
%VdB = y_est;
objfun = @(p) calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, p(1), p(2), rho, rho_l, c_l, cp);

% Perform optimization
%opts = optimset('Display', 'iter');
[popt, fval_opt] = fminsearch(objfun, [G0, mu0]);

G_opt = popt(1);
mu_opt = popt(2);

result_noise_holder(m,n,1) = G_opt;
result_noise_holder(m,n,2) = mu_opt;
disp(sprintf('Iter %d, Noise %d', m, n));
disp(sprintf('G: %d', G_opt));
disp(sprintf('Mu: %d', mu_opt));
disp(sprintf('Mu/G ratio: %d', mu_opt / G_opt));




cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
% Calculate difference between peaks
total_difference = sum(abs(VdB' - (kfit(f_reduced))));
diff_holder(n) = total_difference;
disp(sprintf('total Difference: %.2f \n', total_difference));
plot(f, kfit*1e-3, 'LineWidth', 2);


end

figure;
hold on;

cfit_0 = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
kfit_0 = f(:)./cfit_0(:);

plot(f, kfit*1e-3, 'r-', 'LineWidth', 2);

plot(f(f_reduced), VdB*1e-3, 'b-o', 'LineWidth', 2);
plot(f, kfit_0*1e-3, 'g--', 'LineWidth', 2)
%legend( 'Fitted Curve','VdB','Original Guess');


hold off;

% Assuming result_noise_holder is defined, with dimensions [noise_levels, 2]

% Extract G and mu values from result_noise_holder
G_noise_set = squeeze(result_noise_holder(m,:, 1));
mu_noise_set = squeeze(result_noise_holder(m,:, 2));
% Bar plot for G and mu values with different colors
figure;

% Create grouped bar plot
b = bar([G_noise_set', mu_noise_set'], 'grouped');

% Set different colors for G and mu bars
b(1).FaceColor = 'b'; % Blue for G values
b(2).FaceColor = 'r'; % Red for mu values

% Add legend, labels, and title
legend({'G Values', '\mu Values'}, 'Location', 'Best');
xlabel('Noise Level Index');
ylabel('Values (Pa)');
title('G and \mu Values for Each Noise Level');
grid on;


% Calculate G-coefficients
G_coeffs = mu_noise_set ./ G_noise_set;

% Line plot for G-coefficients
figure;
plot(G_coeffs, '-o', 'LineWidth', 2, 'MarkerSize', 8); % Line plot with markers
xlabel('Noise Level Index');
ylabel('G Coefficients (\mu / G)');
title('G Coefficients for Each Noise Level');
grid on;

end
result_noise_holder(result_noise_holder > 1e6) = 1e6;
% Calculate mean and standard deviation for G and mu across iterations
G_noise_mean = mean(result_noise_holder(:, :, 1), 1);
mu_noise_mean = mean(result_noise_holder(:, :, 2), 1);

G_noise_std = std(result_noise_holder(:, :, 1), 0, 1);
mu_noise_std = std(result_noise_holder(:, :, 2), 0, 1);

% Calculate G-coefficients and their standard deviations
G_coeffs = mu_noise_mean ./ G_noise_mean;
G_coeffs_std = G_coeffs .* sqrt((mu_noise_std ./ mu_noise_mean).^2 + (G_noise_std ./ G_noise_mean).^2);


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


figure;
errorbar(1:num_noise, G_coeffs, G_coeffs_std, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Noise Level Index');
ylabel('G Coefficients (\mu / G)');
title('G Coefficients with Error Bars Across Noise Levels');
grid on;


% Calculate the differences from initial guesses
G_diff = abs(G_noise_mean - G0);
mu_diff = abs(mu_noise_mean - mu0);

% Plot the differences
figure;
hold on;

% Plot G differences (blue)
plot(1:num_noise, G_diff, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'b', 'DisplayName', '|G_{mean} - G_0|');

% Plot Mu differences (red)
plot(1:num_noise, mu_diff, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'r', 'DisplayName', '|\mu_{mean} - \mu_0|');

% Add legend, labels, and title
legend('show', 'Location', 'Best');
xlabel('Noise Level Index');
ylabel('Absolute Difference (Pa)');
title('Differences Between Calculated Mean Values and Initial Guesses');
grid on;
hold off;


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

