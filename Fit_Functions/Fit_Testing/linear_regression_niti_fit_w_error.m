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


G0 = 50E3;     % Initial guess for G (Pa)
mu0 = 50E3;    % Initial guess for mu (Pa)

num_inits = 3;

G0_inits = linspace(10E3, 200E3, num_inits);
mu0_inits = linspace(10E3, 200E3, num_inits);

diff_holder = zeros(num_inits);
result_holder = zeros(num_inits,2);

figure;
hold on;

noise_levels = [0,0.1,0.2,0.5,1];

for n = 1:num_inits

mu0 = mu0_inits(n);
G0 = G0_inits(n);

lambda = rho * cp^2 - 2 * mu0; % Lame's constant (Pa)

% Define the objective function
objfun = @(p) calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, p(1), p(2), rho, rho_l, c_l, cp);

% Perform optimization
opts = optimset('Display', 'iter');
[popt, fval_opt] = fminsearch(objfun, [G0, mu0], opts);

G_opt = popt(1);
mu_opt = popt(2);

result_holder(n,1) = G_opt;
result_holder(n,2) = mu_opt;
%% Compute and Compare Fitted Curve

%G_opt = 20E3;
%mu_opt = 20E3;
cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
% Calculate difference between peaks
total_difference = sum(abs(VdB' - (kfit(f_reduced))));
disp(sprintf('total Difference: %.2f', total_difference));
diff_holder(n) = total_difference;
plot(f, kfit*1e-3, 'r-', 'LineWidth', 2);

end


cfit_0 = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
kfit_0 = f(:)./cfit_0(:);

plot(f(f_reduced), VdB*1e-3, 'b-o', 'LineWidth', 2);
plot(f, kfit_0*1e-3, 'g--', 'LineWidth', 2)
%legend( 'Fitted Curve','VdB','Original Guess');


hold off;

%% Display Results
disp(sprintf('G: %.2f Pa', G_opt));
disp(sprintf('Mu: %.2f Pa', mu_opt));
disp(sprintf('Mu/G ratio: %.2f', mu_opt / G_opt));

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
end

