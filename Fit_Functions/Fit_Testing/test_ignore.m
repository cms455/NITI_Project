data = readmatrix('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Experimental_Data/RabbitCornea_IOP12.xlsx');

selected_freq = data(:,1);
VdB = data(:,2);
stdv = data(:,3);


% Generate F and K Ranges:
f_max = ceil(max(selected_freq));
f_min = floor(min(selected_freq));
f = linspace(f_min,f_max,(f_max-f_min)+1);

VdB = VdB *1e3;
k_max = ceil(max(VdB));
k_min = floor(min(VdB));
k = linspace(k_min, k_max, k_max-k_min+1);

num_freq = length(selected_freq);

% Initialize f_reduced as an array of selected frequencies
f_reduced = zeros(1, num_freq);

% Loop over each frequency in select_freq
for i = 1:num_freq
    [~, idx] = min(abs(f - selected_freq(i)));
    f_reduced(i) = idx;
end

%% Optimization Setup
% Constants and initial guesses
rho = 1000;    % Density of material (kg/m^3)
rho_l = 1000;  % Density of liquid (kg/m^3)
c_l = 1480;    % Speed of sound in liquid (m/s)
cp = 1540;     % Phase velocity (m/s)
h = 710E-6;    % Thickness (m)

f = linspace(600,2400,(2400-600));
k = linspace(200,1000, (1000-200));

G0 = 6E3;     % Initial guess for G (Pa)
mu0 = 50E3;    % Initial guess for mu (Pa)

lambda = rho * cp^2 - 2 * mu0; % Lame's constant (Pa)


cfit = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);

figure;
hold on;
plot(f,kfit);
plot(f(f_reduced_idx),VdB);

% Define the objective function
objfun = @(p) calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, p(1), p(2), rho, rho_l, c_l, cp);

% Perform optimization
opts = optimset('Display', 'iter');
[popt, fval_opt] = fminsearch(objfun, [G0, mu0], opts);

G_opt = popt(1);
mu_opt = popt(2);

%% Compute and Compare Fitted Curve

%G_opt = 20E3;
%mu_opt = 20E3;
cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
% Calculate difference between peaks
peak_difference = abs(max(VdB) - max(cfit));
disp(sprintf('Peak Difference: %.2f', peak_difference));


cfit_0 = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
kfit_0 = f(:)./cfit_0(:);

figure;
hold on;
plot(f(f_reduced), VdB*1e-3, 'b-', 'LineWidth', 2);
scatter(f(f_reduced), VdB*1e-3, 'b');
plot(f, kfit*1e-3, 'r-', 'LineWidth', 2);
plot(f, kfit_0*1e-3, 'g--', 'LineWidth', 2)
legend('VdB', 'Fitted Curve','Original Guess');


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

