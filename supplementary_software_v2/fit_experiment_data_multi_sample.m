
%load the data
data = readmatrix('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Experimental_Data/kMeans_subject6.csv');

%recover the selected frequences, k_values at selected frequencies, and
%stdv
selected_freq = [944, 1335, 1660, 2115];
VdB_samples = data();
VdB_samples(5,:) = [];




% Generate F and K Ranges:
f_max = ceil(max(selected_freq));
f_min = floor(min(selected_freq));
f = linspace(f_min- 400,f_max + 400,(f_max-f_min)+1);

VdB_samples = VdB_samples *1e3;
k_max = ceil(max(VdB_samples,[],'all')) + 10e3;
k_min = floor(min(VdB_samples,[],'all'));
k = linspace(0, k_max, k_max-k_min+1);

% find the indexes of the selected frequences for linspace(f)
num_freq = length(selected_freq);
% initialize f_reduced as an array of selected frequencies
f_reduced = zeros(1, num_freq);
% Loop over each frequency in select_freq
for i = 1:num_freq
    [~, idx] = min(abs(f - selected_freq(i)));
    f_reduced(i) = idx;
end
f_reduced_idx = f_reduced;
num_samples = size(VdB_samples,1);

%PARAMS
rho = 1000;    % Density of material (kg/m^3)
rho_l = 1000;  % Density of liquid (kg/m^3)
c_l = 1480;    % Speed of sound in liquid (m/s)
cp = 1540;     % Phase velocity (m/s)
h = 800E-6;    % Thickness (m)
G0 = 40E3;     % Initial guess for G (Pa)
mu_factor = 110; % mu = mu_factor*G;
mu0 = G0*mu_factor;    % Initial guess for mu (Pa)

kfit_holder = cell(num_samples);
G_opt_holder = zeros(1,num_samples);
mu_opt_holder = zeros(1,num_samples);

for s = 1:num_samples
    VdB = squeeze(VdB_samples(s,:));
    [G_opt, mu_opt, mu_opt_div_G, kfit_final]= fit_data_to_curve(rho, rho_l, c_l, ...
        cp, h, G0, mu0, mu_factor, f,k,VdB,f_reduced_idx);
    kfit_holder{s} = kfit_final;
    G_opt_holder(s) = G_opt;
    mu_opt_holder(s) = mu_opt;
end

figure;
hold on;
for s = 1:num_samples
    kfit_final = kfit_holder{s};
    plot(kfit_final, 'DisplayName',['Sample ', num2str(s)]);
end
title('Dispersion Curves for Subject')
xlabel('Frequency (1/s)')
ylabel('Wave-number (1/m)')
set(gca, 'FontSize', 14);
