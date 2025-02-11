%load the data
data = readmatrix('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Experimental_Data/kMeans_Rab.csv');
subject ='6';

%recover the selected frequences, k_values at selected frequencies, and
%stdv
selected_freq = [944, 1335, 1660, 2115];

VdB_samples = data;
%VdB_samples(5,:) = [];
num_samples = size(VdB_samples,1);

VdB_samples_left = (2/(num_samples))*sum(VdB_samples(1:4,:),1);
VdB_samples_right = (2/(num_samples))*sum(VdB_samples(5:8,:),1);
VdB_samples = cat(1, VdB_samples_left, VdB_samples_right);

error_data = readmatrix('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Experimental_Data/kStdv_Rab.csv');
%error_data(5,:) = [];
error_data_left = (2/(num_samples))*sum(error_data(1:4,:),1);
error_data_right = (2/(num_samples))*sum(error_data(5:8,:),1);
error_data = cat(1, error_data_left, error_data_right);

% Generate F and K Ranges:
f_max = ceil(max(selected_freq));
f_min = floor(min(selected_freq));
f = linspace(0,f_max,(f_max)+1);

VdB_samples = VdB_samples *1e3;
k_max = ceil(max(VdB_samples,[],'all')) + 10e3;
k_min = floor(min(VdB_samples,[],'all'));
k = linspace(0, k_max, k_max);

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
mu_factor = 110/3;

kfit_holder = cell(num_samples);
G_opt_holder = zeros(1,num_samples);
mu_opt_holder = zeros(1,num_samples);
confidence_interval_holder = cell(1,num_samples);

for s = 1:num_samples
    VdB = squeeze(VdB_samples(s,:));
    VdB_error = squeeze(error_data(s,:));
    G_list = [5e3, 10e3, 15e3, 20e3, 25e3, 30e3, 35e3, 40e3,50e3, 60e3];
    G_loss = zeros(1,length(G_list));
    
    for i = 1:length(G_list)
        G_loss(i)  = calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, G_list(i), mu_factor*G_list(i) , rho, rho_l, c_l, cp);
        disp(sprintf('Loss of %d: %d', G_list(i), G_loss(i)));
    end
    
    [~, G_idx] = min(G_loss);
    G0 = G_list(G_idx);
    mu0 = G0*mu_factor;    % Initial guess for mu (Pa)

    [G_opt, mu_opt, mu_opt_div_G, kfit_final, confidence_interval]= fit_data_to_curve_lsq(rho, rho_l, c_l, ...
        cp, h, G0, mu0, mu_factor, f,k,VdB,f_reduced_idx, VdB_error);
    kfit_holder{s} = kfit_final;
    G_opt_holder(s) = G_opt;
    mu_opt_holder(s) = mu_opt;
    confidence_interval_holder{s} = confidence_interval;
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

id = round(1000*rand(1));

save_folder  = '/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Saved_Results/Final_Results/';

save_k_name = ['kfit_holder_',subject,'_', num2str(id)];
save_k_path = [save_folder,save_k_name];

save_G_name = ['G_opt_holder_',subject, '_', num2str(id)];
save_G_path = [save_folder,save_G_name];

save_ci_name = ['confidence_interval_holder_',subject, '_', num2str(id)];
save_ci_path = [save_folder,save_ci_name];

save(save_ci_path, 'confidence_interval_holder');
save(save_k_path, 'kfit_holder');
save(save_G_path,'G_opt_holder');


function loss = calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, G, mu, rho, rho_l, c_l, cp)
% Calculate lambda
lambda = rho * cp^2 - 2 * mu;
% Compute fitted values
cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
% Calculate loss (mean squared error)
loss = 0;
%Loop through each point of data and calculate the loss
for n = 1:length(VdB)
    f_idx = f_reduced_idx(n);
    loss = loss + (VdB(n) - kfit(f_idx))^2;
end
end


plot_fit_data_v2(kfit_holder, G_opt_holder, confidence_interval_holder, mu_factor, subject, VdB_samples,f,k, f_reduced_idx);