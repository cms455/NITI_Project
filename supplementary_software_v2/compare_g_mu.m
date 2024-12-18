
% Load example data, which contains the following arrays
% t  : [nt, 1] vector of time values in seconds 
% x  : [nx, 1] vector of lateral positions meters
% v  : [nx, nt] array giving the vertical surface vibration velocity, as measured with OCE
% dt : [scalar] temporal sampling interval in seconds
% dx : [scalar] spatial sampling interval in meters 
% h  : [scalar] corneal thickness (in meters) used in the analytical forward model
%load_data = importdata('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Cornea_NITI_water_mueq10G_fwP.txt');
load_data = importdata('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/550 um thickness/Cornea_NITI_water_mueq5G_fwP.txt');

data = load_data.data;

%load_ex_data = load('example_NITI_data.mat');
%test_data = load_ex_data.data_holder;



f_duplicates = data(:,1);
k_duplicates = data(:,2);
Vdb_flat = data(:,3);

f = unique(f_duplicates);
k = unique(k_duplicates);
VdB = reshape(Vdb_flat, length(f), length(k));
VdB = flipud(VdB');

f = f*1e3;

k = k*1e3;

fmax = 5000; % Hz
kmax = 2000; % 1/m
fmask = (f <= fmax);
kmask = (k <= kmax);
f = f(fmask);
k = k(kmask);
VdB = VdB(kmask, fmask);

%Input Parameters
h = 710e-6;
G0 = 1E3;
mu0 = 2E3;
reduce_flag = true;

G0_start = 10E3;
G0_end = 30E3;

mu0_start = 50E3;
mu0_end = 150E3;

G0_expected = 20E3;
mu0_expected = 100E3;

num_points = 5;

G0_set = linspace(G0_start, G0_end, num_points);
mu0_set = linspace(mu0_start, mu0_end, num_points);

results = zeros(num_points, num_points,3);

if reduce_flag
    select_freq = [600, 900, 1200, 1500];
    num_freq = length(select_freq);
    
    % Initialize f_reduced as an array of selected frequencies
    f_reg = f;
    f_reduced = zeros(1, num_freq);
    
    % Loop over each frequency in select_freq
    for i = 1:num_freq
        [~, idx] = min(abs(f - select_freq(i)));  
        f_reduced(i) = idx;  
    end
    
    disp(f(f_reduced));
    
    f = f(f_reduced);
    VdB= VdB(:,f_reduced);
end
    
for n =1:num_points
    G0 = G0_set(n);
    for m =1:num_points
        mu0 = mu0_set(m);
        [G, mu, fval] = fit_spectrum_niti(f, k, VdB, h, G0, mu0);
        [mu_iso, fval_iso] = fit_spectrum_iso(f, k, VdB, h, G0);
        results(n,m,:) = [G, mu, mu_iso];
    end
end

% Define the coolwarm colormap as a set called 'cmap'
n = 256; % Number of colors
r = [linspace(0, 1, n/2), ones(1, n/2)]';       % Red channel
g = [linspace(0, 1, n/2), linspace(1, 0, n/2)]'; % Green channel
b = [ones(1, n/2), linspace(1, 0, n/2)]';       % Blue channel

% Combine RGB channels into a colormap
cmap = [r, g, b];

% Plotting the results
[G0_mesh, mu0_mesh] = meshgrid(G0_set, mu0_set);
% Extract G, mu, mu_iso from results
G = results(:, :, 1);
mu = results(:, :, 2);
mu_iso = results(:, :, 3);

G_difference = G - G0_expected;
mu_difference = mu - mu0_expected;
mu_iso_difference = mu_iso - mu0_expected;


result_difference = [G_difference, mu_iso_difference, mu_difference];

max_val = max(results,[],'all');
max_val_difference = max(result_difference,[],'all');
min_val_difference = min(result_difference,[],'all');
% Plotting G as a heat map
figure;
imagesc(G0_set, mu0_set, G); % Use G0_set and mu0_set for axes
colorbar;
xlabel('G0');
ylabel('mu0');
clim([0, max_val]);
title('Heatmap of G');
set(gca, 'YDir', 'normal'); % Ensure the y-axis is oriented normally

% Plotting mu as a heat map
figure;
imagesc(G0_set, mu0_set, mu); % Use G0_set and mu0_set for axes
colorbar;
xlabel('G0');
ylabel('mu0');
clim([0, max_val]);
title('Heatmap of mu');
set(gca, 'YDir', 'normal');

% Plotting mu_iso as a heat map
figure;
imagesc(G0_set, mu0_set, mu_iso); % Use G0_set and mu0_set for axes
colorbar;
xlabel('G0');
ylabel('mu0');
clim([0, max_val]);
title('Heatmap of mu_{iso}');
set(gca, 'YDir', 'normal');

% Plotting G difference as a heat map
figure;
imagesc(G0_set, mu0_set, G_difference); % Use G0_set and mu0_set for axes
colorbar();
colormap(cmap)
xlabel('G0');
ylabel('mu0');
clim([min_val_difference, max_val_difference]);
title('Heatmap of G difference');
set(gca, 'YDir', 'normal');

% Plotting mu difference as a heat map
figure;
imagesc(G0_set, mu0_set, mu_difference); % Use G0_set and mu0_set for axes
colorbar();
colormap(cmap)
xlabel('G0');
ylabel('mu0');
clim([min_val_difference, max_val_difference]);
title('Heatmap of mu difference');
set(gca, 'YDir', 'normal');

% Plotting mu_iso difference as a heat map
figure;
imagesc(G0_set, mu0_set, mu_iso_difference); % Use G0_set and mu0_set for axes
colorbar;
colormap(cmap)
xlabel('G0');
ylabel('mu0');
clim([min_val_difference, max_val_difference]);
title('Heatmap of mu_{iso} difference');
set(gca, 'YDir', 'normal');
