%DESCRIPTION
% Takes in a matrix of data of shape(f,k,stdv), and fits a dispersion
% curve, finding the optimized G and Mu. Then it plots the curve with the
% data. With the constraints that G must be positive > 1e3.

%load the data
data = readmatrix('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Experimental_Data/RabbitCornea_IOP12.xlsx');

%recover the selected frequences, k_values at selected frequencies, and
%stdv
selected_freq = data(:,1);
VdB = data(:,2);
stdv = data(:,3);


% Generate F and K Ranges:
f_max = ceil(max(selected_freq));
f_min = floor(min(selected_freq));
f = linspace(f_min,f_max,(f_max-f_min)+1);

VdB = VdB *1e3;
k_max = ceil(max(VdB)) + 5e3;
k_min = floor(min(VdB));
k = linspace(k_min, k_max, k_max-k_min+1);

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

%% Optimization Setup
% Constants and initial guesses
rho = 1000;    % Density of material (kg/m^3)
rho_l = 1000;  % Density of liquid (kg/m^3)
c_l = 1480;    % Speed of sound in liquid (m/s)
cp = 1540;     % Phase velocity (m/s)
h = 800E-6;    % Thickness (m)

G0 = 20E3;     % Initial guess for G (Pa)
mu_factor = 110; % mu = mu_factor*G;
mu0 = G0*mu_factor;    % Initial guess for mu (Pa)
lambda = rho * cp^2 - 2 * mu0; % Lame's constant (Pa)

% Define the objective function
lb = [1e3]; % G must be positive (Weird stuff happens for G < 1e3)
ub = []; % No upper bound in this case
objfun = @(p) calc_niti_amode_loss(f, k, VdB, h, f_reduced_idx, p(1), p(1)*mu_factor, rho, rho_l, c_l, cp);

% Perform optimization + Optimization Constratins
opts = optimoptions('fmincon', ...
    'Display', 'iter', ...
    'Algorithm', 'interior-point', ...
    'StepTolerance', 1e-6, ...
    'TolFun', 1e-6, ...
    'TolX', 1e-6);
figure; 
hold on;
plot(f(f_reduced_idx),VdB);
[popt, fval_opt] = fmincon(objfun, G0, [], [], [], [], lb, ub, [], opts);

G_opt = popt(1);
mu_opt = G_opt*mu_factor;

%% Compute and Compare Fitted Curve
%Plot the optimized curve
cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
%Calculate original guess curve
cfit_0 = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
kfit_0 = f(:)./cfit_0(:);

%Create a plot comparing the fitted curve with original guess and the data
figure;
hold on;
plot(f(f_reduced), VdB*1e-3, 'b-', 'LineWidth', 2);
scatter(f(f_reduced), VdB*1e-3, 'b');
plot(f, kfit*1e-3, 'r-', 'LineWidth', 2);
plot(f, kfit_0*1e-3, 'g--', 'LineWidth', 2)
legend('Measured Data', 'Data Points','Fitted Curve','Original Guess');
xlabel('Frequency (1/s)')
ylabel('Wavenumber (1/m)')
set(gca, 'FontSize',20);


hold off;

%% Display Results
disp(sprintf('G: %.2f Pa', G_opt));
disp(sprintf('Mu: %.2f Pa', mu_opt));
disp(sprintf('Mu/G ratio: %.2f', mu_opt / G_opt));

%% Functions
%Optimization Function 
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

%plot each line to confirm that it is converging correctly. 
plot(f,kfit);
end

