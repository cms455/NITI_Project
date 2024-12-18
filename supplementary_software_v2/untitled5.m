
load_data = importdata('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/550 um thickness/Cornea_NITI_water_mueqG_fwP.txt');

data = load_data.data;

% Build VdB
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



figure;
plot(VdB);
G  = 20E3;
mu = 20E3;
h = 710E-6;

f =  linspace(1,5000,50);
k = linspace(1,2000,50);
% Optimization using the simplex method
opts = optimset('Display', 'iter');
objfun = @(p) calc_niti_amode_loss(f,k, VdB,h, f_reduced_idx, p(1), p(2));
[popt, fval_opt] = fminsearch(objfun, [G0, mu0], opts);
G_opt = popt(1);
mu_opt = popt(2);


cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);

figure;
plot(cfit);

disp(sprintf('G: %d', G_opt));
disp(sprintf('Mu: %d', mu_opt));
disp(sprintf('mu = %d * G', mu_opt/G_opt));




function loss = calc_niti_amode_loss(f,k, VdB,h, f_reduced_idx, G,mu)
rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;

cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);

loss = 0;
num_freq = length(f);

for n = 1:num_freq
   f_idx = f_reduced_idx(n);
   loss = loss + (VdB(n) - cfit(f_idx))^2;
end


end

