%
%
% example_run_niti_fit
% 
% This example script runs the NITI fitting routine on a set of example OCE data to 
% estimate the moduli G and mu. The comments in the script and the called functions
% explain their inputs, outputs, and behavior.
%
% The output of this script is an f-k plot with overlaid best-fit NITI and isotropic
% A0 dispersion curves, similar to Figure 5b in the paper (though this is a different
% measurement than the one shown in that figure).
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
%


%
% Load example data, which contains the following arrays
% t  : [nt, 1] vector of time values in seconds 
% x  : [nx, 1] vector of lateral positions meters
% v  : [nx, nt] array giving the vertical surface vibration velocity, as measured with OCE
% dt : [scalar] temporal sampling interval in seconds
% dx : [scalar] spatial sampling interval in meters 
% h  : [scalar] corneal thickness (in meters) used in the analytical forward model
load_data = importdata('/Users/calvinsmith/Bouma_lab/NITI_project/Dispersion_Data/Cornea_iso_vacuum_Geqmu_fwP.txt');
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
%}

%f = test_data{1};
%k = test_data{2};
%VdB_flat = test_data{3};
%VdB = reshape(VdB_flat, length(k), length(f));





%
% Compute the 2D Fourier power spectrum
% We need to use a large number of bins here because the bin size affects our fitting routine
% Smaller bins yield a better fit. Alternatively, one could downsample the data and use fewer
% bins.
%
%nfft = 16384;
%[VdB, f, k] = xttools_power_spectrum(v, dt, dx, nfft);


%
% Fit only to data over a limited frequency range
%



%
% Run optimization - TI and models
% G0 and mu0 are initial guesses for the fit
% If either of these are not reasonably of the same order as the true moduli, the fit may not
% converge. If this happens, update the initial guess and re-run the fit function.
% When we call the isotropic fit, we use G0 as an initial guess since the A0 mode of a NITI
% material scales with G.
%
%f = f';
f = f*1e3;
%k= k';
k = k*1e3;

%VdB = fliplr(VdB);

fmax = 5000; % Hz
kmax = 2000; % 1/m
fmask = (f <= fmax);
kmask = (k <= kmax);
f = f(fmask);
k = k(kmask);
VdB = VdB(kmask,fmask);
%}
%Vdb = Vdb';


h = 710e-6;

G0 = 14E3;
mu0 = 14E3;
[G, mu, fval] = fit_spectrum_niti(f, k, VdB, h, G0, mu0);
[mu_iso, fval_iso] = fit_spectrum_iso(f, k, VdB, h, G0);

[phi_max, ~] = global_mode_energy(f, k, VdB, 1);
gof_ti = abs(fval)/phi_max;
gof_iso = abs(fval_iso)/phi_max;


%
% Re-calculate NITI best fit dispersion curve
% The following assumptions are used inside the compute_niti_amode function.
% rho   : [scalar] assumed corneal tissue density 
% rho_l : [scalar] density of the liquid bounding the cornea (assume water)
% c_l   : [scalar] P-wave speed of the liquid bounding the cornea (assume water)
% cp    : [scalar] P-wave speed of the corneal tissue
% lambda : [scalar] stiffness tensor entry that enforces incompressibility, lambda >> mu > G
%
rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;
lambda = rho*cp^2 - 2*mu;
cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);
kfit = reshape(kfit, 1, []);


%
% Re-calculate isotropic best fit dispersion curve
% The following assumptions are used inside the fit_spectrum_niti function.
% rho   : [scalar] assumed corneal tissue density 
% rho_l : [scalar] density of the liquid bounding the cornea (assume water)
% c_l   : [scalar] P-wave speed of the liquid bounding the cornea (assume water)
% cp    : [scalar] P-wave speed of the corneal tissue
% lambda : [scalar] stiffness tensor entry that enforces incompressibility, lambda >> mu > G
%
rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;
lambda = rho*cp^2 - 2*mu_iso;
cfit_iso = compute_niti_amode(f, k, h, mu_iso, mu_iso, lambda, rho, rho_l, c_l);
kfit_iso = f(:)./cfit_iso(:);
kfit_iso = reshape(kfit_iso, 1, []);


%
% Plot results
%

save('kfit_data.mat','kfit');
save('kfit_iso_data.mat','kfit_iso');
%{
load_kfit = load('kfit_data.mat');
k_fit = load_kfit.kfit;
load_kfit_iso = load('kfit_iso_data.mat');
k_fit_iso = load_kfit_iso.kfit_iso;
%}

figure
pcolor(f/1000,k/1000, (VdB));
shading flat
%caxis([-70, -45])
cb = colorbar;
ylabel(cb, 'Power Spectrum (dB)')
xlabel('Frequency (kHz)')
ylabel('Wavenumber (1/mm)')
hold on
h_ti = plot(f(:)/1000, kfit(:)/1000, '-r', 'LineWidth', 1);
h_iso = plot(f(:)/1000, kfit_iso(:)/1000, '-k', 'LineWidth', 1);
%ylim()
legend([h_ti, h_iso], 'NITI', 'Iso', 'Location', 'NorthWest')
legend('boxoff')
colormap;
