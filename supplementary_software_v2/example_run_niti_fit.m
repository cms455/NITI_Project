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
%
load('example_onscale_results_niti_guided.mat')


%
% Compute the 2D Fourier power spectrum
% We need to use a large number of bins here because the bin size affects our fitting routine
% Smaller bins yield a better fit. Alternatively, one could downsample the data and use fewer
% bins.
%
nfft = 16384;
[VdB, f, k] = xttools_power_spectrum(v, dt, dx, nfft);


%
% Fit only to data over a limited frequency range
%
fmax = 3000; % Hz
kmax = 2000; % 1/m
fmask = (f <= fmax);
kmask = (k <= kmax);
f = f(fmask);
k = k(kmask);
VdB = VdB(kmask, fmask);

num_VdB = length(f)* length(k);
VdB_flat = reshape(VdB, 1, num_VdB);

data_holder = {f,k, VdB_flat};


%
% Run optimization - TI and models
% G0 and mu0 are initial guesses for the fit
% If either of these are not reasonably of the same order as the true moduli, the fit may not
% converge. If this happens, update the initial guess and re-run the fit function.
% When we call the isotropic fit, we use G0 as an initial guess since the A0 mode of a NITI
% material scales with G.
%
G0 = 7E3;
mu0 = 1E6;
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
figure
pcolor(f/1000, k/1000, VdB)
shading flat
caxis([-20, 0])
cb = colorbar;
ylabel(cb, 'Power Spectrum (dB)')
xlabel('Frequency (kHz)')
ylabel('Wavenumber (1/mm)')
hold on
h_ti = plot(f(:)/1000, kfit(:)/1000, '-r', 'LineWidth', 1);
h_iso = plot(f(:)/1000, kfit_iso(:)/1000, '-k', 'LineWidth', 1);
ylim([0, 1])
legend([h_ti, h_iso], 'NITI', 'Iso', 'Location', 'NorthWest')
legend('boxoff')