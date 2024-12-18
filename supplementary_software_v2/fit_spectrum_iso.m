function [mu_opt, fval_opt] = fit_spectrum_iso(f, k, VdB, h, mu0, w)
% FIT_SPECTRUM_ISO - Find the isotropic A0 mode that matches the measured spectrum
%
% [mu_opt, fval_opt] = fit_spectrum_iso(f, k, VdB, h, mu0, w)
%
% This is identical to fit_spectrum_niti, except it enforces the condition
% G = mu. This function fits ONLY the A0 mode. For isotropic phantoms, we typically
% see both A0 and S0 in the spectrum. To fit to truly isotropic data, this function
% is a poor choice, as both A0 and S0 will give a more accurate result.
%
% Parameters
% ----------
% f   : [double, nf] array of frequency values in Hz
% k   : [double, nk] array of wavenumber values in 1/m
% VdB : [double, nk x nf] 2D power spectrum in dB (rows = k, cols = f)
% h   : [double] phantom or cornea thickness
% G0  : [double] initial guess of G
% mu0 : [double] initial guess of mu
% w   : [double] (optional) array of weights to apply to power spectrum
%
% Returns
% -------
% G_opt    : [double] optimal estimate of G
% mu_opt   : [double] optimal estimate of mu
% fval_opt : [double] final value of objective function
%
% Also outputs a figure showing the fit overlaid on the normalized power spectrum (linear scale).
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
%
% ---

    if nargin < 6
        w = 1;
    end

    % Optimization using the simplex method
    opts = optimset('Display', 'iter');
    objfun = @(p) -ti_mode_energy(f, k, VdB, w, h, p);
    [popt, fval_opt] = fminsearch(objfun, mu0, opts);
    mu_opt = popt;
    
    % Generate plot of fit result in the f-k domain
    % plot_fit_result(f, k, VdB, w, h, mu_opt)
    
    
    
function cost = ti_mode_energy(f, k, VdB, w, h, mu)
% Compute the energy covered by the A0 mode (weighted by a 7-point Gaussian window)
% This is the value we want to maximize for fitting the A0 mode to the 2D Fourier 
% spectrum.
%
% ---
    % These are assumed constants
    rho = 1000; % Density of the corneal tissue in kg/m^3
    rho_l = 1000; % Density of the liquid bounding the cornea in kg/m^3 (assume water)
    c_l = 1480; % P-wave speed of the liquid bounding the cornea in m/s (assume water)
    cp = 1540; % P-wave speed of the corneal tissue in m/w
    lambda = rho*cp^2 - 2*mu; % stiffness matrix modulus that enforces incompressibility
    
    % Gaussian window
    gwin = gausswin(7);
    gwin = gwin./sum(gwin);

    % Convert VdB to normalized energy
    energy = w.*(10.^(VdB/10));
    dk = k(2) - k(1);
    
    % Forward model solve for the A0 mode phase velocity, c(f)
    % Call NITI function with G = mu to give isotropic response
    c = compute_niti_amode(f, k, h, mu, mu, lambda, rho, rho_l, c_l);

    % Compute energy in the f-k space under the A0 mode
    ki = f(:)./c(:);
    ki = reshape(ki, 1, []);
    idcs = transpose(-3:3) + round(ki/dk);
    energy_sum = 0;
    sum_num = 0;
    for idx = 1:length(f)
        if all(isfinite(idcs(:,idx))) && all(idcs(:,idx) > 0) && all(idcs(:,idx) <= length(k))
            energy_sum = energy_sum + sum(gwin.*abs(energy(idcs(:,idx), idx)));
            sum_num = sum_num + 1;
        end
    end
    beta = 1;
    regval = abs(mu/lambda);
    cost = energy_sum/sum_num - beta*regval;
    
    
function plot_fit_result(f, k, VdB, w, h, mu)
    rho = 1000;
    rho_l = 1000;
    c_l = 1480;
    cp = 1540;
    lambda = rho*cp^2 - 2*mu;
    cfit = compute_niti_amode(f, k, h, mu, mu, lambda, rho, rho_l, c_l);
    kfit = f(:)./cfit(:);
    kfit = reshape(kfit, 1, []);
    
    figure
    set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
    imagesc(f(:)', k(:), w.*(10.^(VdB/10)))
    hold on
    plot(f, kfit, '-r', 'LineWidth', 1)
    xlabel('Frequency (Hz)')
    ylabel('Wavenumber (1/m)')
    cb = colorbar;
    ylabel(cb, 'Normalized Spectral Power')
    caxis([0, 1])
