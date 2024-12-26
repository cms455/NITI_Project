function [G_opt, mu_opt, fval_opt] = fit_spectrum_niti(f, k, VdB, h, G0, mu0, w)
% FIT_SPECTRUM_NITI - Find the NITI A0 mode that matches measured spectrum
%
% [G_opt, mu_opt, fval_opt] = fit_spectrum_niti(f, k, VdB, h, G0, mu0, w)
%
% Given the 2D Fourier power spectrum VdB(f,k) for a cornea/phantom of thickness h, 
% estimate the best-fit values of G and mu based on the A0 mode of a NITI material. 
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
% Also outputs a figure showing the fit overlaid on  the normalized power spectrum (linear scale).
%
% Notes
% -----
% - Optimization is performed using a simplex method.
% - The sampling used in computing VdB (number of FFT points) will affect the accuracy
%   of the fit. Be sure that enough points are used. For the NITI paper, I used nfft = 16384.
%   Alternatively, downsample the XT plot prior to computing the FFT.
% - The fit will fail if the initial guess is too far from the true value.
% - It is mathematically possible for the fit to return negative values of the constants, but
%   these are non-physical. Re-run with a different inital guess if this happens.
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
%
% ---

    if nargin < 7
        w = 1;
    end

    % Optimization using the simplex method
    opts = optimset('Display', 'iter');
    objfun = @(p) -ti_mode_energy(f, k, VdB, w, h, p(1), p(2));
    [popt, fval_opt] = fminsearch(objfun, [G0, mu0], opts);
    G_opt = popt(1);
    mu_opt = popt(2);
    
    % Generate plot of fit result in the f-k domain
    %plot_fit_result(f, k, VdB, w, h, G_opt, mu_opt)
    
    
    
function cost = ti_mode_energy(f, k, VdB, w, h, G, mu)
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
    c = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);

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
    
    
function plot_fit_result(f, k, VdB, w, h, G, mu)
    rho = 1000;
    rho_l = 1000;
    c_l = 1480;
    cp = 1540;
    lambda = rho*cp^2 - 2*mu;
    cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
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
