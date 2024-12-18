function [cost, kopt] = global_mode_energy(f, k, VdB, w)
% GLOBAL_MODE_ENERGY - Used to compute goodness of fit measure
%
% Compute the energy covered by an arbitrary curve (weighted by a 7-point 
% Gaussian window) that maximizes the Fourier spectrum energy at each 
% frequency. This curve is not constrained by any model and may not follow
% a mode or even be smooth. However, this means it provides an upper bound
% on the A0 mode fitting routine and can be used to compute the goodness
% of fit measure.
%
%    GOF = cost(TI or Iso fit) / cost(global mode energy)
%
% Parameters
% ----------
% f   : [double, nf] array of frequency values in Hz
% k   : [double, nk] array of wavenumber values in 1/m
% VdB : [double, nk x nf] 2D power spectrum in dB (rows = k, cols = f)
% w   : [double] (optional) array of weights to apply to power spectrum
%
% Returns
% -------
% cost : [scalar] upper bound estimate on the A0 mode fit, used in
%        computing goodness of fit (Phi_max in Supplementary Note 7)
% kopt : [double, nf] array of wavenumber values in 1/m defining the
%        arbitrary f-k curve used to compute the global mode energy
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
% 
% ---

    if nargin < 4
        w = 1;
    end
    
    % Gaussian window
    gwin = gausswin(7);
    gwin = gwin./sum(gwin);

    % Convert VdB to normalized energy
    energy = w.*(10.^(VdB/10));
    
    % Compute energy in the f-k space for all window locations
    % We do this by convolving the energy with the gaussian window,
    % returning a same size array. The peak value of the convolution at
    % each frequency is added to the energy_sum and its position is returned.
    energy_sum = 0;
    sum_num = 0;
    kopt = zeros(size(f));
    for idx = 1:length(f)
        gconv = conv(energy(:,idx), gwin, 'same');
        [maxval, maxidx] = max(gconv);
        energy_sum = energy_sum + maxval;
        kopt(idx) = k(maxidx);
        sum_num = sum_num + 1;
    end
    cost = energy_sum/sum_num;