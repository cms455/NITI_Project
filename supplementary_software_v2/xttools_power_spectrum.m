function [VdB, f, k] = xttools_power_spectrum(v, dt, dx, nfft)
% XTTOOLS_POWER_SPECTRUM - Compute 2D Fourier power spectrum for OCE
%
% Computes the 2D Fourier power spectrum in decibels given the OCE-measured
% surface vertical vibration velocity field v sampled with temporal spacing
% dt and lateral spacing dx. The FFT is computed using the nearest power of 
% 2 to the input variable nfft.
%
% Parameters
% ----------
% v    : [nx, nt] surface vertical vibration velocity field
% dt   : [scalar] temporal sampling interval in seconds
% dx   : [scalar] spatial sampling interval in meters
% nfft : [scalar] requested number of points for FFT
%
% Returns
% -------
% VdB : [nfft/2, nfft/2] 2D Fourier power spectrum in decibels
% f   : [1, nfft/2] vector of frequency bins in Hz
% k   : [1, nfft/2] vector of wavenumber bins in 1/m
%
% Notes
% -----
% - The FFT is actually computed using the nearest power of 2 above the 
%   input variable nfft.
% - This function assumes that the wave is propagating to the right as a
%   convention. To track leftward-moving shear waves, pass in flipud(v).
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
% 
% ---

    % Construct frequency/wavenumber vectors
    nfft = 2*ceil(nfft/2);
    fs = 1/dt;
    ks = 1/dx;
    f = (0:(nfft-1))/nfft*fs - fs/2;
    k = (0:(nfft-1))/nfft*ks - ks/2;
    
    % 2D FFT and conversion to dB
    vspect = fft2(v, nfft, nfft);
    VdB = 20*log10(abs(vspect)/max(abs(vspect(:))));
    clear vspect
    
    % Extract the right-moving wave components and rearrange
    VdB = VdB(1:(nfft/2), (nfft/2+1):nfft);
    VdB = fliplr(VdB);
    
    % Extract only the positive frequency-wavenumber components
    f = fftshift(f);
    k = fftshift(k);
    f = f(1:nfft/2);
    k = k(1:nfft/2);