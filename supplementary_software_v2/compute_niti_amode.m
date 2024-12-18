function c = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l)
% COMPUTE_NITI_AMODE - Improved function for computing A0 mode
% 
% c = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l)
%
% This function uses an improved mode-tracing algorithm to compute the A0
% mode of a NITI material in the frequency-wavenumber range given
% by f and k. The model assumes a uniform, NITI plate of thickness
% h bounded above by air and below by a liquid. 
%
% The mode tracing algorithm uses an exhaustive search to find three
% high-frequency points on the mode which are then improved using a bounded
% minimization of abs(detM). These serve as initial points. New trial
% points along the mode are estimated using parabolic extrapolation with
% the three points above it in frequency. The trial points are then refined
% with bounded minimization. 
%
% Parameters
% ----------
% f      : [vector] frequency in Hz
% k      : [vector] wavenumber in 1/m
% h      : [scalar] plate thickness in m
% G      : [scalar] second shear modulus
% mu     : [scalar] first shear modulus
% lambda : [scalar] Lame parameter
% rho    : [scalar] density of the plate
% rho_l  : [scalar] density of the liquid
% c_l    : [scalar] acoustic wave speed of the liquid
%
% Returns
% -------
% c : [vector, size of f] phase velocity of the A0 mode
% 
% Notes/TODO
% ----------
% - To obtain an isotropic A0 mode, set G = mu.
% - This does not work with complex moduli.
% - The sampling of f and k may affect the performance of the mode tracing.
% - The wavenumber vector is not really used and can be replaced with a
%   tolerance (taking the place of dk) and a valid range. Alternatively, a
%   maximum wave speed can be assigned based on sqrt(G/rho) so that the k
%   range is not needed.
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
% 
% ---

    % Init
    opts = optimoptions('fmincon', 'StepTolerance', 1E-6, 'Display', 'off');
    nfreq = length(f);    
    kresult = zeros(size(f));
   
 
    % Start at largest frequency, perform a coarse exhaustive search over
    % the given k range, upsampled by a factor of 4. We know that there
    % will be multiple minima of abs(detM), but we only want the A0 mode,
    % which is the slowest mode, c = f/k. After finding the coarse root,
    % refine the estimate with a minimization algorithm. Repeat this so
    % that we have three points on the mode.
    dk = k(2) - k(1);
    ktest = k(1):(dk/4):k(end);
    for n = 0:2
        % Coarse exhaustive search
        [~, detM] = compute_niti_kappa(f(nfreq - n), ktest, h, G, mu, lambda, rho, rho_l, c_l);
        [~, locs] = findpeaks(-abs(detM));
        cpeaks = f(end)./ktest(locs);
        kinit = f(end)/min(cpeaks);
        
        % Refine estimate
        objfun = @(ki) abs_detM(f(nfreq - n), ki, h, G, mu, lambda, rho, rho_l, c_l);
        kopt = fminbnd(objfun, max(kinit - 2*dk, 0), kinit + 5*dk);

        
        % Store result
        kresult(nfreq - n) = kopt;
    end
    
    
    % Work our way down along the frequencies using only the
    % optimization-based search. To find an initial guess, use the three
    % previous points to extrapolate a parabola.
    for i = (nfreq - 3):(-1):1
        if f(i) > 0 % If f = 0, then we must have k = 0
            % Assume the zero falls on a parabolic arc fit through the next
            % three points. Use Lagrange interpolation for this, as the
            % Vandermonde approach performs poorly. Since the sampling is
            % regular in f, the terms in the interpolant are simplified
            kinit = 3*kresult(i+1) - 3*kresult(i+2) + kresult(i+3);

            % If kinit is negative, it is out of bounds. Set it to a small
            % positive number. If kinit is larger than the previous k, set
            % it to the previous k
            kinit = max(kinit, 1E-6);
            kinit = min(kinit, kresult(i+1));

            % Refine estimate with optimization. We know that the A0 mode
            % increases monotonically, so limit our search bounds
            % accordingly. We keep the search limited to a local range to
            % ensure that we do not stray too far from our initial point
            % and hop modes. Set the lower bound so that k must remain
            % positive. 
            objfun = @(ki) abs_detM(f(i), ki, h, G, mu, lambda, rho, rho_l, c_l);
            finit = objfun(kinit);
            [kopt, fopt] = fminbnd(objfun, max(kinit - 2*dk, 0), kinit + 5*dk);
            
            % Check that we have actually improved the result. If not,
            % return the initial guess.
            if fopt > finit
                kopt = kinit;
            end

            % Store result
            kresult(i) = kopt;
        end
    end
    
    c = f./kresult;
    c(kresult == 0) = 0;
    
    
function v = abs_detM(f, k, h, G, mu, lambda, rho, rho_l, c_l)
    [~, detM] = compute_niti_kappa(f, k, h, G, mu, lambda, rho, rho_l, c_l);
    v = abs(detM);
    