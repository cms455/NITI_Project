% 
%
% Run this example script to generate a series of A0 mode dispersion curves
% for NITI materials of varying degrees of anisotropy. We keep G = 20 kPa 
% constant and vary mu. 
%
% The output plots are A0 mode dispersion curves in the k-f and c-f 
% representations, similar to Figure 2a and 2c, respectively.
%
% Author: John J. Pitre, Jr.
%
% Pitre, JJ, MA Kirby, DS Li, TT Shen, RK Wang, M O'Donnell, and I Pelivanov.
%    Nearly-incompressible transverse isotropy (NITI) of cornea elasticity: 
%    model and experiments with acoustic micro-tapping OCE. Scientific Reports 
%    (2020).
%


%
% Parameters
%
% rho    : [scalar] assumed corneal tissue density 
% rho_l  : [scalar] density of the liquid bounding the cornea (assume water)
% c_l    : [scalar] P-wave speed of the liquid bounding the cornea (assume water)
% cp     : [scalar] P-wave speed of the corneal tissue
% lambda : [scalar] stiffness tensor entry that enforces incompressibility, lambda >> mu > G
% h      : [scalar] thickness of the NITI layer (cornea)
% G      : [scalar] NITI out-of-plane shear modulus (typically order 10 kPa)
% mu     : [vector] NITI in-plan shear modulus (typically order 1 MPa)
% k      : [vector] wavenumber vector (1/m)
% f      : [vector] frequency vector (Hz)
rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;
h = 0.00071;
G = 26000;
mu = [G];
k = linspace(0, 2000, 2048);
f = linspace(0, 7000, 1024);


%
% Compute A0 mode for each value of mu - this represents increasing anisotropy
% 
result = cell(1, length(mu));
for i = 1:length(mu)  
    fprintf(1, 'Computing A0 mode for mu = %6.2f kPa\n', mu(i)/1000);
    
    % NITI stiffness modulus that enforces incompressibility
    % Set so that the P-wave speed is always roughly the same.
    lambda = rho*cp^2 - 2*mu(i);  
        
    % Compute A0 mode
    c = compute_niti_amode(f, k, h, G, mu(i), lambda, rho, rho_l, c_l);
    result{i}.c = c(:);
    result{i}.k = f(:)./c(:);
    
end


%
%
%

figure
subplot(1,2,1)
hold on
for i = 1:length(mu)
    plot(f/1000, result{i}.k/1000, '-', 'LineWidth', 1)
end
xlabel('Frequency (kHz)')
ylabel('Wavenumber (1/mm)')
box off
xlim([0, 7])
ylim([0, 1.5])

subplot(1,2,2)
hold on
cols = lines(length(mu));
for i = 1:length(mu)
    plot(f/1000, result{i}.c, '-', 'LineWidth', 1)
end
plot(f/1000, sqrt((f/1000*h/2)/(2*(rho))))
xlabel('Frequency (kHz)')
ylabel('Phase Velocity (m/s)')
box off
xlim([0, 7])
ylim([0, 4.5])
legend('\mu = G', '\mu = 2G', '\mu = 5G', '\mu = 10G', 'Location', 'SouthEast')
legend('boxoff')


