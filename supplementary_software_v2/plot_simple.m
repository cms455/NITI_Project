


rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;
h = 0.00055;
G = 26000;
mu = 4e6;
k = linspace(0, 2000, 2048);
f = linspace(0, 7000, 1024);

lambda = rho * cp^2 - 2 * mu;
cfit = compute_niti_amode(f, k, h, G, mu, lambda, rho, rho_l, c_l);
kfit = f(:)./cfit(:);

figure; 
plot(f,kfit,'Color', [0.5,0,0.5],'LineWidth', 2,'DisplayName', 'Dispersion Curve');
xlabel('Frequency [Hz]');
ylabel('Wavenumber [1/mm]');
set(gca, 'FontSize', 14);