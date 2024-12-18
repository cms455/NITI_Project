% Example G and mu values (replace with your data)
load_G = load('G_values_data.mat', 'G_values');
G_values = load_G.G_values;

load_mu = load('mu_values_data.mat', 'mu_values');
mu_values = load_mu.mu_values;


all_G = G_values(:); 
all_mu = mu_values(:); 

G_min = 1e3;
G_max = 1e7;
mu_min = 1e3;
mu_max = 1e7;


[G_grid, mu_grid] = meshgrid(logspace(log10(G_min), log10(G_max), 100), ...
                             logspace(log10(mu_min), log10(mu_max), 100));

bandwidth = 0.15; 

density = zeros(size(G_grid));
for i = 1:length(all_G)
    density = density + exp(-((log10(G_grid) - log10(all_G(i))).^2 + ...
                              (log10(mu_grid) - log10(all_mu(i))).^2) / (2 * bandwidth^2));
end


density = density / (2 * pi * bandwidth^2 * length(all_G));

% Plot density contour
figure;
hold on;
contourf(G_grid, mu_grid, density, 20, 'LineColor', 'none'); % Filled contour plot
scatter(G_values, mu_values, 10, 'filled', 'DisplayName', sprintf('Noise Level %.2f', noise_levels(n)), 'MarkerFaceColor', 'r');
ylim([mu_min,mu_max]);
xlim([G_min, G_max])
set(gca, 'XScale', 'log', 'YScale', 'log'); % Logarithmic axes
colormap parula ;
cb = colorbar; 
clim([0,1])
ylabel(cb,'Density','FontSize',16,'Rotation',270)

xlabel('G (a.u.)');
ylabel('\mu (a.u.)');
%title('Gaussian Density Contour Plot of G vs. \mu');

grid on;
