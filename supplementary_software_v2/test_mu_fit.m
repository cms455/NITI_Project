% Define mu values explicitly (in kPa, converted to Pa)
mu_values = [ 20, 40, 80, 100, 200, 400, 800, 1000, 2000] * 1e3; % Convert kPa to Pa
num_mu = length(mu_values);

% Fix G0 value
G0 = 25.6e3; % Pa

figure;
hold on;

% Initialize legend entries and plot handles
legend_entries = cell(1, num_mu); 
plot_handles = gobjects(1, num_mu);

% Define a colormap for gradient coloring
colors = colormap(parula(num_mu)); % Generate gradient colors for all mu values

% Loop over mu values and plot
for n = 1:num_mu
    mu0 = mu_values(n); % Current mu value
    disp(mu0);
    
    % Compute cfit and kfit
    cfit = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
    kfit = f(:) ./ cfit(:);
    
    % Assign color from colormap
    color_index = round((n / num_mu) * size(colors, 1));
    plot_handles(n) = plot(f, kfit, 'Color', colors(color_index, :), 'LineWidth', 2);
    
    % Add entry to legend
    legend_entries{n} = sprintf('\\mu = %d kPa', mu0 / 1e3); % Convert back to kPa for legend
end

% Add the legend
legend(plot_handles, legend_entries, 'Location', 'best', 'FontSize', 10);

% Set axis labels and title
xlabel('Frequency [Hz]', 'FontSize', 20);
ylabel('Wavenumber [1/mm)]', 'FontSize', 20);
%title('Gradient Color Plot of kfit vs Frequency for Different \mu Values');

% Finalize the figure
set(gca, 'FontSize',14);
hold off;
