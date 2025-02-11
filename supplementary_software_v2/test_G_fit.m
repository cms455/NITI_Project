% Define G0 values explicitly
G0_values = [ 20, 40, 80, 100, 200, 400, 800, 1000, 2000] * 1e3; % Convert kPa to Pa
num_G = length(G0_values);

h = 550e-6;

figure;
hold on;

% Initialize legend entries and plot handles
legend_entries = cell(1, num_G); 
plot_handles = gobjects(1, num_G); 

% Define a colormap for gradient coloring
colors = colormap(parula(num_G)); % Generate gradient colors for all G0 values

% Loop over G0 values and plot
for n = 1:num_G
    G0 = G0_values(n); % Current G0 value
    disp(G0)
    
    % Fixed parameters
    mu0 = 4e6; 
    lambda = rho * cp^2 - 2 * mu0;
    cfit = compute_niti_amode(f, k, h, G0, mu0, lambda, rho, rho_l, c_l);
    kfit = f(:) ./ cfit(:);
    
    % Assign color from colormap
    color_index = round((n / num_G) * size(colors, 1));
    plot_handles(n) = plot(f, kfit, 'Color', colors(color_index, :), 'LineWidth', 2);
    
    % Add entry to legend
    legend_entries{n} = sprintf('G = %d kPa', G0 / 1e3); % Convert back to kPa for legend
end

% Add the legend
legend(plot_handles, legend_entries, 'Location', 'best', 'FontSize', 10);

% Set axis labels and title
xlabel('Frequency [Hz]', 'FontSize', 20);
ylabel('Wavenumber [1/mm]', 'FontSize', 20);
%title('Gradient Color Plot of kfit vs Frequency for Different G0 Values');

% Finalize the figure
set(gca,'FontSize',14);
hold off;
