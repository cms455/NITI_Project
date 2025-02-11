function plot_fit_data_v2(kfit_holder, G_opt_holder, confidence_interval_holder, mu_factor, subject, VdB_samples, f, k, f_reduced_idx)

num_samples = length(kfit_holder);

rho = 1000;    % Density of material (kg/m^3)
rho_l = 1000;  % Density of liquid (kg/m^3)
c_l = 1480;    % Speed of sound in liquid (m/s)
cp = 1540;     % Phase velocity (m/s)
h = 800E-6;    % Thickness (m)
mu_factor = 110 / 3;

% Generate a colormap with enough unique colors for all samples
colors = lines(num_samples); % Use 'lines' colormap for distinct colors

figure;
hold on;

for s = 1:num_samples
    kfit_final = kfit_holder{s};
    VdB = VdB_samples(s, :);

    % Retrieve the confidence interval for the current sample
    ci = confidence_interval_holder{s};
    ci_lower = ci(1); % Lower bound for G
    ci_upper = ci(2); % Upper bound for G

    % Calculate mu and lambda for confidence bounds
    G_lower = ci_lower;
    G_upper = ci_upper;
    mu_lower = G_lower * mu_factor;
    mu_upper = G_upper * mu_factor;
    lambda_lower = rho * cp^2 - 2 * mu_lower;
    lambda_upper = rho * cp^2 - 2 * mu_upper;

    % Compute the curves for the confidence intervals
    cfit_lower = compute_niti_amode(f, k, h, G_lower, mu_lower, lambda_lower, rho, rho_l, c_l);
    kfit_lower = f(:) ./ cfit_lower(:);

    cfit_upper = compute_niti_amode(f, k, h, G_upper, mu_upper, lambda_upper, rho, rho_l, c_l);
    kfit_upper = f(:) ./ cfit_upper(:);

    % Plot the fitted curve
    color = colors(s, :);
    plot(f, kfit_final, '-', 'Color', color, 'LineWidth', 2, ...
        'DisplayName', ['Kfit Sample ', num2str(s)]);

    % Plot confidence intervals as shaded regions
    fill([[1:length(kfit_lower)]; fliplr([1:kfit_upper])], [kfit_lower'; fliplr(kfit_upper')], color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
        'DisplayName', ['CI Sample ', num2str(s)]);
    
    % Plot VdB data points
    plot(f(f_reduced_idx), 1e-3 * VdB, 'o', 'Color', color, ...
        'DisplayName', ['VdB Sample ', num2str(s)], 'LineWidth', 1);
end

title(sprintf('Dispersion Curves with Confidence Intervals for %s', subject));
xlabel('Frequency (1/s)');
ylabel('Wave-number (1/mm)');
set(gca, 'FontSize', 14);
legend('show', 'Location', 'best'); % Add a legend to distinguish samples
hold off;

% Plot averaged sample data with confidence intervals
figure;
hold on;

kfit_sum = 0;
for s = 1:num_samples
    kfit_final = kfit_holder{s};
    kfit_sum = kfit_sum + kfit_final;
end

kfit_avg = kfit_sum / num_samples;
G_avg = mean(G_opt_holder);
ci_avg = mean(cell2mat(confidence_interval_holder'), 1); % Average CI bounds

mu_avg = G_avg * mu_factor;
lambda_avg = rho * cp^2 - 2 * mu_avg;

cfit_avg = compute_niti_amode(f, k, h, G_avg, mu_avg, lambda_avg, rho, rho_l, c_l);
kfit_avg_curve = f(:) ./ cfit_avg(:);

% Compute CI bounds for the average
G_lower_avg = ci_avg(1);
G_upper_avg = ci_avg(2);
mu_lower_avg = G_lower_avg * mu_factor;
mu_upper_avg = G_upper_avg * mu_factor;
lambda_lower_avg = rho * cp^2 - 2 * mu_lower_avg;
lambda_upper_avg = rho * cp^2 - 2 * mu_upper_avg;

cfit_lower_avg = compute_niti_amode(f, k, h, G_lower_avg, mu_lower_avg, lambda_lower_avg, rho, rho_l, c_l);
kfit_lower_avg = f(:) ./ cfit_lower_avg(:);

cfit_upper_avg = compute_niti_amode(f, k, h, G_upper_avg, mu_upper_avg, lambda_upper_avg, rho, rho_l, c_l);
kfit_upper_avg = f(:) ./ cfit_upper_avg(:);

% Plot the averaged curve
plot(f, kfit_avg, '-', 'LineWidth', 2, 'DisplayName', 'Average Kfit');

% Plot confidence intervals for the average
fill([f; flip(f)], [kfit_lower_avg; flip(kfit_upper_avg)], 'cyan', ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'DisplayName', 'CI Average');

plot(f, kfit_avg_curve, '--', 'LineWidth', 2, 'DisplayName', 'Fitted Average Curve');

title(sprintf('Averaged Dispersion Curves for %s', subject));
xlabel('Frequency (1/s)');
ylabel('Wave-number (1/mm)');
set(gca, 'FontSize', 14);
legend('show', 'Location', 'best');
hold off;

disp(sprintf('Average G: %.2f', G_avg));
disp(sprintf('Confidence Interval for G (Avg): [%.2f, %.2f]', ci_avg(1), ci_avg(2)));
end
