num_samples = length(kfit_holder);
mu_factor = 110/3;
figure;
hold on;

% Generate a colormap with enough unique colors for all samples
colors = lines(num_samples); % Use 'lines' colormap for distinct colors

for s = 1:num_samples
    kfit_final = kfit_holder{s};
    VdB = VdB_samples(s, :);

    % Choose color for this sample
    color = colors(s, :);

    % Plot kfit_final with unique color
    plot(f, kfit_final, '-', 'Color', color, 'LineWidth', 2, ...
        'DisplayName', ['Kfit Sample ', num2str(s)]);

    % Plot VdB data points with the same color but with markers
    plot(f(f_reduced_idx), 1e-3 * VdB, 'o-', 'Color', color, ...
        'DisplayName', ['K Sample ', num2str(s)],'LineWidth',1);
end

title('Dispersion Curves for Subject 11');
xlabel('Frequency (1/s)');
ylabel('Wave-number (1/mm)');
set(gca, 'FontSize', 14);
legend('show','Location','best'); % Add a legend to distinguish samples
hold off;

figure;
hold on;


kfit_sum = 0;
for s = 1:num_samples
    kfit_final = kfit_holder{s};
    VdB = VdB_samples(s, :);

    % Choose color for this sample
    color = colors(s, :);

    % Plot kfit_final with unique color
    plot(f(1,1:end),kfit_final, '-', 'Color', color, 'LineWidth', 2, ...
        'DisplayName', ['Kfit Sample ', num2str(s)]);

    % Plot VdB data points with the same color but with markers
    %plot(f(f_reduced_idx), 1e-3 * VdB, 'o-', 'Color', color, ...
     %   'DisplayName', ['VdB Sample ', num2str(s)],'LineWidth',1);
    kfit_sum = kfit_sum + kfit_final;
end

title('Dispersion Curves for Subject 11');
xlabel('Frequency (1/s)');
ylabel('Wave-number (1/mm)');
set(gca, 'FontSize', 14);
legend('show','Location','best'); % Add a legend to distinguish samples
hold off;

kfit_avg = 1e3*kfit_sum/num_samples;

G_avg = sum(G_opt_holder)/num_samples;
mu = G_avg*mu_factor;

lambda = rho * cp^2 - 2 * mu;
cfit = compute_niti_amode(f, k, h, G_avg, mu, lambda, rho, rho_l, c_l);
kfit_G_avg = f(:)./cfit(:);


G0 = 30e3;
mu0 = mu_factor * G0;
%[G_opt, mu_opt, mu_opt_div_G, kfit_final] = fit_data_to_curve_no_constraint(rho, rho_l, c_l, cp, h, G0, mu0, mu_factor, f, k, kfit_avg(2:end), [1:length(f)-1]);

lambda = rho * cp^2 - 2 * mu_opt;
cfit = compute_niti_amode(f, k, h, G_opt, mu_opt, lambda, rho, rho_l, c_l);
kfit_avg_opt = f(:)./cfit(:);


figure;
hold on;
plot(kfit_avg*1e-3,'Linewidth',2,'DisplayName','Sample K Average');
plot(kfit_G_avg*1e-3,'Linewidth',2,'DisplayName','Kfit G Averge');
plot(kfit_avg_opt*1e-3,'Linewidth',2,'DisplayName','Sample K Average Fit');

title('Averaged Sample K for Subject 11');
xlabel('Frequency (1/s)');
ylabel('Wave-number (1/mm)');
set(gca, 'FontSize', 14);
legend('show','Location','best'); % Add a legend to distinguish samples
hold off;

disp(sprintf('Average G: %d', G_avg));
disp(sprintf('Averaged Ks fit G : %d', G_opt));
