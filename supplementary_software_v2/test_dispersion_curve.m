% Constants
rho = 1000;
rho_l = 1000;
c_l = 1480;
cp = 1540;
h = 0.00071;
G = 26000;
mu = G;
lambda = rho * cp^2 - 2 * mu;

v_l = sqrt((lambda + 2 * mu) / rho);
v_t = sqrt(mu / rho);

% Frequency and angular frequency
f = linspace(1, 7000, 1024); % Avoid f = 0 to prevent division by zero
w = 2 * pi * f;             % Angular frequency

% Initialize array to store k values
k_values = zeros(size(w));

% Equation to solve
eqn = @(k, w) ((k.^2 - (w / v_t).^2).^2 .* (w / v_l) * h + ...
               4 * k.^2 .* (w / v_l) .* (w / v_t).^2 * h);

% Solve for k at each w
for i = 1:length(w)
    omega = w(i);
    k_initial_guess = omega / v_l; % Reasonable initial guess for k
    % Numerical solver with error handling
    try
        k_solution = fzero(@(k) eqn(k, omega), k_initial_guess);
        k_values(i) = k_solution;
    catch
        k_values(i) = NaN; % Assign NaN if solver fails
        warning('Solver failed at w = %.2f rad/s', omega);
    end
end

% Filter out NaN values for plotting
valid_idx = ~isnan(k_values);
w_valid = w(valid_idx);
k_valid = k_values(valid_idx);

% Plot k as a function of w
figure;
plot(w_valid, k_valid, 'LineWidth', 1.5);
xlabel('Angular Frequency (rad/s)', 'FontSize', 12);
ylabel('Wave Number (k)', 'FontSize', 12);
title('Wave Number (k) as a Function of Angular Frequency (\omega)', 'FontSize', 14);
grid on;
