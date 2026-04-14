clear; clc;

%% Material properties
E1 = 1.0;    % stiff phase (material)
E2 = 1e-4;   % soft phase (void)
nu = 0.3;

%% Volume fractions
vf = linspace(0, 1, 500);

%% --- Hashin-Shtrikman Bounds ---
% Bulk moduli
K1 = E1 / (3*(1 - 2*nu));
K2 = E2 / (3*(1 - 2*nu));

% Shear moduli
G1 = E1 / (2*(1 + nu));
G2 = E2 / (2*(1 + nu));

% HS upper bound (soft phase as matrix)
K_HS_upper = K1 + (1 - vf) ./ (1./(K2 - K1) + vf ./ (K1 + 4/3*G1));
G_HS_upper = G1 + (1 - vf) ./ (1./(G2 - G1) + 2*vf.*(K1 + 2*G1) ./ (5*G1.*(K1 + 4/3*G1)));

% HS lower bound (stiff phase as matrix)
K_HS_lower = K2 + vf ./ (1./(K1 - K2) + (1 - vf) ./ (K2 + 4/3*G2));
G_HS_lower = G2 + vf ./ (1./(G1 - G2) + 2*(1-vf).*(K2 + 2*G2) ./ (5*G2.*(K2 + 4/3*G2)));

% Convert back to effective Young's modulus
E_HS_upper = 9*K_HS_upper.*G_HS_upper ./ (3*K_HS_upper + G_HS_upper);
E_HS_lower = 9*K_HS_lower.*G_HS_lower ./ (3*K_HS_lower + G_HS_lower);

%% --- Voigt and Reuss bounds (for reference) ---
E_voigt = vf*E1 + (1-vf)*E2;
E_reuss = 1 ./ (vf/E1 + (1-vf)/E2);

%% --- SIMP interpolation for different p values ---
p_values = [1, 2, 3];
colors = lines(length(p_values));

%% --- Plot ---
figure('Position', [100, 100, 850, 580]);
hold on;

% Shade HS feasible region
fill([vf, fliplr(vf)], [E_HS_upper, fliplr(E_HS_lower)], ...
    [0.85 0.92 1.0], 'EdgeColor', 'none', 'DisplayName', 'H-S Feasible Region');

% Shade Voigt-Reuss region outside HS
fill([vf, fliplr(vf)], [E_voigt, fliplr(E_HS_upper)], ...
    [0.95 0.95 0.95], 'EdgeColor', 'none', 'DisplayName', 'Voigt-Reuss Region');
fill([vf, fliplr(vf)], [E_HS_lower, fliplr(E_reuss)], ...
    [0.95 0.95 0.95], 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Voigt and Reuss bounds
plot(vf, E_voigt, 'k--',  'LineWidth', 1.2, 'DisplayName', 'Voigt (upper)');
plot(vf, E_reuss, 'k:',   'LineWidth', 1.2, 'DisplayName', 'Reuss (lower)');

% HS bounds
plot(vf, E_HS_upper, 'b-', 'LineWidth', 2.0, 'DisplayName', 'H-S Upper');
plot(vf, E_HS_lower, 'r-', 'LineWidth', 2.0, 'DisplayName', 'H-S Lower');

% SIMP curves
for i = 1:length(p_values)
    p = p_values(i);
    E_SIMP = E2 + (E1 - E2) * vf.^p;
    plot(vf, E_SIMP, '-', 'Color', colors(i,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('SIMP p = %d', p));
end

hold off;

xlabel('Volume Fraction \rho', 'FontSize', 13);
ylabel('Normalised Young''s Modulus E/E_1', 'FontSize', 13);
title('SIMP Interpolation vs Hashin-Shtrikman Bounds', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 10);
grid on;
box on;
xlim([0 1]);
ylim([0 1.05]);