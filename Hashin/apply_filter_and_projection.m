function xphy_out = apply_filter_and_projection(xval_in, numele, H, Hs, beta, eta)
% Applies density filter + Heaviside projection to density variables,
% and direct linear filter to angle variables.
% Returns physical design variable vector xphy.

    xphy_out = xval_in;

    % Density: filter then project
    x_tilde = (H * xval_in(1:numele)) ./ Hs;
    [x_proj, ~] = heavisideProjection(x_tilde, beta, eta);
    xphy_out(1:numele) = x_proj;

    % Angle: direct linear filter
    xphy_out(numele+1:end) = (H * xval_in(numele+1:end)) ./ Hs;
end