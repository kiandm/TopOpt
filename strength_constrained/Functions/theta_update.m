function theta_new = theta_update(theta, dtheta, step)
% Steepest descent update for fibre angle
    if nargin < 3, step = 0.1; end   % tune this
    theta_new = theta - step * dtheta;
    % Optional: clamp to [-pi/2, pi/2]
    theta_new = max(-pi/2, min(pi/2, theta_new));
end