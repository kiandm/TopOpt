function dcf = applySensitivityFilter(dc, x, W, WA, Ae)
% dcf_i = A_i * sum_j w_ij x_j dc_j / ( max(x_i,1e-3) * sum_j w_ij A_j )
    dcf = Ae .* (W * (x .* dc)) ./ (max(1e-3, x) .* WA);
end