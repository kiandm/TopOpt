function x_tilde = applyDensityFilter(x, Wi, Wj, Ww, Wsum)
    nel = numel(x);
    W = sparse(Wi, Wj, Ww, nel, nel);
    x_tilde = W * x ./ Wsum;
end