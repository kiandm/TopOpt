function dcf = applySensitivityFilter(dc, x, Wi, Wj, Ww, Wsum)
% dĉi = [sumj wij * xj * dcj] / [xi * sumj wij]
    nel = numel(x);
    num = zeros(nel,1);
    for k = 1:numel(Ww)
        i = Wi(k); j = Wj(k); w = Ww(k);
        num(i) = num(i) + w * x(j) * dc(j);
    end
    dcf = num ./ (x .* Wsum);
end
