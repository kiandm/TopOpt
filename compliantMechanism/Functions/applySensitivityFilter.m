function dcf = applySensitivityFilter(dc, Wi, Wj, Ww, Wsum)
% dĉi = [sumj wij * xj * dcj] / [xi * sumj wij]
    %nel = numel(x);
    %num = zeros(nel,1);
    %for k = 1:numel(Ww)
        %i = Wi(k); j = Wj(k); w = Ww(k);
        %num(i) = num(i) + w * x(j) * dc(j);
    %end
    %dcf = num ./ (x .* Wsum);

    % Correct adjoint of linear density filter — no x needed
    nel = length(dc);
    dcf = zeros(nel, 1);
    for k = 1:length(Wi)
        dcf(Wj(k)) = dcf(Wj(k)) + Ww(k) * dc(Wi(k)) / Wsum(Wi(k));
    end
end
