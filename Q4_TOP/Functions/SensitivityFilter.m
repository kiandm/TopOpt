function [W, WA] = SensitivityFilter(ecent, rmin, Ae)
% Cone weights w_ij = max(0, rmin - d_ij); WA = sum_j w_ij*A_j
    nel = size(ecent,1);
    ii = cell(nel,1); jj = cell(nel,1); ww = cell(nel,1);
    for i = 1:nel
        d = hypot(ecent(:,1)-ecent(i,1), ecent(:,2)-ecent(i,2));
        j = find(d < rmin);
        ii{i} = repmat(i, numel(j), 1);  jj{i} = j;  ww{i} = rmin - d(j);
    end
    W  = sparse(vertcat(ii{:}), vertcat(jj{:}), vertcat(ww{:}), nel, nel);
    WA = W * Ae;
end