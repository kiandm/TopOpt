function [Wi, Wj, Ww, Wsum] = SensitivityFilter(ecent, rmin)
% Returns sparse weight triplets (i,j,wij) with fac = rmin - dist >= 0
% and Wsum(i) = sumj wij for normalization.
    nel = size(ecent,1);
    rij = []; cij = []; wij = [];
    for i = 1:nel
        pi = ecent(i,:);
        for j = 1:nel
            pj = ecent(j,:);
            dist = norm(pi - pj);
            fac  = rmin - dist;
            if fac > 0
                rij(end+1,1) = i; 
                cij(end+1,1) = j; 
                wij(end+1,1) = fac; 
            end
        end
    end
    Wi = rij; Wj = cij; Ww = wij;
    % sum of weights per row i
    Wsum = accumarray(Wi, Ww, [nel,1], @sum);
end
