function ratio = aggregationPenalty(r, Ps)
% AGGREGATIONPENALTY  How much stricter is the single smooth-max p-norm than the four mode p-norms, on THIS design?
%   r  = runCase(...) result.  ratio(i) = ||H||_P / max_k ||FI_k||_P  at P = Ps(i)  (default [2 4 8 16 32 64]).
% Deterministic: no optimiser involved.  Two effects are mixed in it:
%   union effect         <= M^(1/P) when M modes are critical in DIFFERENT elements (falls to 1 as P grows)
%   averaging-order      the smooth-max index takes the max over modes at each Gauss point BEFORE the element average, the
%                        mode-separated index averages each mode FIRST; this part does not vanish as P grows and is NOT a benefit
%                        (check it against gpEval: does the design with the lower index have the worse Gauss-point peaks?)
if nargin < 2, Ps = [2 4 8 16 32 64]; end
solid = r.xphy(1:numel(r.xphy)/2) > 0.3;
FI4 = r.FI4(solid,:);  FI1 = r.FI1(solid);
ratio = zeros(size(Ps));
for i = 1:numel(Ps)
    P = Ps(i);  n4 = sum(FI4.^P, 1).^(1/P);  n1 = sum(FI1.^P)^(1/P);
    ratio(i) = n1 / max(n4);
end
n8 = sum(FI4.^8, 1).^(1/8);
fprintf('P           : %s\n', sprintf('%7g', Ps));
fprintf('||H||/max||FI_k||: %s\n', sprintf('%7.3f', ratio));
fprintf('M^(1/P), M=%d  : %s   (M = modes with a P=8 norm > 0.9 of the largest)\n', sum(n8 > 0.9*max(n8)), sprintf('%7.3f', sum(n8 > 0.9*max(n8)).^(1./Ps)));
end