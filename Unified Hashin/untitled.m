clear; clc;
p = 8;  rho = 8;
MT = 1 + [-0.5435 -0.0175 -0.6408 -0.7643];
MC = 1 + [-0.6945 -0.7643 -0.6432 -0.0175];
FT = 1 + [-0.0175 -0.6408 -0.7643 -0.5435];
FC = 1 + [-0.7643 -0.6432 -0.0175 -0.6945];
smax = 1 + [-0.0175 -0.0175 -0.0175 -0.0175];
N = numel(smax);

aggP  = @(FI) sum(FI.^p)^(1/p) - 1;                                     % what the code does
aggPN = @(FI) (sum(FI.^p)/N)^(1/p) - 1;                                 % with 1/N
aggK  = @(FI) max(FI) + log(sum(exp(rho*(FI-max(FI)))))/rho - 1;        % KS, overflow-safe
aggKN = @(FI) max(FI) + log(sum(exp(rho*(FI-max(FI))))/N)/rho - 1;      % KS with 1/N

labels = {'p-norm','p-mean (1/N)','KS','KS (1/N)'};
fns    = {aggP, aggPN, aggK, aggKN};
fprintf('%-14s %8s %8s %8s %8s | %8s | %s\n', '', 'FT','FC','MT','MC', 'smax', 'smax - max(sep)');
for k = 1:4
    f = fns{k};
    gs = [f(FT) f(FC) f(MT) f(MC)];
    fprintf('%-14s %8.4f %8.4f %8.4f %8.4f | %8.4f | %8.4f\n', labels{k}, gs, f(smax), f(smax) - max(gs));
end