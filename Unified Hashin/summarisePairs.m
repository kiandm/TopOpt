function S = summarisePairs(ndiv, aggs, doPlot)
% SUMMARISEPAIRS  One row per load from the files scenarioPairs.m saved: Results/scenarioPairs_<agg>_s<load>_ndiv<ndiv>.mat
%   summarisePairs(50)         also after the Slurm tasks have finished (copy Results/ back first)
% own g  = worst of the compliance constraint and the failure constraint(s) the run optimised; feasible if <= 0.02
% gpAll / gpP99 = worst / 99th-percentile Gauss-point index of the final design, REAL strengths (> 1 means a Gauss point fails)
% active = modes whose MMA shadow price exceeds 1e-3 (trust it only in a converged or nearly converged run)
% norm ratio = ||smooth-max||_P / max_k ||mode k||_P on the SEPARATED design: how much stricter the single norm is on that design
if nargin < 2 || isempty(aggs), aggs = {'separated', 'smoothmax'}; end
if nargin < 3, doPlot = true; end
tolfeas = 0.02;   modes = {'ft','fc','mt','mc'};
loads = [];
for ia = 1:2
    d = dir(fullfile('Results', sprintf('scenarioPairs_%s_s*_ndiv%d.mat', aggs{ia}, ndiv)));
    for k = 1:numel(d)
        t = regexp(d(k).name, '_s([0-9.]+)_ndiv', 'tokens', 'once');
        loads(end+1) = str2double(t{1});                                    %#ok<AGROW>
    end
end
loads = unique(loads);
if isempty(loads), error('summarisePairs:none', 'No scenarioPairs files for ndiv = %d in Results/.', ndiv); end
fprintf('ndiv = %d, feasible = own g <= %.2f;  A = %s, B = %s\n', ndiv, tolfeas, aggs{1}, aggs{2});
fprintf('%-5s | %-9s %7s %6s %6s %6s %-13s | %-9s %7s %6s %6s %6s %-13s | %7s | %s\n', 'load', 'A:exit', 'V', 'own g', 'gpP99', 'gpAll', 'active', ...
        'B:exit', 'V', 'own g', 'gpP99', 'gpAll', 'active', 'dV[%]', 'norm ratio');
V = nan(numel(loads), 2);  feas = false(numel(loads), 2);
for il = 1:numel(loads)
    txt = cell(1, 2);  ratio = NaN;
    for ia = 1:2
        f = fullfile('Results', sprintf('scenarioPairs_%s_s%.2f_ndiv%d.mat', aggs{ia}, loads(il), ndiv));
        if ~exist(f, 'file'), txt{ia} = sprintf('%-9s %7s %6s %6s %6s %-13s', 'missing', '-', '-', '-', '-', '');  continue; end
        F = load(f);  r = F.r;
        own = max([r.g_c, r.g_hs]);
        V(il, ia) = r.V;  feas(il, ia) = own <= tolfeas;
        if strcmp(r.opts.agg, 'separated')
            act = strjoin(modes(r.lam(2:5) > 1e-3), '+');
            P = r.opts.p;  solid = r.xphy(1:numel(r.xphy)/2) > 0.3;
            n4 = sum(r.FI4(solid,:).^P, 1).^(1/P);  ratio = (sum(r.FI1(solid).^P)^(1/P)) / max(n4);
        else
            act = '';  if r.lam(2) > 1e-3, act = 'smax'; end
        end
        txt{ia} = sprintf('%-9s %7.4f %6.3f %6.2f %6.2f %-13s', r.exitflag, r.V, own, r.gp.gpP99, r.gp.gpMaxAll, act);
    end
    fprintf('%-5.2f | %s | %s | %7.1f | %.3f\n', loads(il), txt{1}, txt{2}, 100 * (V(il,2) - V(il,1)) / V(il,1), ratio);
end
S = struct('loads', loads, 'V', V, 'feasible', feas);
both = feas(:,1) & feas(:,2);
if any(both)
    fprintf('both feasible in %d of %d loads; mean dV = %.1f%% over those (dV = (V_B - V_A)/V_A)\n', nnz(both), numel(loads), mean(100 * (V(both,2) - V(both,1)) ./ V(both,1)));
else
    fprintf('no load where both runs are feasible\n');
end
if doPlot
    figure;  hold on;  grid on
    plot(loads, V(:,1), '-o', 'DisplayName', aggs{1});  plot(loads, V(:,2), '-s', 'DisplayName', aggs{2});
    for ia = 1:2, plot(loads(~feas(:,ia)), V(~feas(:,ia), ia), 'kx', 'MarkerSize', 12, 'HandleVisibility', 'off'); end
    xlabel('load\_scale');  ylabel('volume fraction V');  legend('Location', 'northwest');
    title(sprintf('four-mode cantilever, ndiv = %d  (x = own constraint violated by > %.2f)', ndiv, tolfeas));
end
end