
function S = plotThetaSweep(T, opt)
% plotThetaSweep  Diagnose the local-optima problem from a sweep over the initial fibre angle theta_0.
%
%   plotThetaSweep(T)                          T = table of the sweep, one row per run (or a path to a CSV of it)
%   plotThetaSweep('Results/theta_sweep.csv')
%   S = plotThetaSweep(T, c_max=250, tol_g=1e-2, tol_plateau=1e-3, tol_M=inf, theta_ref=90, save=true)
%
% Needed columns: theta_deg, V, C0, M, iter, dV_win, exitflag, g_ft g_fc g_mt g_mc g_c, maxFI_ft maxFI_fc maxFI_mt maxFI_mc.
% C0 is read as the compliance of the starting design (same units as c_max).
%
% Every run is classified by re-applying the stopping rule of main_heaviside.m to the table (feasible = max g <= tol_g,
% plateau = dV_win < tol_plateau, greyness M <= tol_M), so each run that hit maxiter is explained rather than lumped
% together. The defaults are those of Functions/setOpts.m. A warning is raised if the re-computed rule disagrees with
% the exitflag in the table (wrong tolerances).
%
%   marker              meaning (same in every panel)
%   blue filled circle  converged: a genuine local optimum, V can be compared
%   hollow markers      hit maxiter: V is NOT a converged optimum. vermillion = constraint still violated (max g > tol_g),
%                       green = volume still moving (dV_win >= tol_plateau), purple = both
%
% Figure 1  Outcome vs start angle   V, C0 and greyness M against theta_0 (stacked, shared axis). Does V depend on theta_0,
%                                    and does it follow how good the starting design is?
% Figure 2  Can the runs be trusted? iterations, plateau measure, constraint values and max failure indices against theta_0.
% Figure 3  What decides the result? V against C0, V against the worst constraint value, and the runs ranked by V
%                                    (plateaus in the ranking = distinct local optima; a smooth ramp = a continuum).
%
% S (optional output): S.T is T plus gmax and status columns, S.best the row of the best converged run, and the rank
% correlations S.rho_C0_all / S.rho_C0_conv / S.rho_g_all.
arguments
    T
    opt.c_max       (1,1) double  = 250     % compliance limit [N mm]
    opt.tol_g       (1,1) double  = 1e-2    % opts.tol_g in setOpts.m
    opt.tol_plateau (1,1) double  = 1e-3    % opts.tol_plateau in setOpts.m
    opt.tol_M       (1,1) double  = inf     % opts.tol_M in setOpts.m (inf = greyness is not part of the stopping rule)
    opt.maxiter     (1,1) double  = NaN     % NaN: taken as max(T.iter)
    opt.theta_ref   (1,1) double  = 90      % default starting angle [deg] (opts.theta_init = pi/2)
    opt.save        (1,1) logical = false   % write the figures (pdf + png) to opt.folder
    opt.folder      (1,:) char    = 'figures'
end

if ischar(T) || isstring(T), T = readtable(T); end
need    = {'theta_deg','V','C0','M','iter','dV_win','exitflag','g_ft','g_fc','g_mt','g_mc','g_c', ...
           'maxFI_ft','maxFI_fc','maxFI_mt','maxFI_mc'};
missing = need(~ismember(need, T.Properties.VariableNames));
if ~isempty(missing), error('plotThetaSweep:columns', 'Table is missing column(s): %s', strjoin(missing, ', ')); end

T  = sortrows(T, 'theta_deg');
n  = height(T);
th = T.theta_deg;
if isnan(opt.maxiter), opt.maxiter = max(T.iter); end
floorVal = 1e-5;                                         % constraint values <= this are drawn on the floor of the log axes

%% Classify every run with the stopping rule of main_heaviside.m
gHashin  = max([T.g_ft T.g_fc T.g_mt T.g_mc], [], 2);    % worst Hashin mode constraint
gmax     = max([gHashin T.g_c], [], 2);                  % what the 'feasible' test sees
isConv   = lower(string(T.exitflag)) == "converged";
feasible = gmax <= opt.tol_g;
plateau  = T.dV_win < opt.tol_plateau;
greyOK   = T.M <= opt.tol_M;

status = ones(n, 1);                                     % 1 converged
status(~isConv & ~feasible &  plateau) = 2;              % maxiter, constraint still violated
status(~isConv &  feasible & ~plateau) = 3;              % maxiter, volume still moving
status(~isConv & ~feasible & ~plateau) = 4;              % maxiter, both
status(~isConv &  feasible &  plateau) = 5;              % maxiter, neither (greyness / unexplained)

bad = (feasible & plateau & greyOK) ~= isConv;
if any(bad)
    warning('plotThetaSweep:rule', ['The stopping rule re-computed from the table disagrees with exitflag for ' ...
            'theta_0 = %s deg. Check tol_g / tol_plateau / tol_M.'], mat2str(th(bad).'));
end

names  = {'converged', 'maxiter: constraint violated', 'maxiter: volume still moving', ...
          'maxiter: violated + moving', 'maxiter: other'};
colors = [0.000 0.447 0.698      % blue         (Okabe-Ito, colour-blind safe)
          0.835 0.369 0.000      % vermillion
          0.000 0.620 0.451      % bluish green
          0.800 0.475 0.655      % reddish purple
          0.500 0.500 0.500];    % grey
marks   = {'o', 'd', 's', '^', 'v'};
present = ismember(1:5, status);
grey    = [0.80 0.80 0.80];

%% Reference run: best converged one (an unconverged V is not a local optimum)
if any(isConv)
    idxC = find(isConv);
    [Vbest, k] = min(T.V(idxC));
    ibest = idxC(k);
else
    [Vbest, ibest] = min(T.V);
    warning('plotThetaSweep:noneConverged', 'No run converged: the reference is the best unconverged run.');
end
Vconv  = T.V(isConv);
if isempty(Vconv), Vconv = T.V; end
Vworst = max(Vconv);
rho_all  = rankCorr(T.V, T.C0);
rho_conv = rankCorr(T.V(isConv), T.C0(isConv));
rho_g    = rankCorr(T.V, gmax);
nWithin5 = nnz(T.V <= 1.05 * Vbest);

% The V axes are zoomed to +20% of the best run, otherwise the few very bad runs squash the interesting range. Worse
% runs are drawn on the top edge and listed in a note under the panel.
Vtop  = 1.2 * Vbest;
Vclip = min(T.V, Vtop);
over  = T.V > Vtop;
if any(over)
    note = sprintf('off scale, drawn on the edge: \\theta_0 = %s  (V up to %.3f)', strjoin(string(th(over)).', ', '), max(T.V));
else
    note = '';
end
Vlim = [0.98 * min(T.V), 1.01 * Vtop];
%% Figure 1: final design against start angle ------------------------------------------------------------------------
f1 = figure('Name', 'Theta sweep 1: outcome vs start angle', 'Color', 'w', 'Position', [100 80 900 900]);
tl = tiledlayout(f1, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Final design against initial fibre angle  (%d runs, dotted line = default start %g\\circ)', n, opt.theta_ref));

ax1 = nexttile(tl);  setupAxes(ax1);
plot(ax1, th, T.V, '-', 'Color', grey, 'HandleVisibility', 'off');
h = byStatus(ax1, th, Vclip, status, colors, marks, names);
thetaAxis(ax1, opt.theta_ref);  ylim(ax1, Vlim);
yline(ax1, Vbest, '--', sprintf('best converged, V = %.4f', Vbest));
title(ax1, '(a) final volume fraction V (objective)');
subtitle(ax1, withNote(sprintf('converged runs span V = %.4f to %.4f (+%.0f%%);  hollow markers are not converged optima', ...
                               Vbest, Vworst, 100 * (Vworst / Vbest - 1)), note));
ylabel(ax1, 'V');  ax1.XTickLabel = {};
lgd = legend(ax1, h(present), names(present), 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'north';

ax2 = nexttile(tl);  setupAxes(ax2);
plot(ax2, th, T.C0, '-', 'Color', grey, 'HandleVisibility', 'off');
byStatus(ax2, th, T.C0, status, colors, marks, names);
thetaAxis(ax2, opt.theta_ref);
yline(ax2, opt.c_max, '--', sprintf('c_{max} = %g', opt.c_max));
title(ax2, '(b) compliance of the starting design C_0');
subtitle(ax2, 'below the dashed line the start is already compliance-feasible;  compare the dip with the dip in (a)');
ylabel(ax2, 'C_0');  ax2.XTickLabel = {};

ax3 = nexttile(tl);  setupAxes(ax3);
plot(ax3, th, T.M, '-', 'Color', grey, 'HandleVisibility', 'off');
byStatus(ax3, th, T.M, status, colors, marks, names);
thetaAxis(ax3, opt.theta_ref);
if isfinite(opt.tol_M), yline(ax3, opt.tol_M, '--', 'tol_M'); end
title(ax3, '(c) greyness M of the final design [%]');
subtitle(ax3, 'a high M means the design is far from 0/1, so its V is less comparable with the others');
ylabel(ax3, 'M [%]');  xlabel(ax3, '\theta_0 [deg]   (fibre angle repeats every 180\circ, so 175\circ sits next to 5\circ)');
%% Figure 2: can the runs be trusted as local optima? ----------------------------------------------------------------
f2 = figure('Name', 'Theta sweep 2: exit diagnostics', 'Color', 'w', 'Position', [140 100 1250 780]);
tl = tiledlayout(f2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'Why did each run stop, and is it feasible?');

ax = nexttile(tl);  setupAxes(ax);
h = byStatus(ax, th, T.iter, status, colors, marks, names);
thetaAxis(ax, opt.theta_ref);  ylim(ax, [0 1.1 * opt.maxiter]);
yline(ax, opt.maxiter, ':', 'maxiter');
title(ax, '(a) iterations used');  ylabel(ax, 'iterations');  xlabel(ax, '\theta_0 [deg]');
lgd = legend(ax, h(present), names(present), 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'north';

ax = nexttile(tl);  setupAxes(ax);  ax.YScale = 'log';
byStatus(ax, th, T.dV_win, status, colors, marks, names);
thetaAxis(ax, opt.theta_ref);
yline(ax, opt.tol_plateau, '--', sprintf('tol plateau = %g', opt.tol_plateau));
title(ax, '(b) relative volume change over the last window');
subtitle(ax, 'must be below the line to stop');  ylabel(ax, 'dV_{win}');  xlabel(ax, '\theta_0 [deg]');

ax = nexttile(tl);  setupAxes(ax);  ax.YScale = 'log';
line(ax, th, max(T.g_c, floorVal), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 12, 'Color', [0.55 0.55 0.55]);
byStatus(ax, th, max(gHashin, floorVal), status, colors, marks, names);
thetaAxis(ax, opt.theta_ref);
yline(ax, opt.tol_g, '--', sprintf('tol g = %g', opt.tol_g));
title(ax, '(c) constraint values at exit  (g <= 0 feasible)');
subtitle(ax, sprintf('coloured = worst Hashin mode, grey dots = compliance g_c;  values <= %g drawn at %g', floorVal, floorVal));
ylabel(ax, 'g');  xlabel(ax, '\theta_0 [deg]');

ax = nexttile(tl);  setupAxes(ax);
FI    = [T.maxFI_ft T.maxFI_fc T.maxFI_mt T.maxFI_mc];
modes = {'fibre tension', 'fibre compression', 'matrix tension', 'matrix compression'};
hp = plot(ax, th, FI, '-o', 'MarkerSize', 4, 'LineWidth', 1.2);
set(hp, {'DisplayName'}, modes(:));
thetaAxis(ax, opt.theta_ref);
yline(ax, 1, '--', 'failure (FI = 1)', 'HandleVisibility', 'off');
legend(ax, hp, modes, 'Location', 'best', 'Box', 'off');
title(ax, '(d) max failure index per Hashin mode');
subtitle(ax, 'which mode governs, and which runs actually exceed FI = 1');
ylabel(ax, 'max FI');  xlabel(ax, '\theta_0 [deg]');
%% Figure 3: what decides the result? --------------------------------------------------------------------------------
f3 = figure('Name', 'Theta sweep 3: what decides the optimum', 'Color', 'w', 'Position', [180 60 1500 800]);
tl = tiledlayout(f3, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'What decides the final volume?  (labels = \theta_0 in degrees)');

ax = nexttile(tl);  setupAxes(ax);
h = byStatus(ax, T.C0, Vclip, status, colors, marks, names);
text(ax, T.C0(~over), Vclip(~over), "  " + string(th(~over)), 'FontSize', 7, 'Color', [0.35 0.35 0.35], 'Clipping', 'on');
xline(ax, opt.c_max, '--', 'C_0 = c_{max}');  ylim(ax, Vlim);
title(ax, '(a) V against the starting design''s compliance');
subtitle(ax, sprintf('Spearman \\rho = %.2f (all runs), %.2f (converged only)', rho_all, rho_conv));
xlabel(ax, withNote('C_0', note));  ylabel(ax, 'V');
lgd = legend(ax, h(present), names(present), 'Orientation', 'horizontal', 'Box', 'off');
lgd.Layout.Tile = 'north';

ax = nexttile(tl);  setupAxes(ax);  ax.XScale = 'log';
gx = max(gmax, floorVal);
byStatus(ax, gx, Vclip, status, colors, marks, names);
text(ax, gx(~over), Vclip(~over), "  " + string(th(~over)), 'FontSize', 7, 'Color', [0.35 0.35 0.35], 'Clipping', 'on');
xline(ax, opt.tol_g, '--', sprintf('tol g = %g', opt.tol_g));  ylim(ax, Vlim);
title(ax, '(b) V against the worst constraint value');
subtitle(ax, sprintf('Spearman \\rho = %.2f: does a lower V come from violating a constraint?', rho_g));
xlabel(ax, withNote('max g at exit', note));  ylabel(ax, 'V');

ax = nexttile(tl);  setupAxes(ax);
[Vs, order] = sort(T.V);
Vsc = min(Vs, Vtop);
x0  = 0.98 * min(Vs);
plot(ax, [x0 * ones(1, n); Vsc(:).'], [1:n; 1:n], '-', 'Color', grey, 'HandleVisibility', 'off');
byStatus(ax, Vsc, (1:n).', status(order), colors, marks, names);
xline(ax, 1.05 * Vbest, ':', '+5% of best');
set(ax, 'YDir', 'reverse');  ylim(ax, [0.5 n + 0.5]);  xlim(ax, [x0 1.01 * Vtop]);
yticks(ax, 1:n);  yticklabels(ax, string(th(order)));  ax.YAxis.FontSize = 8;
title(ax, sprintf('(c) runs ranked by V:  %d of %d within 5%% of the best', nWithin5, n));
subtitle(ax, 'flat steps = distinct local optima;  a smooth ramp = V varies continuously with \theta_0');
xlabel(ax, withNote('V', note));  ylabel(ax, '\theta_0 [deg]   (best at top)');
%% Console summary ---------------------------------------------------------------------------------------------------
nConv = nnz(isConv);
fprintf('\n=== theta_0 sweep: %d runs, theta_0 = %g ... %g deg ===\n', n, th(1), th(end));
fprintf('Exit: %d converged, %d maxiter\n', nConv, n - nConv);
for k = 2:5
    if any(status == k)
        fprintf('      %2d x %-30s theta_0 = %s\n', nnz(status == k), names{k}, mat2str(th(status == k).'));
    end
end
fprintf('Best converged run: theta_0 = %g deg, V = %.4f, C0 = %.1f\n', th(ibest), Vbest, T.C0(ibest));
i0 = find(abs(th - opt.theta_ref) < 1e-9, 1);
if ~isempty(i0)
    fprintf('Default start (%g deg): V = %.4f (%+.1f%% vs best), %s\n', opt.theta_ref, T.V(i0), ...
            100 * (T.V(i0) / Vbest - 1), char(names{status(i0)}));
end
if nConv > 0
    fprintf('Converged V: min %.4f, max %.4f (+%.0f%% over best), %d of %d converged runs within 5%% of best\n', ...
            Vbest, Vworst, 100 * (Vworst / Vbest - 1), nnz(isConv & T.V <= 1.05 * Vbest), nConv);
end
fprintf('Spearman rho(V, C0): %.2f all runs, %.2f converged only;  rho(V, max g): %.2f all runs\n', rho_all, rho_conv, rho_g);

%% Output and export -------------------------------------------------------------------------------------------------
if opt.save
    if ~isfolder(opt.folder), mkdir(opt.folder); end
    figs  = [f1 f2 f3];
    stems = {'thetaSweep_1_outcome', 'thetaSweep_2_diagnostics', 'thetaSweep_3_drivers'};
    for k = 1:3
        exportgraphics(figs(k), fullfile(opt.folder, [stems{k} '.pdf']), 'ContentType', 'vector');
        exportgraphics(figs(k), fullfile(opt.folder, [stems{k} '.png']), 'Resolution', 200);
    end
    fprintf('Figures written to %s\n', opt.folder);
end
if nargout
    T.gmax   = gmax;
    T.status = categorical(names(status).', names);
    S = struct('T', T, 'best', ibest, 'rho_C0_all', rho_all, 'rho_C0_conv', rho_conv, 'rho_g_all', rho_g);
end
end

%% Local functions ----------------------------------------------------------------------------------------------------
function setupAxes(ax)
hold(ax, 'on');  grid(ax, 'on');  box(ax, 'on');
ax.TitleHorizontalAlignment = 'left';
end

function lab = withNote(txt, note)
% Axis label / subtitle with an optional second line
lab = {txt};
if ~isempty(note), lab{end+1} = note; end
end

function thetaAxis(ax, theta_ref)
% 0..180 deg axis (the fibre angle repeats every 180 deg) with the default start marked
xlim(ax, [0 180]);  xticks(ax, 0:30:180);
xline(ax, theta_ref, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
end
function h = byStatus(ax, x, y, status, colors, marks, names)
% One marker-only line per status, so the legend and the colours/shapes are identical in every panel.
% Converged = filled marker, every maxiter reason = hollow marker.
h = gobjects(1, numel(names));
for k = 1:numel(names)
    sel = status == k;
    if k == 1, face = colors(k, :); else, face = ax.Color; end
    xs = x(sel);  ys = y(sel);
    if isempty(xs), xs = NaN;  ys = NaN; end          % line() returns no object for empty data, which breaks the legend handles
    %h(k) = line(ax, x(sel), y(sel), 'LineStyle', 'none', 'Marker', marks{k}, 'MarkerSize', 7, 'LineWidth', 1.3, ...
    h(k) = line(ax, xs, ys, 'LineStyle', 'none', 'Marker', marks{k}, 'MarkerSize', 7, 'LineWidth', 1.3, ...
                'Color', colors(k, :), 'MarkerFaceColor', face, 'DisplayName', names{k});
end
end

function r = rankCorr(x, y)
% Spearman rank correlation without the Statistics Toolbox (ties are broken by order, adequate for continuous data)
if numel(x) < 3, r = NaN; return; end
c = corrcoef(ranks(x), ranks(y));
r = c(1, 2);
end

function r = ranks(x)
[~, i] = sort(x(:));
r = zeros(numel(x), 1);
r(i) = (1:numel(x)).';
end