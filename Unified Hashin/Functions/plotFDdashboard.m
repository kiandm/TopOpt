function plotFDdashboard(fd, xc, yc, k, mode)
% Four views of one finite-difference check.  fd = fdParse(file);  xc, yc = element centroids (x_cen, y_cen in main_heaviside.m).
% k = which section feeds panels 2-4 (default 1); mode = which mode label for those panels (default: all rows).
if nargin < 4, k = 1; end
if nargin < 5, mode = []; end
col = lines(numel(fd));
short = @(s) strtok(s, ':');                               % 'dgtw/dx: Hashin sens...' -> 'dgtw/dx'

% rows of section k that panels 2-4 use
sel = true(size(fd(k).el));
if ~isempty(mode), sel = strcmp(fd(k).mode, mode); end
if numel(unique(fd(k).el(sel))) < nnz(sel), error('plotFDdashboard: section has several modes per element, pass mode (e.g. ''mt'')'); end
el = fd(k).el(sel);  ana = fd(k).ana(sel);  num = fd(k).num(sel);  err = fd(k).abs_err(sel);  scale = fd(k).scale(sel);
err_n = err ./ scale;                                      % error as a fraction of the largest analytic value
flagged = err_n > 1e-2;                                    % >1 % of the scale

% 1. exceedance curve: fraction of entries whose error is >= x, every section on one axis
subplot(2,2,1); hold on
for s = 1:numel(fd)
    e_s = sort(fd(s).abs_err ./ fd(s).scale);
    stairs(max(e_s, 1e-12), (numel(e_s):-1:1)' / numel(e_s), 'Color', col(s,:), 'LineWidth', 1.4);
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on
xline(1e-2, ':'); xline(1e-3, ':'); xlim([1e-10 1e-1]);
xlabel('|error| / max|analytic|'); ylabel('fraction of entries with error \geq x');
legend(cellfun(short, {fd.name}, 'UniformOutput', false), 'Location', 'southwest');
title('exceedance curve (dotted: 1 % and 0.1 %)');

% 2. error against magnitude: points on the dotted lines have that RELATIVE error; a flat floor is FD/round-off noise
subplot(2,2,2);
loglog(abs(ana), max(err, 1e-16), '.', 'Color', [0.65 0.65 0.65]); hold on
loglog(abs(ana(flagged)), err(flagged), 'o', 'MarkerFaceColor', col(k,:), 'MarkerEdgeColor', 'w');
x_ref = [1e-12 1];
for r = [1 1e-2 1e-4], loglog(x_ref, r*x_ref, 'k:'); end   % 100 %, 1 %, 0.01 % relative error
xlim(x_ref); grid on
xlabel('|analytic|'); ylabel('|numeric - analytic|'); title('error vs magnitude');

% 3. parity
subplot(2,2,3);
plot(ana, num, '.', 'Color', [0.65 0.65 0.65]); hold on
plot(ana(flagged), num(flagged), 'o', 'MarkerFaceColor', col(k,:), 'MarkerEdgeColor', 'w');
lim = [-1 1] * 1.05 * max(scale);
plot(lim, lim, 'k-'); axis([lim lim]); axis square; grid on
xlabel('analytic'); ylabel('finite difference'); title('parity');

% 4. signed error at the checked entries, symmetric-log colour scale; light dots = every element (the domain footprint)
subplot(2,2,4);
plot(xc, yc, '.', 'Color', [0.88 0.88 0.87], 'MarkerSize', 4); hold on
lin = 1e-3;   f = @(v) sign(v) .* log10(1 + abs(v) / lin);   % linear inside +-0.1 %, logarithmic outside
v = sign(num - ana) .* err_n;                                % signed, normalised error
scatter(xc(el), yc(el), 36, f(v), 's', 'filled', 'MarkerEdgeColor', [0.4 0.4 0.4], 'LineWidth', 0.4);
axis equal tight
cmap = interp1([0 0.5 1], [0.16 0.47 0.84; 0.94 0.94 0.93; 0.89 0.29 0.28], linspace(0, 1, 256));   % blue - grey - red
colormap(gca, cmap); caxis([-1 1] * f(0.06));
ticks = [-0.05 -0.01 0 0.01 0.05];
colorbar('Ticks', f(ticks), 'TickLabels', {'-5 %', '-1 %', '0', '1 %', '5 %'});
title('signed error / max|analytic|');
ttl = short(fd(k).name);   if ~isempty(mode), ttl = [ttl ' (' mode ')']; end
if exist('sgtitle'), sgtitle(ttl); end
end