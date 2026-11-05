function plotFDsweep(files, hs)
% Error against FD step size, one curve per section.  files = {'finitediff_1e-3.txt', ...}, hs = [1e-3 ...]:
% run the check once per h_fd and keep each file. Error = relative L2 norm over the checked entries.
rel_l2 = nan(numel(files), 0);
for i = 1:numel(files)
    fd = fdParse(files{i});
    for s = 1:numel(fd)
        rel_l2(i,s) = norm(fd(s).abs_err) / norm(fd(s).num);
    end
end
[hs, order] = sort(hs(:));   rel_l2 = rel_l2(order,:);
loglog(hs, rel_l2, 'o-', 'LineWidth', 1.4); hold on
loglog(hs, rel_l2(end,1) * (hs / hs(end)).^2, 'k:');       % slope-2 guide through the largest-h point of curve 1
grid on; xlabel('step size h'); ylabel('relative L2 error');
legend([cellfun(@(s) strtok(s, ':'), {fd.name}, 'UniformOutput', false), {'slope 2'}], 'Location', 'best');
title({'correct gradient: falls as h^2, then rises (round-off)', 'flat plateau: error is in the gradient'});
end