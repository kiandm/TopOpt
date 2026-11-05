function T = mergeThetaSweep(agg, folder, expected)
% mergeThetaSweep  Join the per-angle CSVs that thetaSweep.m writes when it runs as a Slurm job array.
%   T = mergeThetaSweep('separated')   or   T = mergeThetaSweep('smoothmax')
%   T = mergeThetaSweep(agg, folder, expected)
% Delete old theta_sweep_part_*.csv files before a new sweep: every part found is merged, whatever its age.
% expected = angles the sweep should contain (default 5:5:175, keep in step with th_deg in thetaSweep.m). A task that was
% killed (e.g. by the walltime) writes no CSV at all, so its angle is reported as missing here.
if nargin < 2, folder = 'Results'; end
if nargin < 3, expected = 5:5:175; end

d = dir(fullfile(folder, sprintf('theta_sweep_part_%s_th*.csv', agg)));
if isempty(d), error('mergeThetaSweep:none', 'No theta_sweep_part_%s_th*.csv in %s.', agg, folder); end
parts = cell(numel(d), 1);
for k = 1:numel(d)
    f    = fullfile(d(k).folder, d(k).name);
    io   = detectImportOptions(f);                    % a failed run is a row of empty fields: fix the types, or the parts cannot be stacked
    isTx = ismember(io.VariableNames, {'exitflag', 'file'});
    io   = setvartype(io, io.VariableNames(~isTx), 'double');
    io   = setvartype(io, io.VariableNames(isTx),  'string');
    parts{k} = readtable(f, io);
end
T = sortrows(vertcat(parts{:}), 'theta_deg');

missing = setdiff(expected, T.theta_deg);
if ~isempty(missing)
    warning('mergeThetaSweep:missing', '%s: no result for theta_0 = %s deg (task killed or not run: check the Slurm log and sacct).', agg, mat2str(missing(:).'));
end
failed = T.theta_deg(T.exitflag == "error");
if ~isempty(failed)
    warning('mergeThetaSweep:failed', '%s: run failed for theta_0 = %s deg (see the warning in that task''s Slurm log).', agg, mat2str(failed(:).'));
end

out = fullfile(folder, sprintf('theta_sweep_%s_%s.csv', agg, datestr(now, 'yyyymmdd_HHMMSS')));
writetable(T, out);
fprintf('%s: %d of %d angles merged -> %s\n', agg, height(T), numel(expected), out);
end