addpath('Functions');  if ~exist('Results','dir'), mkdir('Results'); end
th_deg = 5:5:175;                                          % not 0 or 180: those start exactly on the angle bounds
base   = struct('agg','smoothmax','projection','heaviside','rmin_phys',4,'ndiv',75,'plot',false,'save',true);
vars = {'theta_deg','theta_rad','g_ft','g_fc','g_mt','g_mc','maxFI_ft','maxFI_fc','maxFI_mt','maxFI_mc', ...
        'g_c','c','V','M','iter','dV_win','C0','runtime'};
n = numel(th_deg);  D = nan(n, numel(vars));  flag = repmat({'error'}, n, 1);  file = repmat({''}, n, 1);
for k = 1:n
    o = base;  o.theta_init = deg2rad(th_deg(k));
    try
        r = runCase(o);
        % D(k,:)  = [th_deg(k), o.theta_init, r.g_hs, r.maxFI, r.g_c, r.c, r.V, r.M, r.iter, r.dV_win, r.C0, r.runtime];
        g  = nan(1,4);  g(1:numel(r.g_hs))   = r.g_hs;     % smoothmax has ONE constraint and ONE index: they go in the *_ft columns, the other modes stay NaN
        FI = nan(1,4);  FI(1:numel(r.maxFI)) = r.maxFI;
        D(k,:)  = [th_deg(k), o.theta_init, g, FI, r.g_c, r.c, r.V, r.M, r.iter, r.dV_win, r.C0, r.runtime];
        flag{k} = r.exitflag;
        d = dir(fullfile('Results','opt_results_*.mat'));  [~, i] = max([d.datenum]);  file{k} = d(i).name;
    catch ME
        warning('theta_init = %g deg failed: %s', th_deg(k), ME.message);
        D(k,1:2) = [th_deg(k), o.theta_init];
    end
    fprintf('theta_init %g deg: %s\n', th_deg(k), flag{k});  % main_heaviside calls clc, so print after each run
end
T = [array2table(D, 'VariableNames', vars), table(flag, file, 'VariableNames', {'exitflag','file'})];
disp(T);  writetable(T, fullfile('Results', ['theta_sweep_' datestr(now,'yyyymmdd_HHMMSS') '.csv']));