function r = runCase(opts)
% Runs main_heaviside.m for one options struct in THIS function's workspace and returns what a sweep table needs.
main_heaviside;
solid = xphy(1:numele) > 0.3;                              % same mask as the failure-index plot
w     = min(opts.win, iter);
vlast = iterationHistory(iter-w+1:iter, 3);
r = struct('exitflag', exitflag, 'iter', iter, 'runtime', toc(tStart), ...
           'c', c, 'g_c', c/c_max - 1, 'V', v + 1, 'M', M, ...
           'g_hs', g_hs(:)', 'maxFI', max(FailIdx(solid,:), [], 1), ...    % 1x4 [ft fc mt mc] for 'separated'
           'dV_win', (max(vlast) - min(vlast)) / vlast(end), ...
           'C0', iterationHistory(1,2));
end