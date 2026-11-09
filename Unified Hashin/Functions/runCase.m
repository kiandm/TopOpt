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
% ---- additions for scenario studies (every field above is unchanged: thetaSweep.m still works) ----------------
os = opts;  os.agg = 'separated';   om = opts;  om.agg = 'smoothmax';
[gS, ~, ~, FI4] = Hashin_separated(U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength_opt, freedofs, dphix_ref, dphiy_ref, os);
[gM, ~, ~, FI1] = Hashin_smoothmax(U, dK, KE0, xphy, penal, numele, gs, edofMat, coords, conn, matprop, strength_opt, freedofs, dphix_ref, dphiy_ref, om);
r.gS = gS(:)';  r.gM = gM;                  % BOTH aggregations evaluated on this design, whichever one was optimised
r.FI4 = FI4;  r.FI1 = FI1;                  % element-averaged indices of the final design: [ft fc mt mc] and the smooth-max index
r.gp  = gpEval(U, xphy, numele, gs, edofMat, matprop, strength, opts, dphix_ref, dphiy_ref);   % Gauss-point check, REAL strengths
r.lam = mean(lamHist(max(1,iter-19):iter,:), 1);   % shadow prices dV/dg_k, last 20 iterations: [compliance, failure rows...]
r.xphy = xphy;  r.opts = opts;  r.history = iterationHistory(1:iter,:);
end