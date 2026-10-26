function runCase(opts)
% Runs main_heaviside.m for one options struct in THIS function's workspace, so that a calling loop's
% variables are not wiped by main_heaviside's "clearvars -except opts".
%
%   for p = [8 16 32]
%       runCase(struct('p', p, 'plot', false));
%   end
main_heaviside;
end