function [hash, dirty] = getGitInfo()
% Short commit hash of the code that is running, and whether tracked files have uncommitted changes.
% hash = 'nogit' (and dirty = NaN) when git is not available or this is not a repository.
hash = 'nogit';  dirty = NaN;
[st, out] = system('git rev-parse --short HEAD');
if st == 0
    hash = strtrim(out);
    [st2, ~] = system('git diff --quiet HEAD');
    dirty = (st2 ~= 0);
end
end