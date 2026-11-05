function fd = fdParse(filename)
% Reads a finitediff*.txt written by finitediff_main_heaviside.m into one struct per '--- name ---' section.
% Handles both row layouts:  element ana num abs_err rel_err   and   element mode ana num abs_err rel_err.
% fd(k) has: name, el, mode, ana, num, abs_err, scale (per row: max|analytic| of that mode, used to normalise errors).
% Use abs_err for errors: the ana/num columns are printed with only 6 significant digits, so ana-num is noisy at ~1e-6 relative.
fd = struct('name', {}, 'el', {}, 'mode', {}, 'ana', {}, 'num', {}, 'abs_err', {}, 'scale', {});
fid = fopen(filename, 'r');
if fid < 0, error('fdParse: cannot open %s', filename); end
n_bad = 0;
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if strncmp(line, '---', 3)
        k = numel(fd) + 1;
        fd(k).name = strtrim(strrep(line, '---', ''));
        fd(k).el = zeros(0,1); fd(k).mode = cell(0,1);
        fd(k).ana = zeros(0,1); fd(k).num = zeros(0,1); fd(k).abs_err = zeros(0,1);
        continue
    end
    if isempty(fd), continue; end
    t = strsplit(line);
    v = str2double(t);
    label = '';
    if numel(t) == 6, label = t{2}; v(2) = []; end      % mode-labelled row (ft/fc/mt/mc/H)
    if numel(v) ~= 5 || isnan(v(1)) || v(1) ~= round(v(1)), continue; end   % header or blank line
    if ~all(isfinite(v(2:5))), n_bad = n_bad + 1; continue; end              % NaN/Inf row: count it, don't plot it
    fd(k).el(end+1,1) = v(1);   fd(k).mode{end+1,1} = label;
    fd(k).ana(end+1,1) = v(2);  fd(k).num(end+1,1) = v(3);  fd(k).abs_err(end+1,1) = v(4);
end
fclose(fid);
for k = 1:numel(fd)                                          % scale: largest |analytic| among the rows with the same mode label
    [~, ~, g] = unique(fd(k).mode);   fd(k).scale = zeros(size(fd(k).ana));
    for i = 1:max(g), fd(k).scale(g == i) = max(abs(fd(k).ana(g == i))); end
end
if n_bad > 0, warning('fdParse: %d rows with NaN/Inf were skipped in %s', n_bad, filename); end
end