function p = percentile(x, pct)
    % Replacement for prctile(x, pct) without the Statistics Toolbox
    x = sort(x(:));              % flatten and sort
    n = numel(x);                % total elements
    idx = ceil(pct / 100 * n);   % percentile index
    idx = max(1, min(n, idx));   % clamp within bounds
    p = x(idx);                  % output value
end