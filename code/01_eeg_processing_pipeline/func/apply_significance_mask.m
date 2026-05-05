function [masked_tvals, binary_mask] = apply_significance_mask(pvals, tvals, alpha)

    % set default alpha to 0.05 if not provided
    if nargin < 3
        alpha = 0.05;
    end

    if ~isequal(size(pvals), size(tvals))
        error('The size of p-value and T-value arrays must match.');
    end

    binary_mask = (pvals < alpha) & ~isnan(pvals);
    masked_tvals = tvals .* binary_mask;
    masked_tvals(~binary_mask) = NaN;
end