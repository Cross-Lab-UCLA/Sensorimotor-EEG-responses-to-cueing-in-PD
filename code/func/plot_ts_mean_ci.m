function plot_ts_mean_ci(data, x_vals, color,boot_num)
% plot_mean_ci - Plots mean with shaded 95% confidence interval
%
% Usage:
%   plot_mean_ci(data, x_vals)
%
% Inputs:
%   data    - matrix of size [samples x trials] or [samples x numbers]
%   x_vals  - vector of x-axis values (same length as number of rows in data)

if nargin < 4 || isempty(boot_num)
    boot_num = []; 
end

% organize input matrix
if size(x_vals,1) > 1
    x_vals = x_vals';
end
if size(x_vals,2) ~= size(data,1)
    data = data';
end

if ~isempty(boot_num)

    nSamples = size(data, 2);  % Number of trials
    nTimes = size(data, 1);    % Number of time points
    boot_means = zeros(nTimes, boot_num);  % Preallocate

    for b = 1:boot_num
        % Resample trials with replacement
        idx = randi(nSamples, [nSamples, 1]);
        resampled = data(:, idx);

        % Compute mean across resampled trials
        boot_means(:, b) = mean(resampled, 2);
    end

    lower = prctile(boot_means, 2.5, 2);  % 2.5th percentile (lower bound)
    upper = prctile(boot_means, 97.5, 2); % 97.5th percentile (upper bound)
    mean_data = mean(data, 2);

else
    mean_data = mean(data, 2);
    sem = std(data, 0, 2) / sqrt(size(data, 2));
    CI95 = 1.96 * sem;

    upper = mean_data + CI95;
    lower = mean_data - CI95;
end


% Plot shaded area for 95% CI
mean_data = mean_data(:)';
upper = upper(:)';
lower = lower(:)';

fill([x_vals fliplr(x_vals)], [upper fliplr(lower)], ...
    color, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on;

% Plot mean line
%plot(x_vals, mean_data', 'Color', color, 'LineWidth', 3);
end