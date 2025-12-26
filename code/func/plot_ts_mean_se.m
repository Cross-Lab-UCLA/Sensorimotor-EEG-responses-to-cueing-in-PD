function plot_ts_mean_se(data, x_vals, color)
% plot_mean_ci - Plots mean with shaded 95% confidence interval
%
% Usage:
%   plot_mean_ci(data, x_vals)
%
% Inputs:
%   data    - matrix of size [samples x trials] or [samples x numbers]
%   x_vals  - vector of x-axis values (same length as number of rows in data)

if nargin < 3 || isempty(color)
    color = [0 0 0]; 
end

% organize input matrix
if size(x_vals,1) > 1
    x_vals = x_vals';
end
if size(x_vals,2) ~= size(data,1)
    data = data';
end

mean_data = mean(data, 2);
std_data = std(data, 0, 2)/sqrt(size(data,2)); 
upper = mean_data + std_data;
lower = mean_data - std_data;

% Plot shaded area for 95% CI
mean_data = mean_data(:)';
upper = upper(:)';
lower = lower(:)';
fill([x_vals fliplr(x_vals)], [upper fliplr(lower)], ...
    color, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on;

% Plot mean line
plot(x_vals, mean_data, 'Color', color, 'LineWidth', 3);
end