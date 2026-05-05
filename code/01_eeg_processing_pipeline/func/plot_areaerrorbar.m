% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar plots the mean and standard deviation of a   %
% set of data filling the space between the positive and negative mean    %
% error using a semi-transparent background, completely customizable.     %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix, with rows corresponding to observations  %
%                   and columns to samples.                               %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.handle:       Figure handle to plot the result.     %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez-Cagigal                                      %
%   Date:    30/04/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com
%   Updated: 
%           02/11/22 Leo -added linestyle, added plot median and
%           bootstrapped 95% CI
%           02/14/23 Leo - added vertical feature, remove mean plot
% ----------------------------------------------------------------------- %
function plot_areaerrorbar(data, options)

    % Default options
    if(nargin<2)
        options.handle     = figure(1);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        %options.color_area = [243 169 114]./255;    % Orange theme
        %options.color_line = [236 112  22]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'c95';
        options.line_style = '-';
        options.vertical = 0; % vertical == 1, time axis is vertical
    end
    if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data,2); end
    options.x_axis = options.x_axis(:);

    if ~isfield(options,'vertical')
        options.vertical = 0;
    end
    
    % Computing the mean and standard deviation of the data matrix
    data_mean = mean(data,1);
    data_median = median(data,1);
    data_std  = std(data,0,1);
    
    % Type of error plot
    switch(options.error)
        case 'std', error = data_std;
        case 'sem', error = (data_std./sqrt(size(data,1)));
        case 'var', error = (data_std.^2);
        case 'c95', error = (data_std./sqrt(size(data,1))).*1.96;
        case 'c95b', error = bootci(10000,{@nanmedian,data},'alpha',.05,'type','per'); 
    end
    
    % Plotting the result
    %figure(options.handle);
    if strcmp(options.error,'c95b')
        x_vector = [options.x_axis', fliplr(options.x_axis')];
        patch = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_area);
        set(patch, 'edgecolor', 'none');
        set(patch, 'FaceAlpha', options.alpha);
        hold on;
%         plot(options.x_axis, data_median, 'color', options.color_line, ...
%             'LineWidth', options.line_width , 'LineStyle', options.line_style);
        %hold off;
    else
        x_vector = [options.x_axis', fliplr(options.x_axis')];

        if options.vertical == 1
            patch = fill([data_mean+error,fliplr(data_mean-error)],x_vector, options.color_area);
        else
            patch = fill(x_vector, [data_mean+error, fliplr(data_mean-error)], options.color_area);
        end

    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, data_mean, 'color', options.color_line, ...
       'LineWidth', options.line_width , 'LineStyle', options.line_style);
    %hold off;
    end
end