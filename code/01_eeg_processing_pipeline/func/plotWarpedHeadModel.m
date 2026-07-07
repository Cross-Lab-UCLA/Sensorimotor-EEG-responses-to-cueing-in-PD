function plotWarpedHeadModel(EEG, figName)
%% Check plot: subject electrodes co-registered to the standard BEM head model
% Runs EEGLAB's coregister with the stored EEG.dipfit.coord_transform and saves
% a 3-view (left, front, top) PNG of the head mesh + electrodes, matching how
% dipfit applies the transform.
%   green = subject electrodes (warped by coord_transform)
%   brown = reference 10-5 template electrodes on the head model
%
% INPUT:
%   1. EEG     - EEG struct with dipfit settings applied (pop_dipfit_settings)
%   2. figName - full path of the PNG to save (must include subject id)
%
% LM 062226
%%

% --- resolve head model + reference channel files (fall back to dipfit plugin) ---
dipfitDir = fullfile(fileparts(which('pop_dipfit_settings')), 'standard_BEM');
hdm = EEG.dipfit.hdmfile;
if isempty(hdm) || exist(hdm, 'file') ~= 2
    hdm = fullfile(dipfitDir, 'standard_vol.mat');
end
chanfile = EEG.dipfit.chanfile;
if isempty(chanfile) || exist(chanfile, 'file') ~= 2
    chanfile = fullfile(dipfitDir, 'elec', 'standard_1005.elc');
end

% --- render coregister view non-interactively (no blocking, no manual GUI) ---
vis0 = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');
restoreVis = onCleanup(@() set(0, 'DefaultFigureVisible', vis0));
coregister(EEG.chanlocs, chanfile, 'mesh', hdm, ...
    'transform', EEG.dipfit.coord_transform, 'chaninfo1', EEG.chaninfo, ...
    'manual', 'show', 'helpmsg', 'off');
cgFig = findobj(0, 'name', 'coregister()');
cgFig = cgFig(1);  % most recent

% grab the 3-D axes (mesh patch + electrode markers) that coregister just drew
ax = findobj(cgFig, 'type', 'axes');
[~, mi] = max(arrayfun(@(a) numel(get(a, 'Children')), ax));
src  = get(ax(mi), 'Children');
keep = src(ismember(arrayfun(@(c) c.Type, src, 'uni', 0), {'patch', 'line'}));

% --- assemble a clean 3-view montage and save ---
views = {[-90 0], 'Left'; [180 0], 'Front'; [0 90], 'Top'};
fig = figure('Color', 'w', 'Visible', 'off', 'Position', [100 100 1500 540]);
tlo = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for vw = 1:size(views, 1)
    a = nexttile(tlo);
    copyobj(keep, a);
    axis(a, 'equal', 'off');
    view(a, views{vw, 1});
    camlight(a, 'headlight'); lighting(a, 'gouraud'); material(a, 'dull');
    title(a, views{vw, 2}, 'FontWeight', 'bold');
end
sgtitle(tlo, sprintf('Co-registration check: %s  (green = subject, brown = template)', EEG.subject), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

exportgraphics(fig, figName, 'Resolution', 300);
close(fig);
if ~isempty(cgFig) && ishandle(cgFig), close(cgFig); end
end
