%% EEG Dipole Clustering across multiple ICA runs
%
%   1. For every run, dipfit and collect each IC's z-scored activation,
%      spatial filter, scalp map, and dipole. 
%   2. Pool all ICs from all runs, build an all-pairs |corr| similarity on the
%      scalp maps and activation, and agglomerative-cluster into rank clusters.
%   3. Drop poorly-matched clusters: those recovered in fewer than minRuns runs.
%   4. For each good cluster, take the IC closest to the centroid (centrotype =
%      the member most similar to all others) and place it into the
%      dataset-1 (acts as carrier for the IC)
%   5. The cluster's per-run dipole positions and RVs are saved in cltr, carried 
%       downstream for selection
%
% LM 063026
%%
clear all; clc; close all; warning('off','all');
if ispc
    mainDir    = 'C:\clab\DoD-Gait';
else
    mainDir    = '/home/leo/Documents/DoD-Gait';
end
eeglab; ft_defaults; clc; close all;

funcPath    = fullfile(mainDir,'code','01_eeg_processing_pipeline','func');
addpath(funcPath);
reportDir      = fullfile(mainDir,'reports','data_quality');
dataDir     = fullfile(mainDir,'data');
AMICAed_folders = dir(fullfile(dataDir,'*AMICAed*'));
% make sure it is in ascending order
nums = cellfun(@(s) str2double(regexp(s,'\d+$','match','once')), {AMICAed_folders.name});
[~, sortidx] = sort(nums,'ascend');
AMICAed_folders = AMICAed_folders(sortidx);

% check folders has all the subjects
for f = 1:length(AMICAed_folders)
    if length(dir(fullfile(AMICAed_folders(f).folder, AMICAed_folders(f).name,'*.set'))) ~= 58
        warning('There may be missing subjects in the AMICA folder')
        keyboard
    end
end
% label       = extractAfter(AMICAed_folders(1).name,'AMICAed_');
% label       = extractBefore(label,'_01');
% saveFolder =  fullfile(dataDir,['03_bootstrapped_' label]);
saveFolder =  fullfile(dataDir,'03_clustered');
mkdir(saveFolder)
figDir2 = fullfile(reportDir, 'A3_cluster_matching');
if ~exist(figDir2, 'dir'), mkdir(figDir2); end
figDirHM = fullfile(reportDir, 'A3_headmodel_check');
if ~exist(figDirHM, 'dir'), mkdir(figDirHM); end

% parameters
minRuns = 19;  % keep a cluster only if it is recovered in at least this many of the runs (out of nRunsTot)
memberFloor = 0.9;  % drop an IC from a cluster if its R_comb to the cluster centrotype < this (0 = off)

%% co-register
check = true;
eegFiles = dir(fullfile(AMICAed_folders(1).folder, AMICAed_folders(1).name,'*.set'));
EEG = pop_loadset('filename', eegFiles(end).name,  'filepath', eegFiles(end).folder);
warpFile = fullfile(dataDir,'Subject_warpTransformData.mat');

if isfile(warpFile)
    choice = questdlg('Co-registration warp values are found. Would you like to use them?', ...
        'Select', ...
        'Yes','No','No');
    switch choice
        case 'Yes'
            disp('Leading previously used co-registrations');
            load(warpFile); % struct is called T.
            check = false;
        case 'No'
            disp('Performing co-registration check');
    end
end

if check % if true, perform mannual check of head model warps
    T = struct();
    for s = 1:length(eegFiles)
        EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
        if length(EEG.urchanlocs) == 71
            EEG.urchanlocs(65:end-1) = []; % maintain consistency for STUDY.
        end
        EEG = eeg_checkset(EEG);
        EEG = getHeadModelSettings(EEG); % perform warping and checking here

        % save EEG
        %EEG.setname = [EEG.setname ' - coregistered'];
        %pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);

        % save coordinate transform values
        T(s).subject                  = EEG.subject;
        T(s).group                    = EEG.group;
        T(s).warpTransform            = EEG.dipfit.coord_transform;
        T(s).DigitizedChannels        = EEG.etc.DigitizedChannels;
        T(s).ManuallyFixCoreg         = EEG.etc.ManuallyFixCoreg;
    end

    disp('Saving Warp Transform Table.')
    save(fullfile(dataDir, 'Subject_warpTransformData.mat'), 'T');
else
    % add coord_transform from table
    for s = 1:length(eegFiles)
        EEG = pop_loadset('filename', eegFiles(s).name,  'filepath', eegFiles(s).folder);
        if length(EEG.urchanlocs) == 71
            EEG.urchanlocs(65:end-1) = []; % maintain consistency for STUDY.
        end
        EEG = eeg_checkset(EEG);

        subIdx = strcmp({T.subject}, EEG.subject);
        if ~any(subIdx)
            error('Subject %s not found in table T', EEG.subject);
        end

        % if isfield(EEG,'dipfit') && isfield(EEG.dipfit,'coord_transform') ...
        %         && isequal(EEG.dipfit.coord_transform, T(subIdx).warpTransform)
        %     disp('coord_transform is already imported, skipping...')
        %     continue
        % end

        % set up dipfit
        EEG = pop_dipfit_settings(EEG,...
            'coordformat','MNI',...
            'coord_transform',T(subIdx).warpTransform,...
            'model','standardBEM');

        EEG.etc.DigitizedChannels = T(subIdx).DigitizedChannels;
        EEG.etc.ManuallyFixCoreg = T(subIdx).ManuallyFixCoreg;

        %save EEG
        EEG = eeg_checkset(EEG);
        EEG.setname = [EEG.setname ' - coregistered'];
        pop_saveset(EEG, eegFiles(s).name, eegFiles(s).folder);

        % plot warped electrodes on the head model
        plotWarpedHeadModel(EEG, fullfile(figDirHM, ...
            sprintf('%s_warped_headmodel.png', EEG.subject)));
    end

end


%% main loop
nRunsTot = length(AMICAed_folders);
nBefore = nan(length(T),1);   % candidate components per participant (before clustering)
nAfter  = nan(length(T),1);   % retained components (after clustering)

for sub = 1:length(T)
    subject = T(sub).subject;

    %% collect each run's decomposition
    actCell   = cell(1,nRunsTot);   % z-scored activations
    wCell     = cell(1,nRunsTot);   % icaweights
    sphCell   = cell(1,nRunsTot);   % icasphere
    winvCell  = cell(1,nRunsTot);   % icawinv (scalp maps)
    modelCell = cell(1,nRunsTot);   % dipfit.model
    nameCell   = cell(1,nRunsTot);  % set filename per run
    folderCell = cell(1,nRunsTot);  % set folder per run
    dipfitCell = cell(1,nRunsTot);  % dipfit settings per run
    
    delete(gcp('nocreate'))
    parpool("Processes",8)
    parfor run = 1:nRunsTot
        label2     = extractAfter(AMICAed_folders(run).name,'AMICAed_');
        disp(['...' subject ' running ' label2]);
        workingDir = fullfile(AMICAed_folders(run).folder,AMICAed_folders(run).name);
        eegFile    = dir(fullfile(workingDir,['*' subject '*.set']));

        EEG = pop_loadset('filename', eegFile.name,  'filepath', eegFile.folder);
        Tsub = T(strcmp({T.subject},EEG.subject));
        EEG = pop_dipfit_settings(EEG,...
            'coordformat','MNI',...
            'coord_transform',Tsub.warpTransform,...
            'model','standardBEM');
        EEG.etc.DigitizedChannels = Tsub.DigitizedChannels;
        EEG.etc.ManuallyFixCoreg  = Tsub.ManuallyFixCoreg;
        EEG = eeg_checkset(EEG,'ica');   % guarantee icaact = icaweights*icasphere*data
        EEG = pop_multifit(EEG,[],'threshold',100,'rmout','on','plotopt',{'normlen','on'});

        act           = EEG.icaact(:,:);          % [nIC x time]
        actCell{run}  = zscore(act,0,2);          % per-IC z-score along time
        wCell{run}    = EEG.icaweights;
        sphCell{run}  = EEG.icasphere;
        winvCell{run} = EEG.icawinv;
        modelCell{run}= EEG.dipfit.model;
        nameCell{run}   = eegFile.name;
        folderCell{run} = eegFile.folder;
        dipfitCell{run} = EEG.dipfit;
    end
    delete(gcp('nocreate'))

    %carrier = run 1
    carrierName   = nameCell{1};
    carrierFolder = folderCell{1};      % run 1 == _01 folder; this is the carrier
    carrierDipfit = dipfitCell{1};

    %% pool all runs ICs 
    nIC_run = cellfun(@(a) size(a,1), actCell);
    if any(nIC_run ~= nIC_run(1))
        warning('%s: runs have differing IC counts; clustering target = min count.', subject);
    end
    K     = min(nIC_run); % target # clusters = rank
    nChan = size(sphCell{1},2); % # channels feeding ICA

    allMaps = [winvCell{:}]; % [nChan x Ntot] scalp maps
    Ntot    = size(allMaps,2);

    runOf = zeros(Ntot,1); icOf = zeros(Ntot,1);
    e = 0;
    for run = 1:nRunsTot
        for j = 1:nIC_run(run)
            e = e+1; runOf(e) = run; icOf(e) = j;
        end
    end

    % per-estimate dipole position + rv
    allPos = nan(Ntot,3); 
    allRv = nan(Ntot,1);
    for ee = 1:Ntot
        p = modelCell{runOf(ee)}(icOf(ee)).posxyz;
        if ~isempty(p), allPos(ee,:) = p(1,1:3); end
        rj = modelCell{runOf(ee)}(icOf(ee)).rv;
        if ~isempty(rj), allRv(ee) = rj; end
    end

    %% all-pairs similarity on scalp maps (topos)
    R_topo = abs(corr(double(allMaps)));
    R_topo = min(max(R_topo,0),1); 
    R_topo = (R_topo+R_topo')/2;
    allMaps = []; % save ram

    %% all-pairs similarity on activations
    % Combine with the map similarity (R_comb)
    allAct = vertcat(actCell{:});                  % [Ntot x T], member order matches runOf/icOf
    R_act  = abs(corr(allAct'));                 % [Ntot x Ntot] activation correlation
    R_act  = min(max(R_act,0),1); 
    R_act = (R_act+R_act')/2;
    allAct = []; % save ram 
    
    R_comb = min(R_topo, R_act); % similar only if both topo AND activation agree

    %% agglomerative clustering to K groups (on combined map+activation similarity)
    Dmat = 1 - R_comb; Dmat(1:Ntot+1:end) = 0;     % dissimilarity, zero diagonal
    Z   = linkage(squareform(Dmat),'average');
    clu = cluster(Z,'maxclust',K);

    %% membership floor (with rescue to best-clearing cluster)
    % Keep an IC only if it matches some cluster centrotype at R_comb >=
    % memberFloor, and place it in the cluster it matches best. An IC labelled
    % cluster A but more similar to cluster B's centrotype is moved to B rather
    % than dropped (so a strong match is never lost to a poorer assigned label).
    % An IC that clears no centrotype is dropped (clu=0).
    if memberFloor > 0
        ctr = zeros(K,1); % per-cluster centrotype (global IC index)
        for c = 1:K
            mem = find(clu==c);
            if numel(mem)==1
                ctr(c) = mem;
            else
                S = R_comb(mem,mem); % gets similarities among this cluster's members
                S(1:numel(mem)+1:end) = 0; % zero out the diagonal
                [~,kk] = max(sum(S,2)); %find IC that is most similar to all the other members
                ctr(c) = mem(kk); % the representive
            end
        end
        simC = R_comb(:, ctr);                  % [Ntot x K] each IC vs every centrotype
        simC(simC < memberFloor) = -inf;        % keep only clusters the IC clears
        [best, b] = max(simC, [], 2);           % best-clearing cluster per IC
        clu = b;                                % place IC in its best-clearing cluster
        clu(~isfinite(best)) = 0;               % clears no centrotype -> drop
    end

    %% enforce one IC per run per cluster
    % After combined clustering a run may (rarely) contributes >1 IC to a cluster, 
    % but if it does, keep the IC most consistent with the cluster (max R_comb).

    for c = 1:K
        mem = find(clu==c);
        if numel(mem) < 2, continue; end
        Smap = R_comb(mem,mem); Smap(1:numel(mem)+1:end) = 0;   % combined map+act centrotype (anchor)
        [~,kk] = max(sum(Smap,2));  ctrG = mem(kk);
        runs_m = runOf(mem);
        for r = unique(runs_m)'
            idx = find(runs_m == r);                        % this run's members (positions in mem)
            if numel(idx) > 1
                [~,best] = max(R_comb(mem(idx), ctrG));      % keep run's IC consistent in map AND act
                surplus  = idx; surplus(best) = [];
                clu(mem(surplus)) = 0;                       % drop this run's extra ICs
            end
        end
    end

    %% per cluster find coverage flag + refine repersentive IC (closest to centroid)
    cltr = struct([]);
    newW = []; newWinv = []; newModel = modelCell{1}([]);
    keptClusters = []; droppedClusters = []; ctrRun = [];
    slot = 0;
    for c = 1:K
        mem      = find(clu==c);
        runsRep  = numel(unique(runOf(mem)));

        % centrotype = member most similar to all others (combined map+activation)
        Ssub = R_comb(mem,mem); 
        Ssub(1:numel(mem)+1:end) = 0;
        [~,kk] = max(sum(Ssub,2));
        ctr  = mem(kk);         % IC closest to the centroid in combined map+act space
        crun = runOf(ctr);      % IC closest to the centroid's run
        cic = icOf(ctr);        % IC closest to the centroid's IC index

        if runsRep < minRuns                       % not recovered in enough runs -> drop
            droppedClusters(end+1,:) = [c runsRep numel(mem)];
            continue
        end

        slot = slot + 1;
        newW(slot,:)    = wCell{crun}(cic,:) * sphCell{crun};   % sphere folded into the filter row [1 x nChan]
        newWinv(:,slot) = winvCell{crun}(:,cic);               % the IC's own scalp map
        newModel(slot)  = modelCell{crun}(cic);                % the IC's own dipole

        cltr(slot).pos_runs        = allPos(mem,:);
        cltr(slot).rv_runs         = allRv(mem);
        cltr(slot).run_of          = runOf(mem);        % run index for each member
        cltr(slot).ic_of           = icOf(mem);         % per-member original IC index
        cltr(slot).posxyz          = allPos(ctr,:);     % IC closest to the centroid position (point estimate)
        cltr(slot).rv              = allRv(ctr);
        cltr(slot).icawinv         = winvCell{crun}(:,cic);
        cltr(slot).centrotype_run  = crun;
        cltr(slot).centrotype_ic   = cic;
        cltr(slot).nMembers        = numel(mem);
        cltr(slot).runsRepresented = runsRep;

        % reliability measures
        simWithin = R_comb(mem,mem);              % pairwise combined map+act |corr| among this cluster's members
        od = ~eye(numel(mem),'logical');          % off-diagonal mask (exclude self-pairs, which are always 1)
        cltr(slot).clusterCorr = mean(simWithin(od),'omitnan');   % mean within-cluster combined map+act |corr| (1 = tight/coherent, lower = looser)
        Pm = allPos(mem,:); Pm = Pm(all(~isnan(Pm),2),:);         % valid per-run dipole positions
        if size(Pm,1) >= 2
            cltr(slot).posSpread = mean(sqrt(sum((Pm-mean(Pm,1)).^2,2)));  % mean dist from centroid
        else
            cltr(slot).posSpread = NaN;
        end

        keptClusters(end+1) = c;
        ctrRun(end+1)       = crun;
    end

    %% cluster diagnostic figure
    plotClusterDiagnostic(clu, R_comb, allPos, runOf, minRuns, subject, ...
        fullfile(figDir2, sprintf('%s_clusters.fig', subject)), carrierDipfit.hdmfile);

    close all;
    if slot == 0
        warning('%s: no clusters reached minRuns = %d; skipping.', subject, minRuns);
        keyboard
    end

    %% inject into the dataset-1 carrier (shared data/chanlocs/events)
    EEG = pop_loadset('filename', carrierName, 'filepath', carrierFolder);
    EEG.icaweights   = newW;            % [nGood x nChan], each run's sphere already folded in
    EEG.icasphere    = eye(nChan);      % identity: every row carries its own sphere
    EEG.icawinv      = newWinv;
    EEG.icaact       = [];
    EEG.dipfit       = carrierDipfit;   % dipfit settings from run 1
    EEG.dipfit.model = newModel;
    EEG.dipfit.cltr  = cltr;
    EEG = eeg_checkset(EEG,'ica');      % recompute icaact = newW*I*data (exact per injected IC)

    % diagnostics
    EEG.etc.cluster_info = struct( ...
        'K',               K, ...
        'nKept',           slot, ...
        'keptClusters',    keptClusters, ...
        'droppedClusters', droppedClusters, ...   % [clusterID runsRepresented nMembers]
        'centrotypeRun',   ctrRun, ...
        'minRuns',         minRuns);
    fprintf('%s: %d/%d clusters kept; centrotypes drawn from %d distinct runs (dropped %d).\n', ...
        subject, slot, K, numel(unique(ctrRun)), size(droppedClusters,1));
    nBefore(sub) = K;        % components before clustering (candidate clusters)
    nAfter(sub)  = slot;     % components kept after clustering

    setname = extractBefore(EEG.filename,'.set');
    EEG.setname = [setname ' - dipfitted - clustered'];
    pop_saveset(EEG, EEG.filename, saveFolder);
end

%% summarize component counts across participants and export CSV
nExcluded = nBefore - nAfter;
countsT = table(string({T.subject}'), nBefore, nAfter, nExcluded, ...
    'VariableNames', {'subject','nBefore','nKept','nExcluded'});
writetable(countsT, fullfile(reportDir, 'A3_cluster_component_counts.csv'));
fprintf('Components per participant: before %.1f (SD %.1f) | kept %.1f (SD %.1f) | excluded %.2f (SD %.2f)\n', ...
    mean(nBefore,'omitnan'), std(nBefore,'omitnan'), ...
    mean(nAfter,'omitnan'),  std(nAfter,'omitnan'), ...
    mean(nExcluded,'omitnan'), std(nExcluded,'omitnan'));
disp('done.')

%% functions
function plotClusterDiagnostic(clu, Rsim, allPos, runOf, minRuns, subject, savePath, hdmfile)
% Two-panel diagnostic for IC clustering.
%   Shows, per subject, how the pooled across-run ICs were grouped and which
%   clusters were dropped for low run-coverage.
%
%   Panel 1 (similarity map): 2-D MDS embedding of the dissimilarity
%     1-Rsim, so distance between points approximates 1 - |corr| in the
%     same (combined map+activation) space the clustering and the representatives
%     were built on. Axes are symmetric/origin-centred with a dashed crosshair.
%   Panel 2 (dipole positions): the same estimates in MNI space, inside a
%     translucent brain (drawn from hdmfile).
%
%   Inputs:
%     clu       [Ntot x 1] cluster label per pooled IC estimate (0 = enforcement-dropped)
%     Rsim      [Ntot x Ntot] combined map+activation |corr| similarity the
%               clustering and representative selection were built on (R_comb)
%     allPos    [Ntot x 3] each estimate's dipole position (MNI mm; NaN ok)
%     runOf     [Ntot x 1] run index each estimate came from
%     minRuns   keep a cluster if it appears in >= this many runs
%     subject   char, used in the figure title
%     savePath  (optional) base path; if given, saved as both .png and .fig (extension ignored)
%     hdmfile   (optional) EEG.dipfit.hdmfile; adds the translucent brain

K = max(clu);
runsRep = zeros(K,1); ctr = zeros(K,1);
for c = 1:K
    mem = find(clu==c);
    if isempty(mem), continue; end
    runsRep(c) = numel(unique(runOf(mem)));
    Ssub = Rsim(mem,mem); Ssub(1:numel(mem)+1:end) = 0;
    [~,kk] = max(sum(Ssub,2)); ctr(c) = mem(kk);     % centrotype (combined map+act), matches main loop
end
kept    = runsRep >= minRuns;
keepIdx = ismember(clu, find(kept));                 % per-estimate kept mask
dropIdx = ~keepIdx;

% 2-D MDS embedding (computed here so colors can be ordered by spatial position)
D = double(1 - Rsim);            % dissimilarity in the clustering space
D(~isfinite(D)) = 1;             % guard NaN/Inf (e.g. a flat IC)
D = (D + D')/2;                  % enforce symmetry
D(1:size(D,1)+1:end) = 0;        % zero diagonal
Y = cmdscale(D, 2);
if size(Y,2) < 2, Y(:,end+1:2) = 0; end   % pad if rank-deficient

% set up colors
keptIDs  = find(kept);
nKept    = numel(keptIDs);
% assign neighbors maximally different colors.
ang      = atan2(Y(ctr(keptIDs),2), Y(ctr(keptIDs),1));
[~,ord]  = sort(ang);
keptIDs  = keptIDs(ord);
hues     = mod((0:nKept-1)' * 0.3819, 1);
cmKept   = hsv2rgb([hues, repmat(0.9,nKept,1), repmat(0.95,nKept,1)]);
colorOf  = zeros(K,3);
colorOf(keptIDs,:) = cmKept;
estColor = zeros(numel(clu),3); % black = dropped / enforcement-removed (clu==0)
valid    = clu > 0;
estColor(valid,:) = colorOf(clu(valid),:);

% set up fig
fig = figure('Color','w','Visible','on','Position',[100 100 1100 480]);
tlo = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

% panel 1: similarity map (2-D MDS of combined map+activation dissimilarity)
ax = nexttile(tlo);
hold(ax,'on');

% transparent footprint around each kept cluster (hull, or a disk if degenerate/tight)
m = max(abs([Y(:,1); Y(:,2)]));
if ~isfinite(m) || m==0, m = 1; end
minRad = 0.03*m;   % minimum drawn radius
th = linspace(0,2*pi,40)';
for c = keptIDs'
    P   = Y(clu==c, :);
    col = colorOf(c,:);
    okHull = size(P,1) >= 3 && max(max(P,[],1)-min(P,[],1)) >= minRad;   % enough points and 2-D extent
    if okHull
        try 
            kk = convhull(P(:,1),P(:,2)); 
        catch
            okHull = false; 
        end
    end
    if okHull % normal cluster: true convex hull
        patch(ax, P(kk,1), P(kk,2), col, 'FaceAlpha',0.5, ...
            'EdgeColor',col, 'EdgeAlpha',0.5, 'LineWidth',0.5);
    else
        cen = mean(P,1);
        patch(ax, cen(1)+minRad*cos(th), cen(2)+minRad*sin(th), col, 'FaceAlpha',0.6, ...
            'EdgeColor',col, 'EdgeAlpha',0.6, 'LineWidth',0.5);
    end
end

xline(ax, 0, '--', 'Color',[.7 .7 .7]);    % dashed zero reference (drawn behind points)
yline(ax, 0, '--', 'Color',[.7 .7 .7]);

scatter(ax, Y(dropIdx,1), Y(dropIdx,2), 28, estColor(dropIdx,:), 'v', 'filled', 'MarkerFaceAlpha',0.7);   % black down-triangle = dropped IC
scatter(ax, Y(keepIdx,1), Y(keepIdx,2), 24, estColor(keepIdx,:), 'filled');
scatter(ax, Y(ctr(kept),1), Y(ctr(kept),2), 80, 'k', 'o', 'LineWidth',1.5);
axis(ax,'equal');
xlim(ax,[-m m]); ylim(ax,[-m m]);

xlabel(ax,'MDS dim 1 (a.u.)'); ylabel(ax,'MDS dim 2 (a.u.)');
title(ax,'Similarity - Multidimensional Scaling  Map');
subtitle(ax,'each dot = one IC from one run (combined map+activation similarity)', 'FontSize',10);

% legend via off-screen proxy markers (so colours/markers are explained)
hCtr  = plot(ax, nan, nan, 'ko', 'MarkerFaceColor','none', 'LineWidth',1.5, 'MarkerSize',9);
hKept = plot(ax, nan, nan, 'o',  'MarkerFaceColor',[0 0.45 0.74], 'MarkerEdgeColor','none', 'MarkerSize',8);
hDrop = plot(ax, nan, nan, 'v',  'MarkerFaceColor','k', 'MarkerEdgeColor','none', 'MarkerSize',8);
legend(ax, [hKept hDrop hCtr], {'clustered IC (color = cluster)', 'dropped IC', 'representative IC'}, ...
    'Location','best', 'FontSize',8, 'Box','off');

% panel 2: dipole positions
ax = nexttile(tlo); hold(ax,'on');
if nargin >= 8 && ~isempty(hdmfile) && exist(hdmfile,'file')
    S = load(hdmfile);                       % standard BEM volume
    if isfield(S,'vol'), vol = S.vol; else, vol = S; end
    bnd = vol.bnd;
    [~,bi] = min(arrayfun(@(b) prod(range(b.pnt,1)), bnd));  % innermost surface = brain
    patch(ax, 'Vertices',bnd(bi).pnt, 'Faces',bnd(bi).tri, ...
        'FaceColor',[.8 .8 .85], 'FaceAlpha',0.08, 'EdgeColor','none');
    camlight(ax,'headlight'); lighting(ax,'gouraud'); material(ax,'dull');
end
dp = dropIdx & all(~isnan(allPos),2);
kp = keepIdx & all(~isnan(allPos),2);
scatter3(ax, allPos(dp,1),allPos(dp,2),allPos(dp,3), 28, estColor(dp,:), 'v', 'filled', 'MarkerFaceAlpha',0.55); % black down-triangle = dropped IC
scatter3(ax, allPos(kp,1),allPos(kp,2),allPos(kp,3), 24, estColor(kp,:), 'filled', 'MarkerFaceAlpha',0.55);
ctrKept = ctr(kept);                                   % representative IC per kept cluster (panel 1 rings)
cc = ctrKept(all(~isnan(allPos(ctrKept,:)),2));        % keep only those with a valid dipole (panel 2 rings)
if numel(cc) < numel(ctrKept)                          % panel 2 will show fewer rings than panel 1
    missingCtr  = ctrKept(any(isnan(allPos(ctrKept,:)),2));   % IC indices whose dipole is NaN (out-of-head)
    missingRun  = runOf(missingCtr);                          % run each came from
    missingIC   = icOf(missingCtr);                           % original IC index within that run
    warning('%s: %d of %d representative dipoles are NaN; panel 2 shows fewer rings than panel 1.', ...
        subject, numel(ctrKept)-numel(cc), numel(ctrKept));
    keyboard   % inspect missingCtr / missingRun / missingIC
end
scatter3(ax, allPos(cc,1),allPos(cc,2),allPos(cc,3), 80, 'k', 'o', 'LineWidth',1.5);
axis(ax,'equal','vis3d'); view(ax,135,20); grid(ax,'on');
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z'); title(ax,'dipole positions');

sgtitle(tlo, sprintf('%s — cluster diagnostic  (kept %d / %d)', subject, nKept, K), 'Interpreter','none');
if nargin >= 7 && ~isempty(savePath)
    [pth,nm] = fileparts(savePath);
    base = fullfile(pth, nm);
    exportgraphics(fig, [base '.png'], 'Resolution', 600);
    savefig(fig, [base '.fig']);
end
end
