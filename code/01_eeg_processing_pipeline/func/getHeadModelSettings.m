function EEG = getHeadModelSettings(EEG)
%% function to get the tranform for warping to the standard BEM head model in eeglab
% INPUT:
%   1. EEG data
%   2. channel template to match
%
% OUTPUT:
%   1. EEG data with dipfit settings and transforms for auto fitting
%
% LM 072325
%%

for t = 1:length(EEG)
    fprintf(['processing ... ' EEG(t).subject])
    tmpEEG = EEG(t);
    fids = true;

    try tmpEEG.chaninfo.nodatchans(strcmp({tmpEEG.chaninfo.nodatchans.labels},'nas')).labels = 'Nz'; end
    try tmpEEG.chaninfo.nodatchans(strcmp({tmpEEG.chaninfo.nodatchans.labels},'lhj')).labels = 'LPA'; end
    try tmpEEG.chaninfo.nodatchans(strcmp({tmpEEG.chaninfo.nodatchans.labels},'rhj')).labels = 'RPA'; end
    try tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'4Z')).labels = 'Cz'; end
    try tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'8Z')).labels = 'Oz'; end
    if length({tmpEEG.chaninfo.nodatchans.labels}) ~= 3
        fids = false;
    end

    tmpEEG.etc.DigitizedChannels = true;
    tmpEEG.etc.ManuallyFixCoreg = false;

    if fids
        diflist = {tmpEEG.chaninfo.nodatchans.labels 'Cz','Oz'}; % channels to warp by fiducials and Cz
        [newlocs transform] = coregister(tmpEEG.chanlocs, 'standard_1005.elc', 'manual', 'off', ...
        'mesh', 'standard_vol.mat', ...
        'chaninfo1', tmpEEG.chaninfo, ...
        'warp', diflist, ...
        'warpmethod','traditional');

        % check to make sure Cz is close to the top of the head, if not, user
        % adjust
        new_cz_pos = newlocs.pnt((find(strcmp([newlocs.label], 'Cz'))),:);
        new_Oz_pos = newlocs.pnt((find(strcmp([newlocs.label], 'Oz'))),:);
        template_chans = readlocs('standard_1005.elc');
        Cz_idx = find(strcmp({template_chans.labels}, 'Cz'));
        template_cz_pos = [template_chans(Cz_idx).X template_chans(Cz_idx).Y template_chans(Cz_idx).Z];
        Oz_idx = find(strcmp({template_chans.labels}, 'Oz'));
        template_Oz_pos = [template_chans(Oz_idx).X template_chans(Oz_idx).Y template_chans(Oz_idx).Z];

        %%% check Cz %%%
        Cz_distance = pdist2(new_cz_pos,template_cz_pos);
        Cz_distance_Z = pdist2(new_cz_pos(3),template_cz_pos(3));
        Oz_distance = pdist2(new_Oz_pos,template_Oz_pos);
        Oz_distance_Y = pdist2(new_Oz_pos(2),template_Oz_pos(2));

        disp(['Cz distance from template = ' num2str(Cz_distance)]);
        disp(['Cz vertical distance from template = ' num2str(Cz_distance_Z)]);

        disp(['Oz distance from template = ' num2str(Oz_distance)]);
        disp(['Oz vertical distance from template = ' num2str(Oz_distance_Y)]);

        while pdist2(new_cz_pos(3),template_cz_pos(3)) > 5 || pdist2(new_Oz_pos(2),template_Oz_pos(2)) > 5
            % if new dist is more than 5mm off, user fix

            transform = [];
            disp('WARNING!')
            disp('Channels are too far off. Manually fix')
            disp('WARNING!')
            disp(['Current vertical Cz distance difference = ' num2str(Cz_distance_Z)]);
            disp(['Current AP Oz distance difference = ' num2str(Oz_distance_Y)]);

            [newlocs transform] = coregister(tmpEEG.chanlocs, 'standard_1005.elc', 'manual', 'on', ...
                'mesh', 'standard_vol.mat', ...
                'chaninfo1', tmpEEG.chaninfo, ...
                'warp', diflist, ...
                'warpmethod','traditional');

            disp('Fixing...')
            new_cz_pos = newlocs.pnt((find(strcmp([newlocs.label], 'Cz'))),:);
            Cz_distance = pdist2(new_cz_pos,template_cz_pos);
            Cz_distance_Z = pdist2(new_cz_pos(3),template_cz_pos(3));
            disp(['new Cz distance from template = ' num2str(Cz_distance)]);
            disp(['new Cz vertical distance from template = ' num2str(Cz_distance_Z)]);

            new_Oz_pos = newlocs.pnt((find(strcmp([newlocs.label], 'Oz'))),:);
            Oz_distance = pdist2(new_Oz_pos,template_Oz_pos);
            Oz_distance_Y = pdist2(new_Oz_pos(2),template_Oz_pos(2));
            disp(['new Oz distance from template = ' num2str(Oz_distance)]);
            disp(['new Oz AP distance from template = ' num2str(Oz_distance_Y)]);

            tmpEEG.etc.ManuallyFixCoreg = true;
        end

    else
	% Subject HC-17 does not have digitalized channels. Using "globalrescale" to 
	% just the rotate and scale/align the channels to the head model.
        fprintf(['No digitalization of channels for' EEG(t).subject])
        fprintf('...using approximated standard channels.')
    	
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'1L')).labels = 'Fp1'; % as approximatation for scaling
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'2R')).labels = 'Fp2';
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'3LD')).labels = 'M1';
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'3RD')).labels = 'M2';
        
        [newlocs transform] = coregister(tmpEEG.chanlocs, 'standard_1005.elc', 'manual', 'off', ...
        'mesh', 'standard_vol.mat', ...
        'chaninfo1', tmpEEG.chaninfo, ...
        'warp', 'auto', ...
        'warpmethod','globalrescale');

        tmpEEG.etc.DigitizedChannels = false;

        % change channels back
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'Fp1')).labels = '1L';
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'Fp2')).labels = '2R';
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'M1')).labels = '3LD';
        tmpEEG.chanlocs(strcmp({tmpEEG.chanlocs.labels},'M2')).labels = '3RD';
    end

    tmpEEG.etc.warpTransform = transform;
    
    % Use MNI BEM model
    tmpEEG = pop_dipfit_settings(tmpEEG,'coordformat','MNI', 'coord_transform',transform,...
        'model','standardBEM');
    EEG(t) = tmpEEG;
    
    % % plot aligned channels
    % tmpind = find(~cellfun('isempty', { tmpEEG.chanlocs.X }));
    % subplot(1,2,1)
    % h1 = topoplot([],tmpEEG.chanlocs, 'style', 'blank', 'drawaxis', 'on', 'electrodes', ...
    %     'labelpoint', 'chaninfo', tmpEEG.chaninfo); hold on
    % subplot(1,2,2)
    % plotchans3d([ [ tmpEEG.chanlocs(tmpind).X ]' [ tmpEEG.chanlocs(tmpind).Y ]' [ tmpEEG.chanlocs(tmpind).Z ]'], { tmpEEG.chanlocs(tmpind).labels });
    % %load(mesh_file);
    % %ft_plot_mesh(vol.bnd, 'facealpha', 0.5);
    % sgtitle(['Channels Locations on Head model: ' tmpEEG.subject '_' tmpEEG.condition],'Interpreter', 'none');
    % set(gcf, 'units','normalized','outerposition',[.1 .1 .8 .8])
    % set(gcf,'Color',[1 1 1])
    % saveFolder = fullfile('E:\clab\DoD-Gait\reports\data_quality',tmpEEG.subject);
    % mkdir(saveFolder);
    % 
    % saveName = [tmpEEG.subject '_' tmpEEG.condition '_HeadModel_Eloc.png'];
    % saveas(gcf,fullfile(saveFolder,saveName));
    % saveName = [tmpEEG.subject '_' tmpEEG.condition '_HeadModel_Eloc.fig'];
    % saveas(gcf,fullfile(saveFolder,saveName));
    % close all
end
