function PSD = getPSD(EEG)
%% getPSD
% get dipole PSD info from spectopo and fooof
% save a PSD strucy with variables:
% 1: _db -> PSD from spectopo 
% 2: 'nothing' at the end variable name -> PSD in  uV^2
% 3: _ap_fit -> the fit line from fooof
% 4: _with_fit -> the PSD without the fit removed
% 5: _without_fit -> the PSD with the fit removed

% add fooof path
addpath 'E:\clab\DoD-Gait\code\func\fooof_mat'

for ic_idx = 1:size(EEG.icaact,1)

    PSD(ic_idx).subject      = EEG.subject;
    PSD(ic_idx).group        = EEG.group;
    PSD(ic_idx).chanlocs     = EEG.chanlocs;
    PSD(ic_idx).dipfit       = EEG.dipfit.model(ic_idx);

    event   = EEG.event;
    ic_data = EEG.icaact(ic_idx, :);
    srate   = EEG.srate;
    f_range = [3, 55];
    settings = struct(); % Default FOOOF settings
    settings.peak_width_limits = [2 12]; % set the minimum band-width at twice the freq res (i.e., if bin width is 1Hz, set min to 2)

    stand_start_idx = find(contains({event.type},'standing_start'));
    stand_end_idx = find(contains({event.type},'standing_end'));

    % get standing start and end idxexes
    for s = 1:length(stand_end_idx)
        if any(contains({event(stand_end_idx(s):stand_end_idx(s)+10).type},'walk-run1'))
            stand_nocue1_idx = [stand_start_idx(s) stand_end_idx(s)];
        elseif any(contains({event(stand_end_idx(s):stand_end_idx(s)+10).type},'walk-run2'))
            stand_nocue2_idx = [stand_start_idx(s) stand_end_idx(s)];
        elseif any(contains({event(stand_end_idx(s):stand_end_idx(s)+10).type},'walkAuditory-run1'))
            stand_audi1_idx = [stand_start_idx(s) stand_end_idx(s)];
        elseif any(contains({event(stand_end_idx(s):stand_end_idx(s)+10).type},'walkAuditory-run2'))
            stand_audi2_idx = [stand_start_idx(s) stand_end_idx(s)];
        elseif any(contains({event(stand_end_idx(s):stand_end_idx(s)+10).type},'walkVisual-run1'))
            stand_visual1_idx = [stand_start_idx(s) stand_end_idx(s)];
        elseif any(contains({event(stand_end_idx(s):stand_end_idx(s)+10).type},'walkVisual-run2'))
            stand_visual2_idx = [stand_start_idx(s) stand_end_idx(s)];
        else
            warning('there is match of conditions in the HS')
        end
    end

    % get HS idxexes
    all_idx = find(contains({event.type},'R_HS'));
    [PSD, freqs] = extractPSD(PSD, ic_idx, ic_data, all_idx,...
        event, srate, f_range, settings, 'walk_all');

    walk1_idx = find(strcmp({event.type},'R_HS1_walk-run1'));
    if isempty(walk1_idx)
        warning(['Subject ' EEG.subject ' has no walk-run1 trials'])
    else
        [PSD, freqs] = extractPSD(PSD, ic_idx, ic_data, stand_nocue1_idx,...
            event, srate, f_range, settings, 'stand_nocue1');
        [PSD, freqs] = extractPSD(PSD, ic_idx, ic_data, walk1_idx,...
            event, srate, f_range, settings, 'walk_nocue1');
    end

    walk2_idx = find(strcmp({event.type},'R_HS1_walk-run2'));
    if isempty(walk2_idx)
        warning(['Subject ' EEG.subject ' has no walk-run2 trials'])
    else
        [PSD, freqs] = extractPSD(PSD, ic_idx, ic_data, stand_nocue2_idx,...
            event, srate, f_range, settings, 'stand_nocue2');
        [PSD, freqs] = extractPSD(PSD,ic_idx, ic_data, walk2_idx,...
            event, srate, f_range, settings, 'walk_nocue2');
    end

    audi1_idx = find(strcmp({event.type},'R_HS1_walkAuditory-run1'));
    if isempty(audi1_idx)
        warning(['Subject ' EEG.subject ' has no walkAuditory-run1'])
    else
        [PSD, ~] = extractPSD(PSD, ic_idx, ic_data, stand_audi1_idx,...
            event, srate, f_range, settings, 'stand_audi1');
        [PSD, ~] = extractPSD(PSD,ic_idx, ic_data, audi1_idx,...
            event, srate, f_range, settings, 'walk_audi1');
    end

    audi2_idx = find(strcmp({event.type},'R_HS1_walkAuditory-run2'));
    if isempty(audi2_idx)
        warning(['Subject ' EEG.subject ' has no walkAuditory-run2'])
    else
        [PSD, ~] = extractPSD(PSD, ic_idx, ic_data, stand_audi2_idx,...
            event, srate, f_range, settings, 'stand_audi2');
        [PSD, ~] = extractPSD(PSD,ic_idx, ic_data, audi2_idx,...
            event, srate, f_range, settings, 'walk_audi2');
    end

    visu1_idx = find(strcmp({event.type},'R_HS1_walkVisual-run1'));
    if isempty(visu1_idx)
        warning(['Subject ' EEG.subject ' has no walkVisual-run1'])
        disp('Using walkVisual-run2');
    else
        [PSD, ~] = extractPSD(PSD, ic_idx, ic_data, stand_visual1_idx,...
            event, srate, f_range, settings, 'stand_visual1');
        [PSD, ~] = extractPSD(PSD,ic_idx, ic_data, visu1_idx,...
            event, srate, f_range, settings, 'walk_visual1');
    end

    visu2_idx = find(strcmp({event.type},'R_HS1_walkVisual-run2'));
    if isempty(visu2_idx)
        warning(['Subject ' EEG.subject ' has no walkVisual run-2'])
    else
        [PSD, ~] = extractPSD(PSD, ic_idx, ic_data, stand_visual2_idx,...
            event, srate, f_range, settings, 'stand_visual2');
        [PSD, ~] = extractPSD(PSD,ic_idx, ic_data, visu2_idx,...
            event, srate, f_range, settings, 'walk_visual2');
    end

PSD(ic_idx).freqs = freqs;
end

end

%% functions
function [PSD, f] = extractPSD(PSD, n, ic_data, idx, event, srate, f_range, settings, name_key)

ic_data_stand_nocue = ic_data(...
    floor(event(idx(1)).latency):...
    floor(event(idx(end)).latency));
[spectra,freqs] = spectopo(ic_data_stand_nocue,0,srate,...
    'nfft',srate,'plot','off',...
    'winsize',srate,'overlap',srate*.5); % freq res of 1Hz

lowerF_idx = find(freqs == f_range(1));
higherF_idx = find(freqs == f_range(2));
psd_data = spectra(lowerF_idx:higherF_idx);
f = freqs(lowerF_idx:higherF_idx);

PSD(n).([name_key, '_db']) = psd_data;
PSD(n).(name_key) = 10.^(psd_data / 10); %converting to uV^2
try
    fooof_data = fooof(f, PSD(n).(name_key), f_range, settings, true);
    PSD(n).([name_key, '_with_fit']) = fooof_data.fooofed_spectrum;
    PSD(n).([name_key, '_without_fit']) = fooof_data.fooofed_spectrum - fooof_data.ap_fit;
    PSD(n).([name_key, '_ap_fit']) = fooof_data.ap_fit;
catch
    warning('Unable to Fooof.');
    PSD(n).([name_key, '_with_fit']) = [];
    PSD(n).([name_key, '_without_fit']) = [];
    PSD(n).([name_key, '_ap_fit']) = [];
end
end
