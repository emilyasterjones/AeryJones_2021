function eaj_extractdspikes(directoryname, fileprefix, d)
%based on Dvorak et al, 2021: https://www.biorxiv.org/content/10.1101/2020.07.20.211615v3.full

%filter 5-100Hz
load('\\hub.gladstone.internal\huanglab-lfp\Emily\MATLAB\Filters\ejdsfilter4.mat');

% move to the directory
cd(directoryname);
tmpflist = dir(sprintf('*eeg%02d-*.mat', d));

for i = 1:length(tmpflist) %iterate through all eps, all chans
    load(tmpflist(i).name);
    % get the epoch number
    dash = find(tmpflist(i).name == '-');
    e = str2num(tmpflist(i).name((dash(1)+1):(dash(2)-1)));
    t = str2num(tmpflist(i).name((dash(2)+1:dash(2)+2)));
    
    epoch_eeg = double(eeg{d}{e}{t}.data);
    if sum(abs(epoch_eeg))>0 %if not the reference channel
        
        %filter
        epoch_eeg = double(eeg{d}{e}{t}.data);
        filt_eeg = filtereeg2(eeg{d}{e}{t}, dsfilter);
        filt_eeg = filt_eeg.data(:,1);
        
        %z-score normalize
        norm_eeg = (filt_eeg - nanmean(filt_eeg))/nanstd(filt_eeg);
        
        %find peaks, ignoring first & last 100ms of recording
        peak_inds = findpeaks(norm_eeg);
        peak_inds.loc = peak_inds.loc(peak_inds.loc>100 & peak_inds.loc<length(epoch_eeg)-100);
        
        %find minima
        trough_inds = findpeaks(-1*norm_eeg);
        
        %find closest minimum to each peak
        min_inds = interp1(trough_inds.loc, trough_inds.loc, peak_inds.loc, 'nearest');
        [~,min_trough_ind] = ismember(min_inds, trough_inds.loc);
        
        %calculate amplitudes
        amps = abs(norm_eeg(peak_inds.loc) - norm_eeg(min_inds));
        
        start_ind = NaN(length(peak_inds.loc),1);
        end_ind = NaN(length(peak_inds.loc),1);
        for j=1:length(peak_inds.loc)
            if min_inds(j)>peak_inds.loc(j) %look behind
                start_ind(j) = trough_inds.loc(min_trough_ind(j)-1);
                end_ind(j) = min_inds(j);
            else %look ahead
                start_ind(j) = min_inds(j);
                end_ind(j) = trough_inds.loc(min_trough_ind(j)+1);
            end
        end
        lens = end_ind - start_ind;
        
        %z-score normalize by log transform of amplitudes
        log_amps = log(amps);
        norm_amps = (log_amps - nanmean(log_amps))/nanstd(log_amps);
        
        %pick amps>0.75 & 5ms<len<25ms
        ds_inds = norm_amps>0.75 & lens>=5 & lens<=25;
        
        %save characteristics
        ds.startind = start_ind(ds_inds);
        ds.endind = end_ind(ds_inds);
        ds.peakind = peak_inds.loc(ds_inds);
        ds.starttime = eeg{d}{e}{t}.starttime + ds.startind / eeg{d}{e}{t}.samprate;
        ds.endtime = eeg{d}{e}{t}.starttime + ds.endind / eeg{d}{e}{t}.samprate;
        ds.peaktime = eeg{d}{e}{t}.starttime + ds.peakind / eeg{d}{e}{t}.samprate;
        ds.peak = abs(epoch_eeg(peak_inds.loc(ds_inds)) - epoch_eeg(min_inds(ds_inds)));
        ds.baseline = nanmean(filt_eeg);
        ds.std = nanstd(filt_eeg);
        
    else
        ds.startind = [];
        ds.endind = [];
        ds.peakind = [];
        ds.starttime = [];
        ds.endtime = [];
        ds.peak = [];
        ds.baseline = NaN;
        ds.std = NaN;
    end
    dspikes{d}{e}{t} = ds;
    clear eeg ds;
end
save(sprintf('%s/%sdspikes%02d.mat', directoryname, fileprefix, d), 'dspikes');
end
