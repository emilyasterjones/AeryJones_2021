function out = continuous_specgram(index, excludeperiods, eeg, chinfo, varargin)
% out = calcripspectrum(index, excludeperiods, eeg,ripples,cellinfo, options)
%  Computes the spectrogram around the middle of each decoded event.

%  New version uses non-excluded times to normalize rather than entire
%  epoch

%  Options:
%       fpass-- Determines the frequency range for computing spectrum.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       spectrum_window-- Determines the sliding window used to compute
%           the event triggered spectrogram. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.2]
%       cellfilter--Determines which tetrodes to use for ripple extraction.
%           Default is 'isequal($area, ''CA1'') & numcells>1)'
%  out is a structure with the following fields
%       S-- This is a MxNxT matrix of the spectrogram for each tetrode
%           M is time relative to triggering event, N is frequency, T is event
%       F-- Frequency vector
%       T-- time relative to triggering event
%       fit-- This is the fit based on the spectrum computed for the entire
%           epoch to normalize S. To reconstruct S without normalization,
%           add log10(frequency)*fit(2)+fit(1) to S.
%       index
%from mcarr,  adapted AG 6/13/13 and EJ 10/4/19

%parse the options
params = {};
params.Fs = 1000;
params.fpass = [2 350];
params.trialave = 0;
params.tapers = [3 5];
cwin = [1 1];

% Fs = 1000;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'aveargetrials'
                params.trialave = varargin{option+1};
            case 'probe'
                probecrit = varargin{option+1};
            case 'slidingwin'
                cwin = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% search through chinfo struct to find channels meeting probe criteria
probeindex = evaluatefilter(chinfo,probecrit) ;
out.probeindex = probeindex;

%get valid times for calculation of baseline pwrspectra for zscoring
eprobe = eeg{probeindex(1,1)}{probeindex(1,2)}{probeindex(1,3)};
dur = eprobe.endtime - eprobe.starttime;
Fs = length(eprobe.data)/dur;
fulltimes = eprobe.starttime:1/Fs:eprobe.endtime;
validtimes = ~isExcluded(fulltimes, excludeperiods);
validtimes = validtimes(1:length(eprobe.data));

%iterate through probe channels
parpool('local', 6);
parfor p = 1:size(probeindex,1)
    
    %define EEG for probe
    eprobe = double(eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.data);
    
    % Compute a z-scored spectrogram using the mean and std for the entire session
    [Pprobe,~,freqs{p}] = mtspecgramc(eprobe(validtimes),[cwin(1) cwin(2)],params);
    meanPprobe = mean(Pprobe);
    stdPprobe = std(Pprobe);
    
    Pprobe_tmp{p} = Pprobe;
    meanPprobe_tmp{p} = meanPprobe;
    stdPprobe_tmp{p} = stdPprobe;
end
delete(gcp);

%save to correct structure, separate from parfor loop to prevent
%classification error
for p = 1:size(probeindex,1)
    out.freqs = freqs{p};
    out.Pprobe{p} = Pprobe_tmp{p};
    out.meanPprobe{p} = meanPprobe_tmp{p};
    out.stdPprobe{p} = stdPprobe_tmp{p};
end