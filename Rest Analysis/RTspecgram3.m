function out = RTspecgram3(index, excludeperiods, eeg, ripples,chinfo, varargin)
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
%from mcarr,  adapted AG 6/13/13

%parse the options
params = {};
params.Fs = 1000;
params.fpass = [2 350];
params.trialave = 0;
params.tapers = [3 5];
win = [.4 .4];
cwin = [0.1 0.1];   %with sliding: [.1 .01] OR non overlapping: [.1 .1]
minstd = 5;    %can change minstd!!
interripwin = 0;  %does not eliminate ripples that occur in doublets or triplets

Fs = 1000;

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'aveargetrials'
                params.trialave = varargin{option+1};
            case 'spectrum_window'
                cwin = varargin{option+1};
            case 'event_window'
                win = varargin{option+1};
            case 'minstd'
                minstd = varargin{option+1};
            case 'trig'
                trigcrit = varargin{option+1};
            case 'probe'
                probecrit = varargin{option+1};
            case 'slidingwin'
                cwin = varargin{option+1};
            case 'interripwin'
                cwin = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

% search through chinfo struct to find channels meeting trig and probe
% criteria
trigindex = evaluatefilter(chinfo,trigcrit);
probeindex = evaluatefilter(chinfo,probecrit) ;

out.trigindex = trigindex;
out.probeindex = probeindex;

%get valid times for calculation of baseline pwrspectra for zscoring
fulltimes = eeg{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)}.starttime:1/Fs:eeg{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)}.endtime;
validtimes = ~isExcluded(fulltimes, excludeperiods);


%iterate through probe channels 
for t = 1:size(trigindex,1)
    
    %get ripple trigger times
    r = ripples{trigindex(t,1)}{trigindex(t,2)}{trigindex(t,3)};
    
    [valid, excluded] = getvalidrips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);
    %[valid, excluded] = getnoiserips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);

    ripplestdout = r.maxthresh(valid);
    triggers = r.threshind(valid)/Fs;
    %changed 10/8/13 to use thresh crossing as TRIGGERS instead of start
    trigstats = excluded;
    % Define EEG for trig
    etrig = double(eeg{trigindex(t,1)}{trigindex(t,2)}{trigindex(t,3)}.data);
    
     if (numel(triggers) > 0)% && (length(eeg{trigindex(t,1)}{trigindex(t,2)}{trigindex(t,3)}.data)>=length(fulltimes))
        if (length(etrig)<length(validtimes))
            fprintf('Length issue: %d < %d\n',length(etrig), length(validtimes));
            validtimes = validtimes(1:length(etrig));
        end
        % Calculate the event triggered spectrogram for TRIGGER
        [Strig,times,freqs] = mtspecgramtrigc(etrig,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);
        out.freqs = freqs; 
        out.times = times - win(1); % SHIFTS TO CENTER AROUND 0
        
        % Compute a z-scored spectrogram using the mean and std for the entire
        % session  for TRIGGER
        Ptrig = mtspecgramc(etrig(validtimes),[cwin(1) cwin(1)],params);
        meanPtrig = mean(Ptrig);
        stdPtrig = std(Ptrig);
        
        for i = 1:size(Strig,1) %freqs
            for j = 1:size(Strig,3)
                Strig(i,:,j) = (Strig(i,:,j) - meanPtrig)./stdPtrig;
            end
        end
        
        out.Strig{t} = Strig;
        out.stds{t} = ripplestdout;
        out.meanPtrig{t} = meanPtrig;
        out.stdPtrig{t} = stdPtrig;
        out.riptimes = triggers;
        out.trigstats = trigstats;
        %iterate through probe channels
        
        parfor p = 1:size(probeindex,1)
            
            %define EEG for probe
            eprobe = double(eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.data);
            
            % Calculate the event triggered spectrogram for PROBE
            [Sprobe,~,~] = mtspecgramtrigc(eprobe,triggers,[win(1) win(2)],[cwin(1) cwin(2)],params);
            % Compute a z-scored spectrogram using the mean and std for the entire session
            Pprobe = mtspecgramc(eprobe(validtimes),[cwin(1) cwin(1)],params);
            meanPprobe = mean(Pprobe);
            stdPprobe = std(Pprobe);
            
            for i = 1:size(Sprobe,1)
                for j = 1:size(Sprobe,3)
                    Sprobe(i,:,j) = (Sprobe(i,:,j) - meanPprobe)./stdPprobe;
                end
            end
            
            Pprobe_tmp{p} = Pprobe;
            Sprobe_tmp{p} = Sprobe;
            meanPprobe_tmp{p} = meanPprobe;
            stdPprobe_tmp{p} = stdPprobe;
            
        end
        delete(gcp);
        
        %save to correct structure, separate from parfor loop to prevent
        %classification error
        for p = 1:size(probeindex,1)
            out.Pprobe{t}{p} = Pprobe_tmp{p};
            out.Sprobe{t}{p} = Sprobe_tmp{p};
            out.meanPprobe{t}{p} = meanPprobe_tmp{p};
            out.stdPprobe{t}{p} = stdPprobe_tmp{p};
        end

    else  % trigs = 0, no valid ripples
        out.Strig{t} = 0;
        out.stds{t} = 0;
        out.Sprobe{t} = 0;
        out.freqs = 0;
        out.times = 0;
        out.riptimes = 0;
    end
    
end