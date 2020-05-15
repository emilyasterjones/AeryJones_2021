function out = RTcoh(index, excludeperiods, eeg, ripples, chinfo, varargin)

% from mcarr riptriggeredcoherence(directoryname,fileprefix,days,varargin)
%  Computes and saves the coherence around each ripple. 10/3/13

%  Options:
%       fpass-- Determines the frequency range for computing coherence.
%           Default: [2 350]
%       average_trials-- Determines if events are averaged or not.
%           Default: 0
%       coherence_window-- Determines the sliding window used to compute
%           the event triggered coherence. Default: [0.1 0.01]
%       event_window--Determines the size of the window around each
%           triggering event. Default: [0.2 0.4]

%parse the options
Fs = 1000;
params = {};
params.Fs = Fs;
params.fpass = [2 350];
params.trialave = 0;
params.tapers = [3 5];
win = [0.4 0.4];
cwin = [0.1 0.1];  %smoothed for visualization, [.1 .1] for quantification
minstd= 3;
interripwin = 0;



for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'fpass'
                params.fpass = varargin{option+1};
            case 'averagetrials'
                params.trialave = varargin{option+1};
            case 'cwin'
                cwin = varargin{option+1};
            case 'event_window'
                win = varargin{option+1};
            case 'minstd'
                minstd = varargin{option+1};
            case 'trig'
                trigcrit = varargin{option+1};
            case 'coh'
                cohcrit = varargin{option+1};
            case 'probe'
                probecrit = varargin{option+1};
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
cohindex = evaluatefilter(chinfo, cohcrit);
probeindex = evaluatefilter(chinfo,probecrit);

out.trigindex = trigindex;
out.cohindex = cohindex;
out.probeindex = probeindex;

%iterate through probe channels

r = ripples{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};

if ~isempty(cohindex)  %if coh channels are present in this animal
    
    for t = 1:size(cohindex,1) %iterate through coh channels
        
        [valid, ~] = getvalidrips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);
        if numel(valid)>0 %if valid ripples are detected
            
            ripplestdout = r.maxthresh(valid);
            %changed 10-11-13 to trigger from thresh crossing
            triggers = r.threshind(valid)/Fs;
            starttime = r.timerange(1);
            endtime = r.timerange(2);
            
            % Define EEG for trig
            ecoh = double(eeg{cohindex(t,1)}{cohindex(t,2)}{cohindex(t,3)}.data);%*12500/65536
            
            for p = 1:size(probeindex,1)  %iterate through probe chans
                
                %define EEG for probe
                eprobe = double(eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.data);%*12500/65536
                
                % Calculate the event triggered coherence
                data1 = createdatamatc(ecoh,triggers,params.Fs,[win(1) win(2)]);
                data2 = createdatamatc(eprobe,triggers,params.Fs,[win(1) win(2)]);
                
                [C,phi,~,~,~,times,freqs] = cohgramc(data1,data2,[cwin(1) cwin(2)],params);  % S12 S1 S2 supressed with ~
                
                out.coherence{t}{p} = C; %im2uint8(C);
                out.phase{t}{p} = phi;
                out.times = times-win(1);
                out.freqs = freqs;
                out.ripsizes = ripplestdout;
                
                % Calculate continuous coherence across all nonexcluded periods
                bins = starttime+win(1):cwin(1):endtime-win(2);
                bininds = (~isExcluded(bins-win(1),excludeperiods) & ~isExcluded(bins+win(2),excludeperiods));
                bins = bins(bininds)-starttime;
                data1 = createdatamatc(ecoh,bins,params.Fs,[win(1) win(2)]);
                data2 = createdatamatc(eprobe,bins,params.Fs,[win(1) win(2)]);
                
                C = cohgramc(data1,data2,[cwin(1) cwin(2)],params);
                out.meanCoh{t}{p} = mean(squeeze(mean(C)),2);
                out.stdCoh{t}{p} = std(squeeze(std(C)),0,2)';
            end
            out.riptimes = triggers+starttime;
        else
            out.coherence{t} = 0; %im2uint8(C);
            out.phase{t} = 0;
            out.times = 0;
            out.freqs = 0;
            out.riptimes{t} = 0;
            out.ripsizes = 0;
        end
    end
else
    out.coherence{1} = 0; %im2uint8(C);
    out.phase{1} = 0;
    out.times = 0;
    out.freqs = 0;
    out.riptimes{1} = 0;
    out.ripsizes = 0;
end
