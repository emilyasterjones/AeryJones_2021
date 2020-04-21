function out = calccoherence(index, excludetimes, eeg, chinfo, varargin)
%Adapted from MC calccoherence, AG ag_coherence, and EJ RTcoh
%by EJ 10/18/2019

% function out = calccoherence(index, excludetimes, eeg, varargin)
%
%  Plots the coherence for an eeg tetrode pair. If you use a time filter,
%  excluded times are removed and the includedtimes are averaged together.
%
%   out is a structure with the following fields
%       coherence-- Coherence score between electrode pairs
%       frequency-- Frequency vector
%       index-- Only if appendindex is set to 1 (default)


% parse the options
params = [];
params.trialave = 0;
params.Fs = 1000;
params.fpass = [5 110];
params.tapers = [3 5];

appendindex = 1;
% window =[10 0.5];

for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'Fs'
                params.Fs = varargin{option+1};
            case 'fpass'
                params.fpass = varargin{option+1};
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

%find electrodes
cohindex = evaluatefilter(chinfo,cohcrit);
probeindex = evaluatefilter(chinfo,probecrit);
out.cohindex = cohindex;
out.probeindex = probeindex;

%get valid times for calculation of baseline pwrspectra for zscoring
eprobe = eeg{probeindex(1,1)}{probeindex(1,2)}{probeindex(1,3)};
dur = eprobe.endtime - eprobe.starttime;
Fs = length(eprobe.data)/dur;
fulltimes = eeg{probeindex(1,1)}{probeindex(1,2)}{probeindex(1,3)}.starttime:1/Fs:eeg{probeindex(1,1)}{probeindex(1,2)}{probeindex(1,3)}.endtime;
validtimes = ~isExcluded(fulltimes, excludetimes);
validtimes = validtimes(1:length(eprobe.data));
validtimes = find(validtimes);
timebin = 1:1000:length(validtimes);

%calculate coherence between each pair
parpool('local',6);
if ~isempty(cohindex)  %if coh channels are present in this animal
    parfor t = 1:size(cohindex,1) %iterate through coh channels
        ecoh = double(eeg{cohindex(t,1)}{cohindex(t,2)}{cohindex(t,3)}.data);
        for p = 1:size(probeindex,1)  %iterate through probe chans
            eprobe = double(eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.data);
            [C,~,~,~,~,f] = coherencyc(ecoh(validtimes), eprobe(validtimes), params);
            coherence{t}{p} = C;
            freqs{t}{p} = f;
            
            for b = 1:length(timebin)-1
                validbin = validtimes(timebin(b):timebin(b+1));
                [C,~,~,~,~,f] = coherencyc(ecoh(validbin), eprobe(validbin), params);
                coherence_bin{t}{p}(b,:) = C;
                freqs_bin{t}{p} = f;
            end
            
        end
    end
else
    out.coherence{1} = 0;
    out.freqs = 0;
end
delete(gcp);


%save to correct structure, separate from parfor loop to prevent
%classification error
for t = 1:size(cohindex,1)
    for p = 1:size(probeindex,1)
        out.freqs = freqs{1}{1};
        out.coherence{t}{p} = coherence{t}{p};
        out.freqs_bin{t}{p} = freqs_bin{t}{p};
        out.coherence_bin{t}{p} = coherence_bin{t}{p};
    end
end

if (appendindex)
    out.index = index;
end

end