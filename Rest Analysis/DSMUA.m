function out = DSMUA(index, excludeperiods, mua, dspikes, chinfo, varargin)
%parse the options
win = [-0.4 0.4];
binsize = 0.1;

Fs = 1000;
spikeFs = 30000;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'event_window'
                win = varargin{option+1};
            case 'binsize'
                binsize = varargin{option+1};
            case 'trig'
                trigcrit = varargin{option+1};
            case 'probe'
                probecrit = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% search through chinfo struct to find channels meeting trig and probe criteria
trigindex = evaluatefilter(chinfo,trigcrit);
probeindex = evaluatefilter(chinfo,probecrit) ;
out.trigindex = trigindex;
out.probeindex = probeindex;

fulltimes = mua{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)}.starttime:1/Fs:mua{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)}.endtime;
validtimes = ~isExcluded(fulltimes, excludeperiods);

%iterate over trigger channels
for t = 1:size(trigindex,1)
    
    %get dspike trigger times
    d = dspikes{trigindex(t,1)}{trigindex(t,2)}{trigindex(t,3)};
    valid = find(~isExcluded(d.starttime,excludeperiods) & ~isExcluded(d.endtime,excludeperiods));
    
    d.starttime = d.starttime(valid);
    d.endtime = d.endtime(valid);
    d.peaktime = mua{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)}.starttime + d.peakind / Fs;
    triggers = d.peaktime(valid);
    out.dstimes{t} = triggers;
    
    %iterate through probe channels
    for p = 1:size(probeindex,1)
        pmua = mua{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)};
        spiketimes = pmua.spiketimes/spikeFs;
        bins = fulltimes(1):binsize:fulltimes(end);
        validinds = ~isExcluded(spiketimes,excludeperiods);
        
        bininds = ~isExcluded(bins,excludeperiods);
        out.baserate{t}{p} = sum(validinds)/sum(validtimes)*Fs; %FR over all immob
        out.sdbaserate{t}{p} = std(histcounts(spiketimes(validinds),bins(bininds))/binsize);
        
        bins = win(1):binsize:win(2);
        out.psth{t}{p} = NaN(length(triggers),(length(bins)-1)); %no 10x
        out.peakrate{t}{p} = NaN(length(triggers),1);
        out.DSrate{t}{p} = NaN(length(triggers),1);
        out.normrate{t}{p} = NaN(length(triggers),1);
        
        %iterate over dspikes
        for r = 1:length(triggers)
            bins = win(1)+triggers(r):binsize:triggers(r)+win(2); %binsize
            out.psth{t}{p}(r,:) = histcounts(spiketimes,bins)/binsize; %binsize
            out.peakrate{t}{p}(r) = length(spiketimes(spiketimes>triggers(r) & spiketimes<triggers(r)+binsize))/binsize;
            out.DSrate{t}{p}(r) = length(spiketimes(spiketimes>d.starttime(r) & spiketimes<d.endtime(r)))/(d.endtime(r) - d.starttime(r));
            out.normrate{t}{p}(r) = (out.DSrate{t}{p}(r)-out.baserate{t}{p})/out.sdbaserate{t}{p}';
        end
    end
end