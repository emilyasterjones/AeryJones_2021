function out = contMUA(index, excludeperiods, mua, chinfo, varargin)

%parse the options
binsize = 0.1;
spikeFs = 30000;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'binsize'
                binsize = varargin{option+1};
            case 'probe'
                probecrit = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% search through chinfo struct to find channels meeting probe criteria
probeindex = evaluatefilter(chinfo,probecrit);
out.probeindex = probeindex;

if ~isempty(probeindex)    
    for p = 1:size(probeindex,1)
        pmua = mua{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)};
        spiketimes = pmua.spiketimes/spikeFs;
        bins = pmua.starttime:binsize:pmua.endtime;
        validinds = ~isExcluded(spiketimes,excludeperiods);
        
        %sum over time windows
        bininds = ~isExcluded(bins,excludeperiods);
        out.rate{p} = histcounts(spiketimes(validinds),bins(bininds))/binsize;
    end
end

out.index = index;