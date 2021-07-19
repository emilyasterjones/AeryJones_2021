function out = DSdescript(index, excludeperiods, dspikes, eeg, chinfo, varargin)
% out = ripdescript(index, excludeperiods, ripples, options)
%  Computes basic characteristics of detected ripple events

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'trig'
                trigcrit = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

%open dspikes
trigindex = evaluatefilter(chinfo,trigcrit);
out.trigindex = trigindex;
d = dspikes{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};
e = eeg{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};
valid = find(~isExcluded(d.starttime,excludeperiods) & ~isExcluded(d.endtime,excludeperiods));

%calculate DS rate in 1min bins
Fs = e.samprate;
fulltimes = e.starttime:1/Fs:e.endtime;
validtimeinds = find(~isExcluded(fulltimes, excludeperiods));
timeseq = [validtimeinds(1:60000:end); validtimeinds(end)];
timecounts = histcounts(d.startind(valid),timeseq)/60; %rips per 1min bin

out.dspeakamps = d.peak(valid);
out.dslength = d.endind(valid) - d.startind(valid);   %in ms
out.timecounts = timecounts;