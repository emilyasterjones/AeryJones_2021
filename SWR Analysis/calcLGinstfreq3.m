function out = calcLGinstfreq3(index, excludeperiods, eeg, ripples,chinfo, varargin)

%  Computes the histogram of low gamma freqs around the middle of each ripple event.

%  From maggie's calcinstantaneousfrequency  2/17/15 AG

%  Options:
%       appendindex-- Determines whether index is included in output vector
%           Default: 1
%       minstd-- min event size (in std)
%       bins-- for histogram
%       interripwin-- to eliminate ripples that occur too close together



%parse the options
appendindex = 1;
minstd = 5;    %can change minstd!!
interripwin = 0;%do not eliminate ripples that occur in doublets or triplets
win = [.4 .4]; %exlude rips too close to beginning or end
Fs = 1000;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'minstd'
                minstd = varargin{option+1};
            case 'bins'
                bins = varargin{option+1};
            case 'interripwin'
                cwin = varargin{option+1};
            case 'probe'
                probe = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% search through chinfo struct to find trig and probe chans
trig = '(isequal($area,''ca1'') && contains($layer,''pyr 1''))'; %
trigindex = evaluatefilter(chinfo,trig);
out.trigindex = trigindex;
probeindex = evaluatefilter(chinfo,probe);
out.probeindex = probeindex;

%get ripple trigger times
rip = ripples{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};
[valid, excluded] = getvalidrips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);
out.ripsizes = rip.maxthresh(valid);
startinds = rip.startind(valid);
endinds = rip.endind(valid);

load('aglowgamfilter.mat') 

peakdist = [];
cyclesperrip = [];

%iterate through channels
if ~isempty(probeindex)
    for p = 1:size(probeindex,1)
        %peakdist{p} = [];
        cyclesperrip{p} = [];
        instfreqs{p} = nan(length(valid),20);
        avginstfreqs{p} = [];
        startindout{p} = [];
        %filter raw eeg for low gamma
        e = double(eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.data);
        gamma = filtfilt(lowgamfilter.kernel,1,e);
        
        for r = 1:length(valid)  %iterate through events
            % MAKE SURE CHRONUX NOT ON PATH or wrong findpeaks function
            [pks, pkinds] = findpeaks(gamma(startinds(r):endinds(r)));
            if length(pkinds)>1 && length(pkinds)<=20 %can't pad an empty array
                instfreqs{p}(r,:) = padarray(Fs./diff(pkinds),21-length(pkinds),NaN,'post')'; 
            end
            cyclesperrip{p} = [cyclesperrip{p}; length(pkinds)];
            avginstfreqs{p} = [avginstfreqs{p}; mean(Fs./diff(pkinds))];
            startindout{p} = [startindout{p};  startinds(r)];
        end
    end
    out.instfreqs = instfreqs;
    out.cyclesperrip = cyclesperrip;
    out.avginstfreqs = avginstfreqs;
    out.startind = startindout;
else
    out.instfreqs = [];
    out.cyclesperrip = [];
    out.avginstfreqs = [];
    out.startind = [];
end

if (appendindex)
    out.index = index;
end