function out = calcSWRinstfreq(index, excludeperiods, eeg, ripples,chinfo, varargin)

%  Computes the histogram of low gamma freqs around the middle of each ripple event.

%  From maggie's calcinstantaneousfrequency  2/17/15 AG
%  Updated to calculate SWR inst freq, length, & number of cycles 5/4/17 EJ

%  Options:
%       appendindex-- Determines whether index is included in output vector
%           Default: 1
%       minstd-- min event size (in std)
%       bins-- for histogram
%       interripwin-- to eliminate ripples that occur too close together


%parse the options
appendindex = 1;
minstd = 5;    %can change minstd!!
interripwin = 0; %do not eliminate ripples that occur in doublets or triplets
win = [.4 .4]; %exlude rips too close to beginning or end

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'minstd'
                minstd = varargin{option+1};
            case 'interripwin'
                interripwin = varargin{option+1};
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
% trig = '(isequal($area,''ca3'') && contains($layer,''pyr''))';
trigindex = evaluatefilter(chinfo,trig);
out.trigindex = trigindex;
probeindex = evaluatefilter(chinfo,probe);
out.probeindex = probeindex;


%iterate through channels
if ~isempty(probeindex)
    load('ejmouseripplefilter.mat')
    % load('ej200to600hfofilter.mat')
    % ripplefilter = hfofilter;
    
    %get ripple trigger times
    r = ripples{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};
    Fs = r.samprate;
    valid = getvalidrips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);
    
    %get corresponding CA1 or CA3 rips
    if probeindex(1,3)~=trigindex(1,3)
        r1 = ripples{probeindex(1,1)}{probeindex(1,2)}{probeindex(1,3)};
        ripstarts = r1.starttime;
        ripends = r1.endtime;
        CA1rips = r.midtime(valid);
        validprobe = [];
        for b = 1:length(ripstarts) %iterate over all SWRs detected in CA3
            if find(isExcluded(CA1rips,[ripstarts(b) ripends(b)]))
                validprobe = [validprobe b]; %use SWR if it's also detected in CA1
            end
        end
        r = r1;
        valid = validprobe;
    end
    
    out.ripsizes = r.maxthresh(valid);
    startinds = r.startind(valid);
    endinds = r.endind(valid);
    
    out.rippeakamps = r.peak(valid);
    out.ripenergy = r.energy(valid);
    out.riplength = r.endind(valid) - r.startind(valid);   %in ms
    
    for p = 1:size(probeindex,1)
        cyclesperrip{p} = [];
        instfreqs{p} = [];
        avginstfreqs{p} = [];
        startindout{p} = [];
        %filter raw eeg for low gamma
        e = double(eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.data);
        ripples = filtfilt(ripplefilter.kernel,1,e);
        %don't need this unless you want to use phase to identify peaks
        % for ref, real(hgamma) = gamma
        %hgamma = hilbert(gamma);
        %phase = angle(hgamma);
        
        for r = 1:length(valid)  %iterate through events
            % MAKE SURE CHRONUX NOT ON PATH or wrong findpeaks function
            [pks, pkinds] = findpeaks(ripples(startinds(r):endinds(r)));
            instfreqs{p} = [instfreqs{p}; Fs./diff(pkinds)];
            avginstfreqs{p} = [avginstfreqs{p}; mean(Fs./diff(pkinds))];
            startindout{p} = [startindout{p};  startinds(r)];
            cyclesperrip{p} = [cyclesperrip{p}; length(pks)];
            
            %             %if mod(r,33)==0  %plot a random selection just to see how
            %             %things look
            %
            %             figure
            %             plot(real(ripples(startinds(r):endinds(r))))
            %             hold on
            %             scatter(pkinds, pks)
            %             %end
            
        end
    end
    out.instfreqs = instfreqs;
    out.avginstfreqs = avginstfreqs;
    out.startind = startindout;
    out.cyclesperrip = cyclesperrip;
else
    out.instfreqs = [];
    out.avginstfreqs = [];
    out.startind = [];
    out.cyclesperrip = [];
    out.riplength = [];
end

if (appendindex)
    out.index = index;
end