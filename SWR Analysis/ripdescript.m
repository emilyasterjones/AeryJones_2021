function out = ripdescript(index, excludeperiods, ripples, chinfo, varargin)
% out = ripdescript(index, excludeperiods, ripples, options)
%  Computes basic characteristics of detected ripple events

interripwin = 0;  %eliminate ripples that occur in doublets or triplets
minstd = 5;
win = [.4 .4];   %to exclude rips too close to beginning or end

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'interripwin'
                win = varargin{option+1};
            case 'minstd'
                minstd = varargin{option+1};
            case 'probe'
                probecrit = varargin{option+1};
            case 'trig'
                trigcrit = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% only works for ONE channel (pyr1)

%get valid ripples
trigindex = evaluatefilter(chinfo,trigcrit);
out.trigindex = trigindex;
r = ripples{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};
[valid, excluded, validdur] = getvalidrips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);

%if looking at SWR properties in CA3 channels
if exist('probecrit','var')
    probeindex = evaluatefilter(chinfo,probecrit);
    if ~isempty(probeindex)
        out.probeindex = probeindex;
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
        valid = validprobe;
        r = r1;
    else %if no CA3 channel, return empty 
        out.excluded = NaN;
        out.ripstds = NaN;
        out.rippeakamps = NaN;
        out.ripenergy = NaN;
        out.riplength = NaN;
        out.riplong = NaN;
        out.interriptimes = NaN;
        out.ripchains = NaN;
        out.validdur = NaN;
        return
    end
end

%calculate window to next event
interriptimes = zeros(length(r.starttime(valid))-1,1);
for i = 1:length(r.starttime(valid))-1
    rip_end = r.endtime(valid(i));
    next_rip = r.starttime(valid(i+1));
    %find the next exclusion period
    time_to_excl = excludeperiods(:,1)-rip_end;
    past_excl_count = length(time_to_excl) - sum(time_to_excl>0);
    [~, excl_ind] = min(time_to_excl(time_to_excl>0));
    excl_ind = excl_ind + past_excl_count;
    if isempty(time_to_excl(time_to_excl>0)) || next_rip < excludeperiods(excl_ind,1)
        interriptimes(i) = next_rip-rip_end;
    else
        interriptimes(i) = NaN; %no interriptime for rips right before exclusion period
    end
end
interriptimes(end+1) = NaN; %no interriptime for last recorded ripple

%count # doublets, triplets, etc
ripchains = zeros(1,50);
counter = 0;
for i=1:length(interriptimes)
    if interriptimes(i)<0.2
        counter = counter+1; %add a rip to the chain
    elseif counter>0
        ripchains(counter) = ripchains(counter)+1; %end of chain
        counter = 0;
    end
end

out.excluded = excluded;
out.ripstds = r.maxthresh(valid);
out.rippeakamps = r.peak(valid);
out.ripenergy = r.energy(valid);
out.riplength = r.endind(valid) - r.startind(valid);   %in ms
out.riplong = sum(out.riplength>=100)/length(r.startind);
out.interriptimes = interriptimes;
out.ripchains = ripchains;
out.validdur = validdur;