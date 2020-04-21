function [mostvalid, excluded, validdur] = getvalidrips(ripples,index,trigindex,excludeperiods, win, minstd, interripwin)

% getvalidrips runs within DF2.0
% given ripples for appopriate day and trigindex,
% Mostvalid criteria:
% - larger than threshold SD
% - not too close to beginning or end of recording
% - do not occur during exclude times
% - are not detected on more than 12 channels (likely noise)
% IF interripwin is provided, also excludes rips that happen within 1sec of another ripple (endtime - starttime)

% Returns indices of ripples that meet criteria
% Returns exclusion stats [total, noise events, overlapping events]

%get ripples for trig chan

% EJ 4/6/17 changed to eliminate rips that start or end during
% excludeperiods (instead of just overlapping with midpoint)
r = ripples{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};
Fs = r.samprate;
fulltimes = r.timerange(1):1/Fs:r.timerange(2);
validtimeinds = find(~isExcluded(fulltimes, excludeperiods));
validtimes = fulltimes(validtimeinds);
validdur = length(validtimes)/Fs;
validripples = find((r.maxthresh >= minstd) & ~isExcluded(r.starttime,excludeperiods) & ~isExcluded(r.endtime,excludeperiods));
excluded(1) = length(validripples); %initial number of detected nonexcluded rips

if isempty(validripples) %no valid ripples, so don't look for additional exclusions
    mostvalid = validripples;
    return
end

% try  %will run as long as no error thrown
    %exclude ripples too close to beginning or end
    while r.starttime(validripples(1)) < r.timerange(1)+win(1) %r.threshtime
        validripples(1) = [];
    end
    
    while r.endtime(validripples(end)) > r.timerange(2)-win(2) %r.threshtime
        validripples(end) = [];
    end
    
    ripstarts = r.starttime(validripples);
    ripends = r.endtime(validripples);
    counter = zeros(length(ripstarts),length(ripples{index.epochs(1)}{index.epochs(2)}));
    for b = 1:length(ripstarts)
        for c = 1:length(ripples{index.epochs(1)}{index.epochs(2)})
            r1 = ripples{index.epochs(1)}{index.epochs(2)}{c};
            if ~isempty(r1.startind)
                % EJ 4/6/17 changed to require noise events to occur at 5SD for detection
                vr1 = find(r1.maxthresh*r1.std+r1.baseline > r.std*5+r.baseline);%1:length(r1.midtime);
                if find(isExcluded(r1.midtime(vr1),[ripstarts(b) ripends(b)]))
                    counter(b,c) = 1;
                end
            end
        end
    end
    counter = sum(counter,2);
    % EJ 4/6/17 changed to require 16 or more channels to overlap
    morevalid = validripples(counter<16);
    noiserips = length(validripples)-length(morevalid);
    excluded(2) = noiserips;
    
    
%     if interrip win provided, exclude events that happen within 1sec of another event
    if (nargin>6)
        intervals = r.starttime(morevalid(2:end))-r.starttime(morevalid(1:end-1));
        for i = 1:length(intervals)
            if intervals(i)<interripwin
                morevalid(i) = 0;
                morevalid(i+1) = 0;
            end
        end
    end
    mostvalid = morevalid(morevalid>0);
    overlaprips = length(morevalid)-length(mostvalid);
    excluded(3) = overlaprips;

end