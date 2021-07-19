function out = RTCSD(index, excludeperiods, eeg, ripples, chinfo, varargin)
%parse the options
win = [-0.2 0.2];
minstd = 5;
interripwin = 0;  %does not eliminate ripples that occur in doublets or triplets
Fs = 1000;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'win'
                win = varargin{option+1};
            case 'minstd'
                minstd = varargin{option+1};
            case 'trig'
                trigcrit = varargin{option+1};
            case 'freqband'
                freqs = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

out.win = win;

% search through chinfo struct to find channels meeting trig criteria
trigindex = evaluatefilter(chinfo,trigcrit);
out.trigindex = trigindex;

%get EEG & rest times
epoch_eeg = eeg{trigindex(1,1)}{trigindex(1,2)};
Fs_temp = length(epoch_eeg{1}.data)/(epoch_eeg{1}.endtime-epoch_eeg{1}.starttime);
fulltimes = epoch_eeg{1}.starttime:1/Fs_temp:epoch_eeg{1}.endtime;
validtimes = ~isExcluded(fulltimes, excludeperiods);

%if freq band specified, filter trace
if exist('freqs','var')
    for i=1:length(epoch_eeg)
        coeff = designeegfilt(Fs,freqs(1),freqs(2));
        data = double(epoch_eeg{i}.data);
        epoch_eeg{i}.data = int16(filtfilt(coeff,1,data));
    end
end

%smoothing via triangular kernel
for i = 2:length(epoch_eeg)-1
    %smoothing via triangular kernel
    b = epoch_eeg{i}.data;
    a = epoch_eeg{i-1}.data;
    c = epoch_eeg{i+1}.data;
    smooth_eeg{i}.data = (a+2*b+c)/4;
end
%second spatial derivative
csd = NaN(length(epoch_eeg),length(epoch_eeg{1}.data)); %channels x time
for i = 3:length(epoch_eeg)-2
    %smoothing via triangular kernel
    b = smooth_eeg{i}.data;
    a = smooth_eeg{i-1}.data;
    c = smooth_eeg{i+1}.data;
    csd(i,:) = (a-2*b+c)/4;
end
out.csd = csd;

%mean & SD over rest trace
csd_base = csd(:,validtimes(1:length(csd)));
out.baseline = mean(csd_base,2);
out.sd = std(csd_base,0,2);

%get ripple trigger times
rip = ripples{trigindex(1,1)}{trigindex(1,2)}{trigindex(1,3)};
[valid, ~] = getvalidrips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);
rip.starttime = rip.starttime(valid);
rip.endtime = rip.endtime(valid);
triggers = rip.threshind(valid);
out.ripinds = triggers;

%bin by ripple triggers
ripcsd = NaN(length(triggers),length(epoch_eeg),length(win(1):1/Fs:win(2))); %rips x channels x time
for r = 1:length(triggers)
    ripcsd(r,:,:) = csd(:,win(1)*Fs+triggers(r):triggers(r)+win(2)*Fs);
end
out.ripcsd = ripcsd;

%calculate peaks
out.peakripsource = max(ripcsd,[],3); %rips x channels
out.peakripsink = min(ripcsd,[],3); %rips x channels

%Z-score
out.zripsource = (out.peakripsource - out.baseline') ./ out.sd';
out.zripsink = (out.peakripsink - out.baseline') ./ out.sd';

end