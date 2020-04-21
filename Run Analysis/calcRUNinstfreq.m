function out = calcRUNinstfreq(index, excludeperiods, eeg, chinfo, varargin)

%  Computes the instantaneous frequency of filtered bands across continuous data
% EJ 10/14/2019

%parse the options
appendindex = 1;
Fs = 1000;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'probe'
                probe = varargin{option+1};
            case 'thetafilt'
                thetafilt = varargin{option+1};
            case 'SGfilt'
                SGfilt = varargin{option+1};
            case 'FGfilt'
                FGfilt = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

%design default filters if not specified
if isempty(thetafilt)
    theta = [2 12];
    thetafilt = designeegfilt(1000,theta(1),theta(2));
    SG = [20 65];
    SGfilt = designeegfilt(1000,SG(1),SG(2));
    FG = [50 120];
    FGfilt = designeegfilt(1000,FG(1),FG(2));
end

% search through chinfo struct to find trig and probe chans
probeindex = evaluatefilter(chinfo,probe);
out.probeindex = probeindex;
eprobe = eeg{probeindex(1,1)}{probeindex(1,2)}{probeindex(1,3)};
dur = eprobe.endtime - eprobe.starttime;
Fs = length(eprobe.data)/dur;

%iterate through channels
if ~isempty(probeindex)
    for p = 1:size(probeindex,1)
        thetainstfreq{p} = [];
        thetaavginstfreq{p} = [];
        SGinstfreq{p} = [];
        SGavginstfreq{p} = [];
        FGinstfreq{p} = [];
        FGavginstfreq{p} = [];
        
        %filter raw eeg for freq bands
        e = double(eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.data);
        thetaeeg = filtfilt(thetafilt,1,e);
        SGeeg = filtfilt(SGfilt,1,e);
        FGeeg = filtfilt(FGfilt,1,e);
        
        %exclude exclusion periods
        fulltimes = eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.starttime:1/Fs:eeg{probeindex(p,1)}{probeindex(p,2)}{probeindex(p,3)}.endtime;
        validtimes = ~isExcluded(fulltimes,excludeperiods);
        validtimes = validtimes(1:length(eprobe.data));
        thetaeeg = thetaeeg(validtimes);
        SGeeg = SGeeg(validtimes);
        FGeeg = FGeeg(validtimes);
        
        % MAKE SURE CHRONUX NOT ON PATH or wrong findpeaks function
        [~, pkinds] = findpeaks(thetaeeg);
        thetainstfreq{p} = Fs./diff(pkinds);
        thetaavginstfreq{p} = mean(Fs./diff(pkinds));
        [~, pkinds] = findpeaks(SGeeg);
        SGinstfreq{p} = Fs./diff(pkinds);
        SGavginstfreq{p} = mean(Fs./diff(pkinds));
        [~, pkinds] = findpeaks(FGeeg);
        FGinstfreq{p} = Fs./diff(pkinds);
        FGavginstfreq{p} = mean(Fs./diff(pkinds));
    end
    
    out.thetainstfreq = thetainstfreq;
    out.thetaavginstfreq = thetaavginstfreq;
    out.SGinstfreq = SGinstfreq;
    out.SGavginstfreq = SGavginstfreq;
    out.FGinstfreq = FGinstfreq;
    out.FGavginstfreq = FGavginstfreq;
else
    out.thetainstfreq = [];
    out.thetaavginstfreq = [];
    out.SGinstfreq = [];
    out.SGavginstfreq = [];
    out.FGinstfreq = [];
    out.FGavginstfreq = [];
end

if (appendindex)
    out.index = index;
end