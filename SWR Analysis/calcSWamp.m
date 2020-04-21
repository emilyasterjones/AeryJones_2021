function out = calcSWamp(index, excludeperiods, eeg, ripples, chinfo, varargin)
% extracts highest amplitude of 0.2-30 Hz filtered trace within 200ms around ripple

% Outputs:
% out.peak - peak amplitude of the sharp wave filtered trace around each ripple

% defaults
win = [0.2 0.2];
minstd = 3;
interripwin = 0;
Fs = 1000;
for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'win'
                win = varargin{option+1};
            case 'minstd'
                minstd = varargin{option+1};
            case 'interripwin'
                interripwin = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

load('ejsharpwavefilter2.mat');

location = '(isequal($area,''ca1'')&& contains($layer,''pyr 1''))';%
pyr1chan = evaluatefilter(chinfo,location);
r = ripples{index.epochs(1)}{index.epochs(2)}{pyr1chan(1,3)};
trigindex = [index.epochs(1) index.epochs(2) pyr1chan(1,3)];
[valid, ~] = getvalidrips(ripples,index,trigindex, excludeperiods, win, minstd, interripwin);
r.startind = r.startind(valid);
r.endind = r.endind(valid);

for c = 1:length(index.chinfo)
    eegs = eeg{index.epochs(1)}{index.epochs(2)}{index.chinfo(c)};
    
    % filter it
    filtered_data = filtereeg2(eegs, sharpwavefilter, 'int16', 1);
    renv = double(filtered_data.data(:,3));
    
    out.peak{c} = [];
    for t = 1:length(r.startind)
        sw = renv([r.startind(t) - win(1)*Fs, r.endind(t) + win(2)*Fs]);
        out.peak{c} = [out.peak{c} max(abs(sw))];
    end
end

out.index = index;
end