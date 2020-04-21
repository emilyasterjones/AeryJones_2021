function out = velquant(index, excludeperiods, pos, varargin)
% calculates how much time is spent in various velocity bins vs total time

bins = [0 1 3 5 15];
appendindex = 1;

for option = 1:2:length(varargin)-1   
    if ischar(varargin{option})       
        switch(varargin{option})
            case 'appendindex'
                appendindex = varargin{option+1};
            case 'bins'
                bins = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end        
    else
        error('Options must be strings, followed by the variable');
    end
end

%count distribution of velocities, measured each frame
temp = histc(pos{index.epochs(1)}{index.epochs(2)}.vel,bins);
out.avgspeed = nanmean(pos{index.epochs(1)}{index.epochs(2)}.vel(pos{index.epochs(1)}{index.epochs(2)}.vel>1));
out.maxspeed = max(pos{index.epochs(1)}{index.epochs(2)}.vel);
%convert to seconds
%UPDATED 2/2/2019 to account for different frame rates (eg from timestamps bug)
fps = (pos{index.epochs(1)}{index.epochs(2)}.time(end)-pos{index.epochs(1)}{index.epochs(2)}.time(1))/...
        length(pos{index.epochs(1)}{index.epochs(2)}.vel);
counts = temp(1:length(bins)-1)*fps;
out.rawcounts = counts;
out.totaltime = pos{index.epochs(1)}{index.epochs(2)}.time(end)-pos{index.epochs(1)}{index.epochs(2)}.time(1);
out.relcounts = counts/out.totaltime;  %in fraction of total time

if (appendindex)
    out.index = index;
end