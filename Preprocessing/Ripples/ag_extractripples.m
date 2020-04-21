function ag_extractripples(animaldir, prefix, day, min_suprathresh_duration, nstd, varargin)

%function extractripples(animaldir, prefix, day,
%min_suprathresh_duration, nstd, options)
%
%	Reads in the ripple files from the specified day and tetrodes and
%	extracts all of the ripples from that tetrodes.
%
%	assumes position data stores in pos file in animdirectory
%
%directoryname - example '/data99/user/animaldatafolder/', a folder
%                containing processed matlab data for the animal
%
%fileprefix	- folder name where the day's data is stored
%
%day		- the day to process
%
%min_suprathresh_duration
%		- the time (in seconds) which the signal
%       must remain above threshold to be counted as as ripple; this guards
%       against short spikes in signal (e.g. noise) as being counted as
%       ripples. Set min_suprathreshold_duration to some small value, like
%       0.015 s.
%
%nstd		- the number of standard dev that ripple must be from mean to
%			be detected. Start with 2.
%
%
%options	'stdev', stdev   sets the size of the standard deviation used to
%				allow direct comparison across epochs
%       	'baseline', b   sets the size of the baseline used to
%				allow direct comparison across epochs
%           'maxpeakval, m	- ripples with maximal peaks above this value
%				are exluded.  Use this avoid detecting noise
%				events. Default 1000
%           'samethreshperday' - 0 or 1, default 0,
%               0 calculates baseline and threshold per epoch for each
%               tetrode
%               1 uses same threshold for all sessions on the same day on
%               the same tetrode.  Calculates this threshold based on the
%               baseline and stdev for all sessions (run and sleep) for
%               the entire day.    AS added 1-12-10
%
% Outputs:
%ripples 	- structue with various fields, including the following which
%			describe each ripple.
%	starttime - time of beginning of ripple
%	endtime	  - time of end of ripple
%	midtime   - time of midpoint of energy of event
%	peak	  - peak height of waveform)
%	maxthresh - the largest threshold in stdev units at which this ripple
%			would still be detected.
%	energy	  - total sum squared energy of waveform
%	startind  - index of start time in ripple structure
%	endind    - index of end time in ripple structure
%	midind    - index of middle time in ripple structure
%	posind    - index into pos structure for the midpoint of this ripple
% 	posinterp - interpolation factor value between adjacent position
% 		    elements
%
%

stdev = 0;
baseline = 0;
samethreshperday = 0;
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'stdev'
                stdev = varargin{option+1};
            case 'baseline'
                baseline = varargin{option+1};
            case 'samethreshperday'
                samethreshperday = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

% define the standard deviation for the Gaussian smoother which we
% apply before thresholding (this reduces sensitivity to spurious
% flucutations in the ripple envelope)
smoothing_width = 0.004; % 4 ms

taskname = sprintf('%s/%stask.mat',animaldir,prefix);
load(taskname);
for epoch = 1:length(task{day})
    if(~contains(lower(task{day}{epoch}.descript), '*fail*') && ...
            ~contains(lower(task{day}{epoch}.descript), '*not run*'))
        % go through each tetrode
        tmpflist = dir(sprintf('%s/%sripple%02d-%d-*.mat', animaldir, prefix, day, epoch));
        
        %AS added following 1-12-10 to allow all epochs to use same threshold
        if samethreshperday == 1
            allrenv = [];
            for i = 1:length(tmpflist)
                load(tmpflist(i).name);
                % get the epoch number
                dash = find(tmpflist(i).name == '-');
                endidx = length(dash);
                t = str2num(tmpflist(i).name((dash(endidx)+1):(dash(endidx)+3)));
                
                % convert the ripple envelope field to double
                temprenv = double(ripple{day}{epoch}{t}.envmag);
                allrenv = [allrenv; temprenv];
            end
            %baseline = mean(allrenv);
            %stdev = std(allrenv);
            thresh = baseline + nstd * stdev;
        end
        
        
        for i = 1:length(tmpflist)
            % load the ripple file
            load([tmpflist(i).folder '\' tmpflist(i).name]);
            % get the epoch number
            dash = find(tmpflist(i).name == '-');
            endidx = length(dash);
            t = str2num(tmpflist(i).name((dash(endidx)+1):(dash(endidx)+3)));
            
            % convert the ripple envelope field to double
            renv = double(ripple{day}{epoch}{t}.envmag);
            
            % smooth the envelope:
            samprate = ripple{day}{epoch}{t}.samprate;
            kernel = gaussian(smoothing_width*samprate, ceil(8*smoothing_width*samprate));
            renv = smoothvect(renv, kernel);
            % find the ripples
            % calculate the duration in terms of samples
            mindur = round(min_suprathresh_duration * samprate);
            
            % calculate the threshold in uV units
            if samethreshperday == 0 % is calculating threshold per epoch
                baseline = mean(renv);
                stdev = std(renv);
                thresh = baseline + nstd * stdev;
            end
            
            % extract the events if this is a valid trace
            if (thresh > 0) && any(find(renv<baseline))
                tmprip = extracteventsnew(renv, thresh, baseline, 0, mindur, 0)';
                %AG MK changed 10/8/13 new extract events also returns time and ind of
                %thresh crossing
                
                % Assign the fields
                % start and end indeces
                rip.startind = tmprip(:,1);
                rip.endind = tmprip(:,2);
                % middle of energy index
                rip.midind = tmprip(:,8);
                
                %convert the samples to times for the first three fields
                rip.starttime = ripple{day}{epoch}{t}.starttime + rip.startind / samprate;
                rip.endtime = ripple{day}{epoch}{t}.starttime + rip.endind / samprate;
                rip.midtime = ripple{day}{epoch}{t}.starttime + rip.midind / samprate;
                rip.peak = tmprip(:,3);
                rip.energy = tmprip(:,7);
                rip.maxthresh = (tmprip(:,9) - baseline) / stdev;
                rip.threshind = tmprip(:,10);
                rip.threshtime = ripple{day}{epoch}{t}.starttime + tmprip(:,10) / samprate;
                
            else
                rip.startind = [];
                rip.endind = [];
                rip.midind = [];
                rip.starttime = [];
                rip.endtime = [];
                rip.midtime = [];
                rip.peak = [];
                rip.energy = [];
                %%rip.posind = [];
            end
            
            rip.timerange = [0 length(renv)/samprate] + ripple{day}{epoch}{t}.starttime;
            rip.samprate = ripple{day}{epoch}{t}.samprate;
            rip.threshold = thresh;
            rip.baseline = baseline;
            rip.std = stdev;
            rip.minimum_duration = min_suprathresh_duration;
            
            %ripple
            clear feeg
            ripples{day}{epoch}{t} = rip;
            clear rip;
        end
    end
end

if exist('ripples')
    save(sprintf('%s/%sripples%02d.mat', animaldir, prefix, day), 'ripples');
end
end