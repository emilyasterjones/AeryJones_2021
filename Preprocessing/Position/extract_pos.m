function extract_pos(datadir, animaldir, prefix, day, varargin)

%% load task file and trodes data
%load task file
track = '';
orientation = '';
pix2cm = 40/535;
multtracks = 0;
trackepoch = 0;

taskname = sprintf('%s/%stask.mat',animaldir,prefix);
load(taskname);
nepochs = length(task{day});

%added EJ 1/12/18 to accomodate track runs
for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'trackepoch'
                trackepoch = varargin{option+1};
            case 'track'
                if trackepoch==0
                    error('Must provide trackepoch')
                end
                track = varargin{option+1};
                if strcmpi(track, 'LT')
                    orientation = 'horizontal';
                end
            case 'multtracks' %if running 2 tracks on 1 ECU
                multtracks = varargin{option+1};
            case 'orientation'
                orientation = varargin{option+1};
                if strcmpi(track, 'LT')
                    orientOpts = {'horizontal', 'vertical'};
                    if ~sum(strcmpi(orientation, orientOpts)) %same as strcmpi || strcmpi
                        error('Orientation must be horizontal or vertical')
                    end
                else
                    error('Orientation not recognized unless track specified')
                end
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

%Velocity is calculated from smoothed tracking data to smooth over jitter
%Gaussian filter averages over the 1s before & after 30Hz sampled data
posfilt = gaussian(30*0.5, 60);
multepochs = (nepochs>1);

for epoch = 1:nepochs
    if(~contains(lower(task{day}{epoch}.descript), '*fail*') && ...
            ~contains(lower(task{day}{epoch}.descript), '*not run*'))
        
        %load trodes pos file
        if epoch==trackepoch
            filename = findTrackFilePath(datadir, animaldir, prefix, day, epoch, '.videoPositionTracking', multtracks, multepochs);
        else
            filename = findTrackFilePath(datadir, animaldir, prefix, day, epoch, '.videoPositionTracking', 0, multepochs);
        end
        filename = sprintf('%s/%s%02d.videoPositionTracking',datadir,prefix,day);
        posraw = readTrodesExtractedDataFile(filename);
        starttemp = task{day}{epoch}.start/1000;
        endtemp = task{day}{epoch}.end/1000;
        posraw.clockrate = 30000;
        
        startidx = lookup(starttemp,double(posraw.fields(1).data)/posraw.clockrate);
        endidx = lookup(endtemp,double(posraw.fields(1).data)/posraw.clockrate);
        %store variables and convert from pixels to cm   .data = [time, x, y]
        rawpos{day}{epoch}.data = [double(posraw.fields(1).data(startidx:endidx))/posraw.clockrate, double(posraw.fields(2).data(startidx:endidx))*pix2cm, double(posraw.fields(3).data(startidx:endidx))*pix2cm];
        rawpos{day}{epoch}.fields = 'time xpos ypos';
        pos{day}{epoch} = ag_addvelocity(rawpos{day}{epoch}, posfilt); %calculate velocity and smooth
        
        % convert to DF2.0 format (convert data to individual fields)
        pos{day}{epoch}.time = pos{day}{epoch}.data(:,1);
        pos{day}{epoch}.x = pos{day}{epoch}.data(:,2);
        pos{day}{epoch}.y = pos{day}{epoch}.data(:,3);
        pos{day}{epoch}.vel = pos{day}{epoch}.data(:,4);
        
        
        %added EJ 1/8/18
        %if running on a track, calculates the distance to each port (L/R)
        if epoch==trackepoch
            if strcmpi(track, 'LT')
                if strcmpi(orientation, 'horizontal')
                    leftportpos = min(pos{day}{epoch}.x);
                    rightportpos = max(pos{day}{epoch}.x);
                    pos{day}{epoch}.disttoleft = abs(pos{day}{epoch}.x - leftportpos);
                    pos{day}{epoch}.disttoright = abs(pos{day}{epoch}.x - rightportpos);
                elseif strcmpi(orientation, 'vertical')
                    leftportpos = min(pos{day}{epoch}.y);
                    rightportpos = max(pos{day}{epoch}.y);
                    pos{day}{epoch}.disttoleft = abs(pos{day}{epoch}.y - leftportpos);
                    pos{day}{epoch}.disttoright = abs(pos{day}{epoch}.y - rightportpos);
                else
                    error('Orientation must be specified for tracks.')
                end
            end
        end
        
        pos{day}{epoch} = rmfield(pos{day}{epoch},'data');
        pos{day}{epoch} = rmfield(pos{day}{epoch},'fields');
    end
end

save(sprintf('%s/%spos%02d.mat', animaldir, prefix, day), 'pos');

end