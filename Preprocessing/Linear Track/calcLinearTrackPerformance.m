function calcLinearTrackPerformance(datadir, animaldir, prefix, session, varargin)

% creates a perf structure based on a state script log for a linear track run

trackepoch = 1;
multepochs = 0;
multtracks = 0;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'trackepoch' %if track epoch needs to be specified
                trackepoch = varargin{option+1};
                multepochs = 1;
            case 'multtracks' %if running 2 tracks on 1 ECU
                multtracks = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

load(sprintf('%s/%stask.mat',animaldir, prefix));
if(~contains(lower(task{session}{trackepoch}.descript), 'fail') && ...
        ~contains(lower(task{session}{trackepoch}.descript), 'not run'))
    
    track = task{session}{trackepoch}.track;
    rewards = 0;
    filename = findTrackFilePath(datadir, animaldir, prefix, session, trackepoch, '.stateScriptLog', multtracks, multepochs);
    events = parseStateScriptFile(filename);
    
    %remove all updates that aren't rewards
    e = 1;
    while e <= length(events)
        if isempty(strfind(events(e).message,'Reward'))
            events(e) = [];
        else
            e = e+1; %update the pointer only if the row wasn't deleted
        end
    end
    if multtracks
        if contains(track, 'left')
            maze = 'Left maze';
        else
            maze = 'Right maze';
        end
        e = 1;
        while e <= length(events)
            if isempty(strfind(events(e).message,maze))
                events(e) = [];
            else
                e = e+1;
            end
        end
    end
    
    pokes = length(events);
    timepertrial = 30*60/pokes;
    if pokes >= 1
        timetofirstrew = (events(1).time - task{session}{trackepoch}.start)/1000;
    else
        timetofirstrew = NaN;
    end
    if timetofirstrew < 0 %mouse poked before tracking starts (still removing hand from camera view)
        timetofirstrew = 0;
    end
    if pokes >= 2
        timetofirstalt = (events(2).time - events(1).time)/1000;
    else
        timetofirstalt = NaN;
    end
    
    perf{session}{trackepoch}.pokes = pokes;
    perf{session}{trackepoch}.timepertrial = timepertrial;
    perf{session}{trackepoch}.timetofirstrew = timetofirstrew;
    perf{session}{trackepoch}.timetofirstalt = timetofirstalt;
    save(sprintf('%s/%sperf%02d.mat', animaldir, prefix, session), 'perf');
    
end