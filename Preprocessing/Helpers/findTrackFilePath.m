% Returns a track recording file of any type (e.g. '.videoPositionTracking')
% Helpful for when you may have multiple animals run at once (multtracks=1)
% or multiple epochs per day (multepochs=1)
% Can even handle mixed tracks: sometimes multiple animals at once,
% sometimes not (multtracks=2)
% E.g. running 2 animals at once on a linear track, 1st epoch of the day
% File names and folder name are Animal1_Animal2_track-session-epoch
% Used in get_tracking_start.m, extract_pos.m, extract[LinearTrack|WMaze]Events.m,
% and calc[LinearTrack|WMaze]Performance.m

function filename = findTrackFilePath(datadir, animaldir, prefix, session, trackepoch, filetype, multtracks, multepochs)

if multepochs
    index = sprintf('%02d-%02d',session,trackepoch);
else
    index = sprintf('%02d',session);
end

%animal names are name_task (eg Cardamom_OCC)
suffind = strfind(prefix,'_');
animal = prefix(1:suffind-1);
suffix = prefix(suffind+1:end);

%left vs right for OCC depends on what's in the task struct
if strcmpi(suffix,'OCC')
    taskname = sprintf('%s\%stask.mat',animaldir,prefix);
    load(taskname);
    if strcmpi(task{session}{trackepoch}.track,'A') %context A, left
        videoprefix = '.1';
    else %context B, right
        videoprefix = '.2';
    end
end

%unsure if multiple animals run at once, so try both options
if multtracks==2 
    filemask = sprintf('%s%s',prefix,index);
    findfolder = dir([datadir '\' filemask]);
    if ~isempty(findfolder) %single animal
        multtracks=0;
    else
        multtracks=1;
    end
end

if multtracks
    filemask = sprintf('%s*_%s%s',animal,suffix,index);
    findfolder = dir([datadir '\' filemask]);
    if ~isempty(findfolder) %animal is first in pair
        if strcmpi(filetype,'.videoPositionTracking') %2 cameras, .1 for left & .2 for right
            if ~strcmpi(suffix,'OCC')
                videoprefix = '.1';
            end
            filetype = [videoprefix filetype];
        end
    else %animal is second in pair
        filemask = sprintf('*%s_%s%s',animal,suffix,index);
        if strcmpi(filetype,'.videoPositionTracking') %2 cameras, .1 for left & .2 for right
            if ~strcmpi(suffix,'OCC')
                videoprefix = '.2';
            end
            filetype = [videoprefix filetype];
        end
    end
    
    findfile = dir([datadir '\' filemask '\' filemask filetype]);
    filename = [findfile.folder '\' findfile.name];
else
    filemask = sprintf('%s%s',prefix,index);
    filemaskpath = [datadir '\' filemask '\' filemask];
    if strcmpi(filetype,'.videoPositionTracking') %2 cameras, .1 for left & .2 for right
        if ~strcmpi(suffix,'OCC')
            videoprefix = '.1';
        end
        filetype = [videoprefix filetype];
    end
    filename = [filemaskpath filetype];
end