function getTimeStamps(path, varargin)
%read in trodes timestamps files & displays start and end time
%inputs: directory containing .videoTimeStamps files for this experiment
%requires readCameraModuleTimeStamps.m

datestr = '1-Jan-1970';
suffix = '/*.1.videoTimeStamps';

%specify a subset of a folder to give timestamps for, such as only
%recordings after a certain date or only recordings from later versions of
%Trodes
for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'startdate'
                datestr = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

folders = dir(path);
folders(~[folders.isdir]) = []; %remove non-folders
for f = 1:length(folders)
    if(length(folders(f).name)>2) %remove '.' and '..' from dir list
        files = dir(strcat(path,'\',folders(f).name,suffix)); %file list
        cd(strcat(path,'\',folders(f).name))
        for i = 1:length(files) %print start & end times in order by filename
            if(files(i).datenum > datenum(datestr))
                ts = readCameraModuleTimeStamps(files(i).name);
                if(numel(ts)~=0)
                    fprintf('%s\t%f\t%f\n', files(i).name, ts(1), ts(end));
                end
            end
        end
    end
end

end