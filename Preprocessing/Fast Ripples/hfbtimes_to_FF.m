% Reads start and end times saved as .txt files from Igor pipline and
% stores them in DF2.0 format

datadir = 'E:\pHFO times';
filterstring = '(isequal($area,''ca1'') && contains($layer,''pyr 1''))';
animalbasedir = 'E:\CRCNS\Cohort 2';

epochfilter = [];
datafilter = [];
datafilter{1} = {'chinfo',filterstring};
files = dir(datadir);

for i = 3:length(files) %skip over . and ..
        clearvars phfos
        
        % get animal & session # from file name
        % file name format: [animal]_[datatype]times[session].txt
        underscores = strfind(files(i).name,'_');
        t = strfind(lower(files(i).name),'times');
        animal = files(i).name(1:underscores(end)-1);
        type = files(i).name(underscores(end)+1:t-1);
        dot = strfind(files(i).name,'.');
        epochfilter{1} = {'task',['$sess==' files(i).name(dot-2:dot-1)]};
        
        % get channel from chinfos
        f = createfilter('animal', animal, 'epochs', epochfilter, 'data', datafilter);
        d = f(1).epochs{1}(1,1);
        ch = f(1).data.chinfo{1}{1}(1);
        
        % get session start time
        animaldir = [animalbasedir '/' animal];
        taskname = sprintf('%s/%stask.mat',animaldir,animal);
        load(taskname);
        start = task{d}{1}.start/1000;
        
        events.samprate = 5000;
        
        % read start & end times from file
        try
            times = dlmread([files(i).folder '/' files(i).name]);
            times = times(find(times(:,1)>0),:);
            events.starttime = times(:,1)+start;
            events.endtime = times(:,2)+start;
        catch % if file is empty
            events.starttime = [];
            events.endtime = [];
        end
        
        % save 2 types of HFOs
        phfos{d}{1}{ch} = events;
        save(sprintf('%s/%sphfos%02d.mat', animaldir, animal, d), 'phfos');
end