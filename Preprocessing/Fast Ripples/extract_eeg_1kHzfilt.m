function extract_eeg(datadir, animaldir, prefix, day, epoch, varargin)

trackepoch = 0;
multtracks = 0;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'trackepoch'
                trackepoch = varargin{option+1};
            case 'multtracks' %if running 2 tracks on 1 ECU
                multtracks = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

taskname = sprintf('%s/%stask.mat',animaldir,prefix);
load(taskname);
nepochs = length(task{day});

%original UH32
remap = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32];

%build options for finding the file
if nepochs>1
    multepochs = 1;
else
    multepochs = 0;
end
if epoch==trackepoch
    suffix = '_merge.rec';
else
    suffix = '.rec';
end

multepochs = 1;

%get the file location
filemaskpath = findTrackFilePath(datadir, animaldir, prefix, day, epoch, suffix, multtracks, multepochs);
filemaskpath = filemaskpath(1:end-4);
slashind = strfind(filemaskpath,'/');
fullname = filemaskpath(slashind(end)+1:end);

extractLFPBinaryFiles(filemaskpath) %extraction code adds '.rec' back
LFPpath = [filemaskpath '.LFP/' fullname '.LFP'];
timestampspath = [filemaskpath '.LFP/' fullname '.timestamps.dat'];

timestamps_temp.fields.data = double(timestamps_temp.fields.data)/timestamps_temp.clockrate;

starttemp = task{day}{epoch}.start/1000;
endtemp = task{day}{epoch}.end/1000;
startidx = lookup(starttemp,timestamps_temp.fields.data);
endidx = lookup(endtemp,timestamps_temp.fields.data);
fprintf('%s\tpos:\t%f\t%f\tLFP:\t%f\t%f\n',fullname,starttemp,endtemp,...
    timestamps_temp.fields.data(1),timestamps_temp.fields.data(end))

timestamps_temp.clockrate = 30000; %20kHz for linear track
downSampleFactor = timestamps_temp.clockrate/5000;
timestamps_temp.fields.data = timestamps_temp.fields.data(startidx:endidx);
starttime = timestamps_temp.fields.data(1);
endtime = timestamps_temp.fields.data(end);
bandpassLow = .1;
bandpassHigh = 1000;
[c,d] = butter(2, [bandpassLow bandpassHigh]/(timestamps_temp.clockrate/2));

parpool('local', 6); %run the following loop on 6 cores
parfor i = 1:32 %read, format, filter, & downsample data for 1 channel
    LFPpath_chan = sprintf('%s_nt%dch1.dat',LFPpath, i);
    LFP_chan = readTrodesExtractedDataFile(LFPpath_chan);
    
    channelData = LFP_chan.fields.data*double(12500/65536); %convert to uV
    channelData = int16(filtfilt(c,d,double(channelData)));   %filter first
    channelData = downsample(channelData(startidx:endidx),downSampleFactor);   %then downsample
    
    eeg_temp{i}.timerange = [starttime endtime];
    eeg_temp{i}.starttime = starttime;
    eeg_temp{i}.endtime = endtime;
    eeg_temp{i}.samprate = 1000;
    eeg_temp{i}.nTrode = i;
    eeg_temp{i}.nTrodeChannel = 1;
    eeg_temp{i}.data = channelData;
    eeg_temp{i}.descript = 'No description';
end
delete(gcp); %shut down the parallel pool

%save to correct structure, separate from parfor loop to prevent
%classification error
for i = 1:32
    clearvars eeg;
    eeg{day}{epoch}{i} = eeg_temp{remap(i)};
    savename = sprintf('%s/%seeg%02d-%d-%02d.mat',animaldir,prefix,day,epoch,i);
    save(savename,'-v6','eeg');
end

end