function extract_spikes(datadir, animaldir, prefix, day, epoch, varargin)
cd(datadir)

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
slashind = strfind(filemaskpath,'\');
fullname = filemaskpath(slashind(end-1)+1:end);

chinfoname = sprintf('%s/%schinfo.mat',animaldir,prefix);
load(chinfoname);

ccref = 0;
%find first corpus callosum site & load eeg
for c = 1:length(chinfo{day}{epoch})
    if isequal(chinfo{day}{epoch}{c}.area,'cc')
        ccref = remap(c);
        break
    end
end
if ccref==0
    error('Set corpus callosum reference site in chinfo')
end
configfile = sprintf('ch%02d_ref_config_datalogger_20kHz.trodesconf',ccref);

extractSpikeBinaryFiles(datadir, fullname, configfile) %extraction code adds '.rec' back

%%add this back  for newer Trodes versions
fullname = filemaskpath(slashind(end)+1:end);
Spikepath = [filemaskpath '.spikes/' fullname '.spikes'];
starttime = task{day}{epoch}.start/1000;
endtime = task{day}{epoch}.end/1000;

for i = 1:32 %read & store data for 1 channel
    Spikepath_chan = sprintf('%s_nt%d.dat',Spikepath, i);
    Spike_chan{i} = readTrodesExtractedDataFile(Spikepath_chan);
    if ~isempty(Spike_chan{i}.fields(1).data)%i~=ccref
        s = 1;
        e = length(Spike_chan{i}.fields(1).data);
        while s<=e && double(Spike_chan{i}.fields(1).data(s))/Spike_chan{i}.clockrate<starttime
            s=s+1;
        end
        while e>0 && double(Spike_chan{i}.fields(1).data(e))/Spike_chan{i}.clockrate>endtime
            e=e-1;
        end
        if s>e %no spikes in time window
            Spike_chan{i}.fields(1).data = [];
            Spike_chan{i}.fields(2).data = [];
        else
            Spike_chan{i}.fields(1).data = Spike_chan{i}.fields(1).data(s:e);
            Spike_chan{i}.fields(2).data = Spike_chan{i}.fields(2).data((s:e),:);
        end
    end
end

for i = 1:32 %read & store data for 1 channel
    clearvars mua
    
    mua{day}{epoch}{i}.timerange = [starttime endtime];
    mua{day}{epoch}{i}.starttime = starttime;
    mua{day}{epoch}{i}.endtime = endtime;
    mua{day}{epoch}{i}.samprate = Spike_chan{i}.clockrate;
    mua{day}{epoch}{i}.nTrode = i;
    mua{day}{epoch}{i}.ref = Spike_chan{i}.referencentrode;
    mua{day}{epoch}{i}.spiketimes = double(Spike_chan{remap(i)}.fields(1).data);
    mua{day}{epoch}{i}.waveform = double(Spike_chan{remap(i)}.fields(2).data);
    mua{day}{epoch}{i}.descript = 'MUA extracted with exportSpikes, 600-6000Hz filtered, CC reference';
    
    savename = sprintf('%s/%smua%02d-%d-%02d.mat',animaldir,prefix,day,epoch,i);
    save(savename,'mua','-v7.3');
end

end