function ej_rereference(directoryname,fileprefix,day,epoch)
%REREFERENCE(directoryname,fileprefix,dayss)
%
%Applies corpus callosum as a reference to each channel
%
%directoryname - example '/data99/user/animaldatafolder', a folder
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers

chinfoname = sprintf('%s/%schinfo.mat',directoryname,fileprefix);
load(chinfoname);

ccref = 0;
%find first corpus callosum site & load eeg
for c = 1:length(chinfo{day}{epoch})
    if isequal(chinfo{day}{epoch}{c}.area,'cc')
        ccref = c;
        break
    end
end
if ccref==0
    error('Set corpus callosum reference site in chinfo')
end

refname = sprintf('%s/%seeg%02d-%d-%02d.mat', directoryname, fileprefix, day, epoch, ccref);
load(refname);
refdata = eeg{day}{epoch}{ccref}.data;

% create the list of files for this day that we should filter
flist = dir(sprintf('%s/*eeg%02d-%d-*.mat', directoryname, day, epoch));

% go through each file in flist and reference it
for fnum = 1:length(flist)
    % get the tetrode number and epoch
    dash = find(flist(fnum).name == '-');
    endidx = length(dash);
    tet = str2num(flist(fnum).name((dash(endidx)+1):(dash(endidx)+3)));
    
    %load the eeg file, reference, & save
    load(sprintf('%s/%s', directoryname, flist(fnum).name));
    eeg{day}{epoch}{tet}.data = eeg{day}{epoch}{tet}.data - refdata;
    savename = sprintf('%s/%seeg%02d-%d-%02d.mat',directoryname,fileprefix,day,epoch,tet);
    save(savename,'-v6','eeg');
    clear eeg
end
