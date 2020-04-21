function ej_rippledayprocess(directoryname,fileprefix,day,epoch)
%RIPPLEDAYPROCESS(directoryname,fileprefix,days, options)
%
%Applies a ripple filter to all epochs for each day and saves the data in
%in the EEG subdirectory of the directoryname folder.
%
%directoryname - example '/data99/user/animaldatafolder', a folder
%                containing processed matlab data for the animal
%
%fileprefix -    animal specific prefix for each datafile (e.g. 'fre')
%
%days -          a vector of experiment day numbers

load('C:\Users\HuangShared\Documents\MATLAB\Filters\ejmouseripplefilter.mat'); % 125-200

minint = -32768;

% create the list of files for this day that we should filter
flist = dir(sprintf('%s/*eeg%02d-%d-*.mat', directoryname, day, epoch));

% go through each file in flist and filter it
for fnum = 1:length(flist)
    % get the tetrode number and epoch
    % MODIFIED 9/22 to accomodate paths with dashes
    dash = find(flist(fnum).name == '-');
    endidx = length(dash);
    tet = str2num(flist(fnum).name((dash(endidx)+1):(dash(endidx)+3)));
    
    %load the eeg file
    load(sprintf('%s/%s', directoryname, flist(fnum).name));
    a = find(eeg{day}{epoch}{tet}.data < -30000);
    [lo,hi]= findcontiguous(a);  %find contiguous NaNs
    for i = 1:length(lo)
        if lo(i) > 1 & hi(i) < length(eeg{day}{epoch}{tet}.data)
            fill = linspace(double(eeg{day}{epoch}{tet}.data(lo(i)-1)), ...
                double(eeg{day}{epoch}{tet}.data(hi(i)+1)), hi(i)-lo(i)+1);
            eeg{day}{epoch}{tet}.data(lo(i):hi(i)) = int16(fill);
        end
    end
    % filter it and save the result as int16
    ripple{day}{epoch}{tet} = filtereeg2(eeg{day}{epoch}{tet}, ...
        ripplefilter, 'int16', 1);
    % replace the filtered invalid entries with the minimum int16 value of
    % -32768
    for i = 1:length(lo)
        if lo(i) > 1 && hi(i) < length(ripple{day}{epoch}{tet}.data)
            ripple{day}{epoch}{tet}.data(lo(i):hi(i)) = minint;
        end
    end
    
    % convert to DF2.0 format (keep fields field for more info)
    ripple{day}{epoch}{tet}.filteredamp = ripple{day}{epoch}{tet}.data(:,1);
    ripple{day}{epoch}{tet}.instphase = ripple{day}{epoch}{tet}.data(:,2);
    ripple{day}{epoch}{tet}.envmag = ripple{day}{epoch}{tet}.data(:,3);
    ripple{day}{epoch}{tet} = rmfield(ripple{day}{epoch}{tet}, 'data');
    
    % save the resulting file
    ripplefile = sprintf('%s/%sripple%02d-%d-%02d.mat', ...
        directoryname, fileprefix, day, epoch, tet);
    save(ripplefile, 'ripple');
    clear ripple
end
