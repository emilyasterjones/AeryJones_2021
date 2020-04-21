% Set: animal names, site area & layer, source folders for raw and
% preprocessed data and destination folders for .mat and .csv files, target
% site and data extraction option

% From raw data, extract 5kHz sampled EEG saved as .mat and .csv
%       extracts full traces, sleep traces, & awake traces

clearvars

animals = {'Adobo_LT', 'Anise_LT','Baharat_LT', 'Basil_LT', 'Caraway_LT',...
    'Cardamom_LT', 'Chicory_LT', 'Chives_LT', 'Cinnamon_LT', 'Coriander_LT',...
    'Cumin_LT', 'Dill_LT', 'Fenugreek_LT', 'GaramMasala_LT', 'Ginger_LT',...
    'Jerk_LT', 'Mace_LT', 'Mustard_LT','Nutmeg_LT', 'OldBay_LT', 'Oregano_LT',...
    'Paprika_LT', 'Parsley_LT', 'Pepper_LT', 'Provence_LT', 'Pumpkin_LT',...
    'Rosemary_LT', 'Saffron_LT', 'Sage_LT', 'Salt_LT', 'Sumac_LT', 'Tarragon_LT', 'Thyme_LT', 'Vanilla_LT'};

basedir = 'E:\CRCNS\Cohort 2';

rawdatadir = [basedir 'Trodes Raw Data\'];
animaldir = [basedir 'Preprocessed Data\'];
eegdir = 'E:\CRCNS\Cohort 2\';
csvdir = 'E:\LCRCNS\Cohort 2';
filterstring = '(isequal($area,''ca1'') && contains($layer,''pyr 1''))';
opt = 1;

epochfilter = [];
epochfilter{1} = {'task','(contains($treatment,''Veh'') && ~isequal($descript,''FAIL''))'};
datafilter = [];
timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

for a = 1:length(animals)
    clearvars f
    %get reference channel
    datafilter{1} = {'chinfo','(isequal($area,''cc''))'};
    f = createfilter('animal', [animals{a},''], 'epochs', epochfilter, 'data', datafilter);
    ccref = f(1).data.chinfo{1}{1}(1);
    clearvars f
    
    %get # sessions, site #, & sleep times
    datafilter{1} = {'chinfo',filterstring};
    try
        f = createfilter('animal', [animals{a},''], 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
    catch
        disp(['No CA1 channel found for ',animals{a}]);
    end
    
    if exist('f', 'var') %if there was a channel of the specified type on this animal
        taskname = sprintf('%s/%s/%stask.mat',basedir,animals{a},animals{a});
        load(taskname);
        
        for i = 1:size(f(1).epochs{1},1)
            d = f(1).epochs{1}(i,1); %session # is d, index is i
            ch = f(1).data.chinfo{1}{i}(1); %channel
            
            %extract data at 5kHz sampling, 1kHz filter
            eegref = extract_eeg_1kHzfilt(rawdatadir,...
                sprintf('%s\\%s\\',animaldir,animals{a}),sprintf('%s\\%s\\',eegdir,animals{a}),animals{a},d,ccref);
            eeg = extract_eeg_1kHzfilt(rawdatadir,...
                sprintf('%s\\%s\\',animaldir,animals{a}),sprintf('%s\\%s\\',eegdir,animals{a}),animals{a},d,ch);
            eeg{d}{1}{ch}.data = eeg{d}{1}{ch}.data - eegref{d}{1}{ccref}.data;
            eegfile = sprintf('%s/Referenced/%s/%seeg%02d-1-%02d.mat',eegdir,animals{a},animals{a},d,ch);
            save(eegfile, 'eeg')

            %full trace
            csv_file = sprintf('%s/Full ref cell layer/%s%02d.csv', csvdir, animals{a}, d);
            dlmwrite(csv_file,eeg{d}{1}{ch}.data, '\t');
            
            %calculate sleep times
            starttemp = task{d}{1}.start/1000;
            endtemp = task{d}{1}.end/1000;
            Fs = (length(eeg{d}{1}{ch}.data)-1)/(endtemp-starttemp);
            fulltimes = starttemp:1/Fs:endtemp;
            sleeptimes = ~isExcluded(fulltimes, f(1).excludetime{1}{i});
            
            %sleep trace
            csv_file = sprintf('%s/Sleep ref cell layer/%s%02d_sleeponly.csv', csvdir, animals{a}, d);
            dlmwrite(csv_file,eeg{d}{1}{ch}.data(sleeptimes), '\t');
            %awake trace
            csv_file = sprintf('%s/Awake ref cell layer/%s%02d_awakeonly.csv', csvdir, animals{a}, d);
            dlmwrite(csv_file,eeg{d}{1}{ch}.data(~sleeptimes), '\t');
        end
    end
end