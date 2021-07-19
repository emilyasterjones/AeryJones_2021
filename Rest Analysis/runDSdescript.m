clearvars -except animals epochfilter

animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''*FAIL*'') && $dur==120)'};

datafilter = [];  %only works with single channel specified
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};
trigcrit = '(isequal($area,''dg'') && contains($layer,''mua''))';

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);

%Specify function and run
f = setfilterfunction(f, 'DSdescript', {'dspikes','eeg','chinfo'},'trig',trigcrit);
f = runfilter(f);

%% Parse data
for a = 1:length(f)
    results = f(a).output.DSdescript.results;
    lengths{a} = [];
    amps{a} = [];
    timecounts{a} = [];
    for e = 1:length(results{1})  %iterate through epochs, combine counts
        lengths{a} = [lengths{a}; results{1}{e}.dslength];
        amps{a} = [amps{a}; results{1}{e}.dspeakamps];
        timecounts{a} = [timecounts{a}; results{1}{e}.timecounts'];
    end
end