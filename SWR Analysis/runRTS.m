myCluster = parcluster('local');
myCluster.NumWorkers = 6;
saveProfile(myCluster);

animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',...
    'Scotty', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''10mg/kg'') && ~contains($descript,''*FAIL*'') && $dur==120)'};

datafilter = [];  % load all channels...use embedded filter in run function to get specific ones
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

% set varargin for RTspecgram3 (channel filtering criteria)
trigcrit = '(isequal($area,''ca1'') && contains($layer,''pyr 1''))';
probecrit = '(~isequal($area,''dead''))';

% visualize for each animal
for a = animals
    f = createfilter('animal', a{1}, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
    f = setfilterfunction(f, 'RTspecgram3', {'eeg','ripples','chinfo'},'minstd',3,'trig',trigcrit,'probe',probecrit,'slidingwin',[.1 .01]);
    f = runfilter(f);
    plotRTS
end

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
f = setfilterfunction(f, 'RTspecgram3', {'eeg','ripples','chinfo'},'minstd',5,'trig',trigcrit,'probe',probecrit,'slidingwin',[.1 .1]);
f = runfilter(f);