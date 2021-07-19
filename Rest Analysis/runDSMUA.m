% clearvars -except f location
clearvars -except epochfilter animals

animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''*FAIL*'') && $dur==120)'};

datafilter = [];  %only works with single channel specified
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'}; %sleep

trigcrit = '(isequal($area,''dg'') && contains($layer,''mua''))';
probecrit = '(contains($layer,''mua''))';

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
f = setfilterfunction(f, 'DSMUA', {'mua','dspikes','chinfo'},'trig',trigcrit,'probe',probecrit);
f = runfilter(f);


%% PLOT
location = '(isequal($area,''ca1''))';
clearvars baserate DSrate normrate normpsth
for a = 1:length(f)
    results = f(a).output.DSMUA.results;
    baserate{a} = [];
    DSrate{a} = [];
    normrate{a} = [];
    normpsth{a} = [];
    
    animinfo = animaldef(f(a).animal{1});%f(a).animal;
    infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
    load(infofile)
    temp = evaluatefilter(chinfo,location);
    sites{a} = unique(temp(:,3));
    [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
    
    for e = 1:length(results{1})  %iterate through epochs, combine counts
        if ~isempty(results{1}{e}.probeindex) && ~isempty(results{1}{e}.trigindex)
            for t = 1:length(results{1}{e}.baserate)
                for p = 1:length(siteinds)%1:length(results{1}{e}.baserate{t})
                    [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                    baserate{a} = [baserate{a}; results{1}{e}.baserate{t}{siteinds(p)}];
                    normpsth{a} = [normpsth{a}; (results{1}{e}.psth{t}{p}-results{1}{e}.baserate{t}{p})./results{1}{e}.sdbaserate{t}{p}];
                    DSrate{a} = [DSrate{a}; results{1}{e}.DSrate{t}{siteinds(p)}];
                    normrate{a} = [normrate{a}; results{1}{e}.normrate{t}{siteinds(p)}];                    
                end
            end
        end
    end
end