% clearvars -except animals gnames f figopt

animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',...
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''FAIL'') && $dur==120)'}; %or 10mg/kg

datafilter = [];
datafilter{1} = {'chinfo','((isequal($area,''ca1'') && contains($layer,''pyr 1'')) || (isequal($area,''ca1'') && contains($layer,''sr'')))'};

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime >30)'};

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
%interripwin should at least be as large as win
f = setfilterfunction(f, 'calcSWamp', {'eeg','ripples','chinfo'}, 'win', [0 0], 'interripwin', 0, 'minstd', 5);
f = runfilter(f);

%% Extract summary

animnames = [];
animgeno = [];
avgSWamps = [];
for a = 1:length(f)
    results = f(a).output.calcSWamp.results;
    
    %identify sites in location of interest
    location = '(isequal($area,''ca1'') && contains($layer,''sr''))';
    infofile = sprintf('%s%schinfo.mat',f(a).animal{2},f(a).animal{3});
    load(infofile)
    temp = evaluatefilter(chinfo,location);
    sites{a} = unique(temp(:,3));
    [sites{a},siteinds] = intersect(results{1}{1}.index.chinfo,sites{a});
    
    SWamps = [];
    allSWamps{a} = [];
    for e = 1:length(results{1})  %iterate through epochs, combine counts
        chanSWamps = NaN(length(siteinds),length(results{1}{e}.peak{1}));
        [sites{a},siteinds] = intersect(results{1}{e}.index.chinfo,sites{a});
        for c = 1:length(siteinds)  %selected channels
            chanSWamps(c,:) = results{1}{e}.peak{siteinds(c)};
        end
        allSWamps{a} = [allSWamps{a}; (mean(chanSWamps)*double(12500/65536))'];
        SWamps = [SWamps; mean(squeeze(nanmean(chanSWamps)))];
    end
    avgSWamps(a) = mean(SWamps)*double(12500/65536);
    animnames = [animnames; f(a).animal(3) ];
    animgeno = [animgeno; f(a).animal(4)];
end

avgSWamps'