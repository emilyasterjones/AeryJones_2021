animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',...
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(isequal($treatment, ''veh'') && contains($env,''PRE'') && contains($descript,''thresh'') && ~contains($descript,''FAIL''))'}; %or CNO
datafilter = [];  %only works with single channel specified
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};
timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'}; %sleep

trigcrit = '(isequal($area,''ca1'') && contains($layer,''pyr 1''))';
probecrit = '(contains($layer,''mua''))';

% Veh
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''*FAIL*'') && $dur==120)'};
f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
f = setfilterfunction(f, 'RTMUA', {'mua','ripples','chinfo'},'minstd',3,'trig',trigcrit,'probe',probecrit);
f = runfilter(f);


%% PLOT
location = '(isequal($area,''ca3''))';

set(0,'defaultaxesfontweight','normal')
set(0,'defaultaxeslinewidth',2)
set(0,'defaultaxesfontsize',16)
set(0,'DefaultAxesFontName','Arial')
tfont = 20; % title font
xfont = 20;
yfont = 20;

animnames = [];
animgeno = [];
clearvars baserate SWRrate normrate normpsth
for a = 1:length(f)
    results = f(a).output.RTMUA.results;
    fullbaserate{a} = [];
    baserate{a} = [];
    sdfullbaserate{a} = [];
    sdbaserate{a} = [];
    psth{a} = [];
    timebase{a} = [];
    wave{a} = [];
    peakrate{a} = [];
    SWRrate{a} = [];
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
                    fullbaserate{a} = [fullbaserate{a}; results{1}{e}.fullbaserate{t}{p}];
                    baserate{a} = [baserate{a}; results{1}{e}.baserate{t}{siteinds(p)}];
                    sdfullbaserate{a} = [sdfullbaserate{a}; results{1}{e}.sdfullbaserate{t}{p}];
                    sdbaserate{a} = [sdbaserate{a}; results{1}{e}.sdbaserate{t}{p}];
                    psth{a} = [psth{a}; results{1}{e}.psth{t}{p}];
                    normpsth{a} = [normpsth{a}; (results{1}{e}.psth{t}{p}-results{1}{e}.baserate{t}{p})./results{1}{e}.sdbaserate{t}{p}];
                    timebase{a} = [timebase{a}; results{1}{e}.psth{t}{p}(:,5)];
                    wave{a} = [wave{a}; results{1}{e}.wave{t}{p}];
                    peakrate{a} = [peakrate{a}; results{1}{e}.peakrate{t}{p}];%(:,1)
                    SWRrate{a} = [SWRrate{a}; results{1}{e}.SWRrate{t}{siteinds(p)}];
                    normrate{a} = [normrate{a}; results{1}{e}.normrate{t}{siteinds(p)}];
                end
            end
        end
    end
    normrate{a} = normrate{a}(normrate{a}>0);
    animnames = [animnames; f(a).animal(3)];
    animgeno = [animgeno; f(a).animal(4)];
end