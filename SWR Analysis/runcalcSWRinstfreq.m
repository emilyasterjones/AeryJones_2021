%% Construct filter

%%%  NOTE must have chronux NOT on path otherwise wrong findpeaks function
%%%  is used!

animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',...
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''FAIL'') && $dur==120)'};

datafilter = [];  %specify area of interest. will automatically use ca1pyr1 for rip detection
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);

probe = '(isequal($area,''ca1'')&& contains($layer,''pyr 1''))';

%Specify function and run
f = setfilterfunction(f, 'calcSWRinstfreq', {'eeg','ripples','chinfo'},'appendindex',1,'probe',probe,'minstd',5);
f = runfilter(f);

%% plot results!

% clearvars -except f probe animals gnames

% ------------------------------
% Figure and Font Sizes
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'DefaultAxesFontName','Arial')
tfont = 20; % title font
xfont = 20;
yfont = 20;
% load('C:\Users\HuangShared\Documents\MATLAB\EJ\Gain of Toxic Function\colormap_ekopaper.mat')
% load('C:\Users\HuangShared\Documents\MATLAB\EJ\Gain of Toxic Function\semcolormap_ekopaper.mat')

freqbins = 125:1:200;
cyclebins = 3:1:60;

animnames = [];
animgeno = [];
cols = [];
histfreqs = [];
histcycles = [];

for a = 1:length(f)
    results = f(a).output.calcSWRinstfreq.results;
    if ~isempty(results{1}{1}.instfreqs)
        for p = 1:length(results{1}{1}.instfreqs)
            instfreqs{a}{p} = [];
            cycles{a}{p} = [];
            lengths{a}{p} = [];
            for e = 1:length(results{1})  %iterate through epochs, combine counts
                instfreqs{a}{p} = [instfreqs{a}{p}; results{1}{e}.avginstfreqs{p}];
                cycles{a}{p} = [cycles{a}{p}; results{1}{e}.cyclesperrip{p}];
                lengths{a}{p} = [lengths{a}{p}; results{1}{e}.riplength];
            end
            freqsbychan{a}(p,:) = histc(instfreqs{a}{p},freqbins)/length(instfreqs{a}{p});
            cyclesbychan{a}(p,:) = histc(cycles{a}{p},cyclebins)/length(cycles{a}{p});
            
        end
        meanfreqs(a) = mean(instfreqs{a}{1});
        meancycles(a) = mean(cycles{a}{1});
        meanlengths(a) = mean(lengths{a}{1});
        
        animnames = [animnames; f(a).animal(3) ];
        cols = [cols; f(a).animal{5} ];
        animgeno = [animgeno; f(a).animal(4)];
    end
end

for a = 1:length(f)
    cycles{a} = cycles{a}{1};
    instfreqs{a} = instfreqs{a}{1};
end
