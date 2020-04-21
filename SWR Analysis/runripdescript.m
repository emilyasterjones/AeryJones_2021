%% Construct filter
animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',...
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''FAIL'') && $dur==120)'}; %or 10mg/kg

datafilter = [];  %only works with single channel specified
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};
trigcrit = '(isequal($area,''ca1'') && contains($layer,''pyr 1''))';
probecrit = trigcrit;

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);

%Specify function and run
f = setfilterfunction(f, 'ripdescript', {'ripples','chinfo'},'minstd',5,'trig',trigcrit,'probe',probecrit);
f = runfilter(f);


%% plot results!
%ONLY WORKS WITH ONE CHANNEL PER MOUSE

% Figure and Font Sizes
set(0,'defaultaxesfontweight','normal');
set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'DefaultAxesFontName','Arial')
tfont = 20; % title font
xfont = 20;
yfont = 20;
% Parse data
animnames = [];
animgeno = [];
for a = 1:length(f)
    results = f(a).output.ripdescript.results;
    rates{a} = [];
    excluded{a} = [];
    ripsizes{a} = [];
    riplengths{a} = [];
    rippeaks{a} = [];
    ripenergy{a} = [];
    iri{a} = [];
    for e = 1:length(results{1})  %iterate through epochs, combine counts
        rates{a} = [rates{a}; length(results{1}{e}.riplength)/results{1}{e}.validdur];
        excluded{a} = [excluded{a}; results{1}{e}.excluded(2)];
        ripsizes{a} = [ripsizes{a}; results{1}{e}.ripstds];
        riplengths{a} = [riplengths{a}; results{1}{e}.riplength];
        rippeaks{a} = [rippeaks{a}; results{1}{e}.rippeakamps];
        ripenergy{a} = [ripenergy{a}; mean(results{1}{e}.ripenergy)];
        if size(results{1}{e}.interriptimes,1)<size(results{1}{e}.interriptimes,2)
            results{1}{e}.interriptimes = results{1}{e}.interriptimes';
        end
        iri{a} = [iri{a}; results{1}{e}.interriptimes];
    end
    animnames = [animnames; f(a).animal(3) ];
    animgeno = [animgeno; f(a).animal(4)];
end

for a = 1:length(f)
    ischain{a} = iri{a}<.2;
end

for a = 1:length(f)
    avgrate(a) = mean(rates{a});
    avgsize(a) = mean(ripsizes{a});
    avglength(a) = mean(riplengths{a});
    avgexcluded(a) = mean(excluded{a});
    avgpeak(a) = mean(rippeaks{a});
    periri(a) = sum(iri{a}<.2)/length(iri{a})*100;
    perlong(a) = sum(riplengths{a}>100)/length(riplengths{a})*100;
end

avgexcluded'
avgrate'
avgsize'
avglength'
avgpeak'
periri'
perlong'