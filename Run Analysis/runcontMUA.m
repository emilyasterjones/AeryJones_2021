animals = {'Adobo_LT', 'Caraway_LT', 'Chives_LT', ...
    'Cinnamon_LT', 'Coriander_LT', 'Fenugreek_LT', 'GaramMasala_LT', 'Salt_LT', ...
    'Baharat_LT', 'Cardamom_LT', 'Jerk_LT', 'Mace_LT', 'Mustard_LT', ...
    'Tarragon_LT', 'Vanilla_LT',...
    'Basil_LT', 'Cumin_LT', 'Dill_LT', 'Nutmeg_LT', 'Paprika_LT', 'Parsley_LT', 'Sumac_LT',...
    'Pepper_LT', 'Sage_LT', 'Anise_LT', 'Thyme_LT', 'OldBay_LT', 'Rosemary_LT',...
    'Provence_LT', 'Saffron_LT'};
analysisdir = '\\hub.gladstone.internal\HuangLab-LFP\Emily\DREADDs WMaze\Analysis';
group = 'LT';

epochfilter = [];
epochfilter{1} = {'task','(contains($env,''LT'') && isequal($treatment, ''veh'') && ~contains($descript,''FAIL''))'};

datafilter = [];  % load all channels...use embedded filter in run function to get specific ones
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
timefilter{1} = {'pos','($vel>1)'};

% set varargin for RTspecgram3 (channel filtering criteria)
probecrit = '(contains($layer,''**mua**''))';

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
f = setfilterfunction(f, 'contMUA', {'mua','chinfo'},'probe',probecrit);
f = runfilter(f);


%% PLOT
location1 = '(isequal($area,''ca1''))';
location2 = '(isequal($area,''dg''))';

set(0,'defaultaxesfontweight','normal')
set(0,'defaultaxeslinewidth',2)
set(0,'defaultaxesfontsize',16)
set(0,'DefaultAxesFontName','Arial')
tfont = 20; % title font
xfont = 20;
yfont = 20;

R = NaN(1,length(f));
p = NaN(1,length(f));
for a = 1:length(f)
    results = f(a).output.contMUA.results;
    loc1{a} = [];
    loc2{a} = [];
    
    %load chinfo
    animinfo = animaldef(f(a).animal{1});%f(a).animal;
    infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
    load(infofile)
    
    %location 1
    %get site indices
    temp = evaluatefilter(chinfo,location1);
    sites{a} = unique(temp(:,3));
    [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
    
    %iterate through epochs, combine counts
    for e = 1:length(results{1})
        if ~isempty(results{1}{e}.probeindex)
            for p = 1:length(siteinds)
                [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                loc1{a} = [loc1{a}; results{1}{e}.rate{siteinds(p)}'];
            end
        end
    end
    
    %repeat for location 2
    %get site indices
    temp = evaluatefilter(chinfo,location2);
    sites{a} = unique(temp(:,3));
    [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
    
    %iterate through epochs, combine counts
    for e = 1:length(results{1})
        if ~isempty(results{1}{e}.probeindex)
            for p = 1:length(siteinds)
                [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                loc2{a} = [loc2{a}; results{1}{e}.rate{siteinds(p)}'];
            end
        end
    end
    
    if ~isempty(loc1{a}) && ~isempty(loc2{a})
        [R(a),p(a)] = corr(loc1{a},loc2{a},'rows','complete');
    end
end
R'