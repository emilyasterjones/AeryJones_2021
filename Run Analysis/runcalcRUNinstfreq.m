%% Construct filter
%%%  NOTE must have chronux NOT on path otherwise wrong findpeaks function
%%%  is used!
clearvars -except f

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
epochfilter{1} = {'task','(contains($env,''LT'') && isequal($treatment, ''CNO'') && ~contains($descript,''FAIL''))'};

datafilter = []; 
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
velbins = [1 50];

probecrit = '(~isequal($area,''dead''))';

theta = [5 11];
thetacoeff = designeegfilt(1000,theta(1),theta(2));
SG = [20 50];
SGcoeff = designeegfilt(1000,SG(1),SG(2));
FG = [50 110];
FGcoeff = designeegfilt(1000,FG(1),FG(2));

% Specify function and run
for v = 1:length(velbins)-1
    timefilter{1} = {'pos', sprintf('(($vel > %d) & ($vel < %d))',velbins(v),velbins(v+1))};
    g = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
    g = setfilterfunction(g, 'calcRUNinstfreq', {'eeg','chinfo'},'probe',probecrit,'thetafilt',thetacoeff,'SGfilt',SGcoeff,'FGfilt',FGcoeff);
    g = runfilter(g);
    f{v} = g;
end

%% calculate results

%CA1 pyr, sr, slm
%CA3 pyr/sr
%DG mol, gc/hil
% location = '(isequal($area,''ca1'') && contains($layer,''pyr''))';
location = '(isequal($area,''ca1'') && contains($layer,''sr''))';
% location = '(isequal($area,''ca1'') && contains($layer,''slm''))';
% location = '(isequal($area,''ca3'') && (contains($layer,''pyr'') | contains($layer,''sr'')))';
% location = '( isequal($area,''dg'') && contains($layer,''mol''))';
% location = '( isequal($area,''dg'') && (contains($layer,''gc'')||contains($layer,''hil'')))';

animnames = [];
animgeno = [];
animvirus = [];

%STDout printing defaults
format compact
format shortG

% calculate power & normalized power
for v = 1:length(velbins)-1
    g = f{v};
    nanimals = length(g);
    for a = 1:length(g)
        results = g(a).output.calcRUNinstfreq.results;
        %load chinfo file to lookup site locations
        animinfo = animaldef(g(a).animal{1});
        infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
        load(infofile)
        
        animnames = [animnames; animinfo(3)];
        animgeno = [animgeno; animinfo(4)];
        animvirus = [animvirus; animinfo(5)];
        
        %identify sites in location of interest
        pwrtemp = evaluatefilter(chinfo,location);
        sites{a} = unique(pwrtemp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
        
        thetafreq{a}{v} = [];
        SGfreq{a}{v} = [];
        FGfreq{a}{v} = [];
        for c = 1:length(siteinds)  %iterate thru selected sites
            %combine data over epochs
            thetafreq{a}{v}{c} = [];
            SGfreq{a}{v}{c} = [];
            FGfreq{a}{v}{c} = [];
            for e = 1:length(results{1})
                [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                if ~isempty(results{1}{e}.thetaavginstfreq) && c<=length(siteinds)
                    thetafreq{a}{v}{c} = [thetafreq{a}{v}{c}; results{1}{e}.thetainstfreq{siteinds(c)}];
                    SGfreq{a}{v}{c} = [SGfreq{a}{v}{c}; results{1}{e}.SGinstfreq{siteinds(c)}];
                    FGfreq{a}{v}{c} = [FGfreq{a}{v}{c}; results{1}{e}.FGinstfreq{siteinds(c)}];
                end
            end
        end
        thetameanfreqs(a,v) = mean(cell2mat(thetafreq{a}{v}));
        SGmeanfreqs(a,v) = mean(cell2mat(SGfreq{a}{v}));
        FGmeanfreqs(a,v) = mean(cell2mat(FGfreq{a}{v}));
    end
    
    fprintf('Velocity > %d and < %d cm/s\n',velbins(v),velbins(v+1))
    thetameanfreqs(:,v)
    SGmeanfreqs(:,v)
    FGmeanfreqs(:,v)
end

for a = 1:length(g)
    thetafreq{a} = thetafreq{a}{1};
    trim = min(cellfun(@length,thetafreq{a}));
    for c = 1:length(thetafreq{a})
        thetafreq{a}{c} = thetafreq{a}{c}(1:trim);
    end
    allthetafreqs{a} = cell2mat(thetafreq{a});
    avgthetafreqs{a} = nanmean(allthetafreqs{a},2);
    allthetafreqs{a} = allthetafreqs{a}(:);
    
    SGfreq{a} = SGfreq{a}{1};
    trim = min(cellfun(@length,SGfreq{a}));
    for c = 1:length(SGfreq{a})
        SGfreq{a}{c} = SGfreq{a}{c}(1:trim);
    end
    allSGfreqs{a} = cell2mat(SGfreq{a});
    avgSGfreqs{a} = nanmean(allSGfreqs{a},2);
    allSGfreqs{a} = allSGfreqs{a}(:);
    
    FGfreq{a} = FGfreq{a}{1};
    trim = min(cellfun(@length,FGfreq{a}));
    for c = 1:length(FGfreq{a})
        FGfreq{a}{c} = FGfreq{a}{c}(1:trim);
    end
    allFGfreqs{a} = cell2mat(FGfreq{a});
    avgFGfreqs{a} = nanmean(allFGfreqs{a},2);
    allFGfreqs{a} = allFGfreqs{a}(:);
end