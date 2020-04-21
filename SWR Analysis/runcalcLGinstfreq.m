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
%
%Specify function and run
f = setfilterfunction(f, 'calcLGinstfreq3', {'eeg','ripples','chinfo'},'appendindex',1,'probe',probe,'minstd',5);
f = runfilter(f);

figopt=2;

%% plot results!

% ------------------------------
% Figure and Font Sizes
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'DefaultAxesFontName','Arial')
tfont = 20; % title font
xfont = 20;
yfont = 20;

if figopt==1   % Make histogram of inst freqs
    bins = 15:.5:70;
    cbins = 0:1:20;
    %figure
    %aggregate data over animals
    animnames = [];
    animgeno = [];
    cols = [];
    avgcounts = [];
    avgcycles = [];
    avgfreqorder = [];
    for a = 1:length(f)
        results = f(a).output.calcLGinstfreq3.results;
        if ~isempty(results{1}{1}.instfreqs)
            for p = 1:length(results{1}{1}.instfreqs)
                
                cycles{a}{p} = [];
                instfreqs{a}{p} = [];
                freqorder{a}{p} = nan(20,1)';
                for e = 1:length(results{1})  %iterate through epochs, combine counts
                    instfreqs{a}{p} = [instfreqs{a}{p}; results{1}{e}.avginstfreqs{p}];
                    cycles{a}{p} = [cycles{a}{p}; results{1}{e}.cyclesperrip{p}];
                    freqorder{a}{p} = [freqorder{a}{p}; results{1}{e}.instfreqs{p}];
                end
                countsbychan{a}(p,:) = histc(instfreqs{a}{p},bins)/length(instfreqs{a}{p});
                cyclesbychan{a}(p,:) = histc(cycles{a}{p},cbins)/length(cycles{a}{p});
                freqorderbychan{a}(p,:) = nanmean(freqorder{a}{p});
            end
            avgcounts = [avgcounts; nanmean(countsbychan{a})];
            avgcycles = [avgcycles; nanmean(cyclesbychan{a})];
            avgfreqorder = [avgfreqorder; nanmean(freqorderbychan{a})];
            
            animnames = [animnames; f(a).animal(3) ];
            cols = [cols; f(a).animal{5} ];
            animgeno = [animgeno; f(a).animal(4)];
            
        end
    end
    
    temp = [];
    for g = 1:length(gnames)
        match = regexp(animgeno,gnames(g));
        temp=[temp g*(~cellfun(@isempty,match))];
    end
    groupnum = sum(temp,2);
    
    [gmeans,gstds,gsems] = grpstats(avgcounts,groupnum,{'mean','std','sem'});
    nonz = find(gmeans(1,:)>0);
    % plot by genotype (avg of indiv mice)
    figure
    for g = 1:length(gnames)
        hold on
        h = fill([bins(nonz) fliplr(bins(nonz))], [gmeans(g,nonz)-gsems(g,nonz) fliplr(gmeans(g,nonz)+gsems(g,nonz))],sem_colors(g,:),'FaceAlpha',.5);
        set(h,'EdgeColor','none')
        lines(g) = plot(bins(nonz),gmeans(g,nonz),'Color',colors(g,:),'Linewidth',3,'DisplayName',gnames{g});
    end
    
    ylabel('Proportion of Power')
    xlabel('Instantaneous Frequency')
    xlim([20 60])
    legend(lines)
    legend('boxoff')
    
    % plot avg instfreq by mouse and by geno
elseif figopt ==2
    %     clearvars -except f figopt probe animals gnames
    clearvars instfreqs
    animnames = [];
    animgeno = [];
    meanfreq = [];
    for a = 1:length(f)
        results = f(a).output.calcLGinstfreq3.results;
        if ~isempty(results{1}{1}.instfreqs)
            instfreqs{a} = [];
            for p = 1:length(results{1}{1}.instfreqs)
                lastind = 0;
                instfreqs{a}{p} = [];
                for e = 1:length(results{1})  %iterate through epochs, combine counts
                    if p<=length(results{1}{e}.instfreqs)
                        instfreqs{a}{p} = [instfreqs{a}{p}; results{1}{e}.avginstfreqs{p}];
                    end
                end
            end
            
            animnames = [animnames; f(a).animal(3) ];
            animgeno = [animgeno; f(a).animal(4)];
            
        end
    end
    meanfreq'
    
    for a = 1:length(f)
        allinstfreqs{a} = cell2mat(instfreqs{a});
        avginstfreqs{a} = nanmean(allinstfreqs{a},2);
        avginstfreqs{a} = avginstfreqs{a}(~isnan(avginstfreqs{a}));
        allinstfreqs{a} = allinstfreqs{a}(:);
        allinstfreqs{a} = allinstfreqs{a}(~isnan(allinstfreqs{a}));
    end
end
