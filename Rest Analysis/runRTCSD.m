% clearvars -except f figopt

animals = {'Aphid', 'Fruitfly', 'Grasshopper', 'Honeybee', 'Rolypoly',...
    'Cicada', 'Cricket', 'Ladybug', 'Mantis'};

epochfilter = [];
epochfilter{1} = {'task','(contains($treatment,''CNO'') && ~isequal($descript,''FAIL''))'};

datafilter = [];  % load all channels...use embedded filter in run function to get specific ones
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

trigcrit = '(isequal($area,''ca1'') && contains($layer,''pyr 1''))';

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
f = setfilterfunction(f, 'RTCSD', {'eeg','ripples','chinfo'},'trig',trigcrit,'minstd',5,'win',[0 0.1],'freqband',[20 50]); %,'freqband',[20 50]
f = runfilter(f);
% save('\\hub.gladstone.internal\HuangLab-LFP\Emily\DREADDs Revisions\Analysis\2021 CNO SG RTCSD.mat','f', '-v7.3')

%% plot & calculate
for a = 1:length(f)
    results = f(a).output.RTCSD.results{1};
    nchan = size(results{1}.ripcsd,2);
    win = size(results{1}.ripcsd,3);
    timewin = results{1}.win;
    ripCSD{a} = [];
    for e = 1:length(results)  %iterate through epochs
        ripCSD{a} = [ripCSD{a}; results{e}.ripcsd];
    end
    
    % plot average
    figure
    CSD = squeeze(nanmean(ripCSD{a}));
    M = max(max(abs(CSD))); % abosolute maximum CSD, for the colormap scale
    clims = [-M M];
    im = imagesc(CSD(3:end-2,:),clims);
    colormap(jet); % blue = sink; red = source
    cb = colorbar('EastOutside');
    set(cb,'YTick',[]);
    set(gca,'ycolor','none')
    set(gca,'TickDir','out');
    xticks(1:100:win)
    xticklabels(timewin(1):0.1:timewin(2))
    ylabel('Site');
    xlabel('Time (s)');
    title(sprintf('%s (sink=blue, source=red)',f(a).animal{3}));
    
    %add traces
    hold on
    for i = 3:nchan-2
        plot(-CSD(i,:)/M+(i-2),'k')
        hold on
    end
    
    %average over locations
    location = '(isequal($area,''ca1'') && contains($layer,''pyr''))';
    
    animinfo = animaldef(f(a).animal{1});%f(a).animal;
    infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
    load(infofile)
    temp = evaluatefilter(chinfo,location);
    sites{a} = unique(temp(:,3));
    
    peakripsource{a} = [];
    peakripsink{a} = [];
    zripsource{a} = [];
    zripsink{a} = [];
    
    for e = 1:length(results)
        for p = 1:length(sites{a})            
            peakripsource{a} = [peakripsource{a}; results{e}.peakripsource(:,sites{a}(p))];
            peakripsink{a} = [peakripsink{a}; results{e}.peakripsink(:,sites{a}(p))];
            zripsource{a} = [zripsource{a}; results{e}.zripsource(:,sites{a}(p))];
            zripsink{a} = [zripsink{a}; results{e}.zripsink(:,sites{a}(p))];
        end
    end
end