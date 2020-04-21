clearvars -except f location

animals = {'Adobo_LT', 'Caraway_LT', 'Chives_LT', ...
    'Cinnamon_LT', 'Coriander_LT', 'Fenugreek_LT', 'GaramMasala_LT', 'Salt_LT', ...
    'Baharat_LT', 'Cardamom_LT', 'Jerk_LT', 'Mace_LT', 'Mustard_LT', ...
    'Tarragon_LT', 'Vanilla_LT',...
    'Basil_LT', 'Cumin_LT', 'Dill_LT', 'Nutmeg_LT', 'Paprika_LT', 'Parsley_LT', 'Sumac_LT',...
    'Pepper_LT', 'Sage_LT', 'Anise_LT', 'Thyme_LT', 'OldBay_LT', 'Rosemary_LT',...
    'Provence_LT', 'Saffron_LT'};
analysisdir = '\\hub.gladstone.internal\HuangLab-LFP\Emily\DREADDs WMaze\Analysis';
group = 'PV CRC';

epochfilter = [];
epochfilter{1} = {'task','(contains($env,''LT'') && isequal($treatment, ''CNO'') && ~contains($descript,''FAIL''))'};

datafilter = [];  % load all channels...use embedded filter in run function to get specific ones
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
velbins = [1 50];

probecrit = '(~isequal($area,''dead''))';

calculate baseline & SD across all movement
timefilter{1} = {'pos', '($vel>1)'};
g = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
g = setfilterfunction(g, 'continuous_specgram', {'eeg','chinfo'},'probe',probecrit,'slidingwin',[1 1]);
g = runfilter(g);
f{1} = g;

for v = 1:length(velbins)-1
    timefilter{1} = {'pos', sprintf('(($vel > %d) & ($vel < %d))',velbins(v),velbins(v+1))};
    g = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
    g = setfilterfunction(g, 'continuous_specgram', {'eeg','chinfo'},'probe',probecrit,'slidingwin',[1 1]);
    g = runfilter(g);
    f{v+1} = g;
end

figopt=2;

%% PLOT
%figopts:
%1 = plot average power across all frequencies, raw and Z-scored
%2 = average power in frequency bins at different speeds, raw and Z-scored
%3 = relationship between speed and power in frequency bins, raw and Z-scored
%4 = baseline & SD
%5 = correlation between regions
%6 = distributions of power in frequency bins

%CA1 pyr, sr, slm
%CA3 pyr/sr
%DG mol, gc/hil
% location = '(isequal($area,''ca1'') && contains($layer,''pyr''))';
location = '(isequal($area,''ca1'') && contains($layer,''sr''))';
% location = '(isequal($area,''ca1'') && contains($layer,''slm''))';
% location = '(isequal($area,''ca3'') && (contains($layer,''pyr'') | contains($layer,''sr'')))';
% location = '( isequal($area,''dg'')&& contains($layer,''mol''))';
% location = '( isequal($area,''dg'') && (contains($layer,''gc'')||contains($layer,''hil'')))';

theta = [5 11];
slowgamma = [20 50];
fastgamma = [50 110];

% Figure and Font Sizes
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'DefaultAxesFontName','Arial')
format compact
format shortG

% calculate base & SD across all movement
g = f{1};
for a = 1:length(g)
    nanimals = length(g);
    results = g(a).output.continuous_specgram.results;
    %load chinfo file to lookup site locations
    animinfo = animaldef(g(a).animal{1});
    infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
    load(infofile)
    
    %identify sites in location of interest
    pwrtemp = evaluatefilter(chinfo,location);
    sites{a} = unique(pwrtemp(:,3));
    [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
    for c = 1:length(siteinds)
        base{a}{c} = [];
        sd{a}{c} = [];
        for e = 1:length(results{1})
            [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
            if length(results{1}{e}.freqs)>1 && c<=length(siteinds)
                base{a}{c} = cat(3,base{a}{c},results{1}{e}.meanPprobe{siteinds(c)});
                sd{a}{c} = cat(3,sd{a}{c},results{1}{e}.stdPprobe{siteinds(c)});
            end
        end
    end
end

% calculate power & normalized power
for v = 1:length(velbins)-1
    % % %     g = f{v+1};
    g = f{v};
    for a = 1:length(g)
        results = g(a).output.continuous_specgram.results;
        %load chinfo file to lookup site locations
        animinfo = animaldef(g(a).animal{1});
        infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
        load(infofile)
        
        %identify sites in location of interest
        pwrtemp = evaluatefilter(chinfo,location);
        sites{a} = unique(pwrtemp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
        
        pwr{a}{v} = [];
        distpwr{a}{v} = [];
        normpwr{a}{v} = [];
        for c = 1:length(siteinds)  %iterate thru selected sites
            %combine data over epochs
            pwr{a}{v}{c} = [];
            distpwr{a}{v}{c} = [];
            normpwr{a}{v}{c} = [];
            for e = 1:length(results{1})
                [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                if length(results{1}{e}.freqs)>1 && c<=length(siteinds)
                    chanpwr = results{1}{e}.Pprobe{siteinds(c)};%*double(12500/65536);
                    pwr{a}{v}{c} = cat(1,pwr{a}{v}{c},mean(chanpwr,1));
                    distpwr{a}{v}{c} = cat(1,distpwr{a}{v}{c},chanpwr);
                    normchanpwr = chanpwr;
                    for i = 1:length(results{1}{e}.freqs)
                        normchanpwr(:,i) = (chanpwr(:,i)-base{a}{c}(1,i,e))./sd{a}{c}(1,i,e);
                    end
                    normpwr{a}{v}{c} = cat(3,normpwr{a}{v}{c},mean(normchanpwr,1));
                end
            end
        end
    end
end

% find freqs that match bands
g = f{1};
freqs = g(a).output.continuous_specgram.results{1}{1}.freqs;
bandnames = {'theta', 'slow gamma', 'fast gamma'};
bands = [lookup(theta(1),freqs) lookup(theta(2),freqs); ...
    lookup(slowgamma(1),freqs) lookup(slowgamma(2),freqs); ...
    lookup(fastgamma(1),freqs) lookup(fastgamma(2),freqs)];

if figopt==1 %plot average power across all frequencies, raw and Z-scored
    %Raw Power
    for a = 1:nanimals
        figure
        hold on
        for v = 1:length(velbins)-1
            pwrtemp = [];
            for c = 1:length(pwr{a}{v}) %for each channel
                pwrtemp = [pwrtemp; squeeze(mean(pwr{a}{v}{c}(1,:,:),3))];
            end
            semilogy(freqs,mean(pwrtemp,1))
            xlim([0 125])
            labels{v} = sprintf('%d-%dcm/s',velbins(v),velbins(v+1));
        end
        %         legend(labels{1},labels{2},labels{3},labels{4})
    end
    
    %Normalized Power
    for a = 1:nanimals
        figure
        hold on
        for v = 1:length(velbins)-1
            pwrtemp = [];
            for c = 1:length(pwr{a}{v}) %for each channel
                pwrtemp = [pwrtemp; squeeze(mean(normpwr{a}{v}{c}(1,:,:),3))];
            end
            semilogy(freqs,mean(pwrtemp,1))
            xlim([0 125])
            labels{v} = sprintf('%d-%dcm/s',velbins(v),velbins(v+1));
        end
        legend(labels{1},labels{2},labels{3},labels{4})
    end
    
elseif figopt==2 %average power in frequency bins at different speeds, raw and Z-scored
    for v = 1:length(velbins)-1
        avgpwr = NaN(length(nanimals),length(bands));
        %         avgnormpwr = NaN(length(nanimals),length(bands));
        fprintf('Velocity > %d and < %d cm/s\n',velbins(v),velbins(v+1))
        for b = 1:size(bands,1)
            bandstart = bands(b,1);
            bandend = bands(b,2);
            
            for a = 1:nanimals
                pwrtemp = [];
                normpwrtemp = [];
                for c = 1:length(pwr{a}{v}) %for each channel
                    allpwrtemp{a}(b,c,:) = mean(distpwr{a}{v}{c}(:,bandstart:bandend),2);
                    pwrtemp = [pwrtemp mean(squeeze(mean(pwr{a}{v}{c}(1,bandstart:bandend,:),2)))];
                    normpwrtemp = [normpwrtemp mean(squeeze(mean(normpwr{a}{v}{c}(1,bandstart:bandend,:),2)))];
                end
                
                if ~isempty(pwr{a}{v})
                    avgpwr(a,b) = mean(pwrtemp);
                    avgnormpwr(a,b) = mean(normpwrtemp);
                end
            end
        end
        avgpwr
        avgnormpwr
        disp('Slow vs Fast Gamma Ratio')
        avgpwr(:,2)/avgpwr(:,3)'
        avgnormpwr(:,2)/avgnormpwr(:,3)'
        
        for a = 1:nanimals
            allthetapwr{a} = squeeze(allpwrtemp{a}(1,:,:));
            allthetapwr{a} = allthetapwr{a}(:);
            avgthetapwr{a} = mean(squeeze(allpwrtemp{a}(1,:,:)))';
            allSGpwr{a} = squeeze(allpwrtemp{a}(2,:,:));
            avgSGpwr{a} = mean(squeeze(allpwrtemp{a}(2,:,:)))';
            allFGpwr{a} = squeeze(allpwrtemp{a}(3,:,:));
            avgFGpwr{a} = mean(squeeze(allpwrtemp{a}(3,:,:)))';
            avgSGFGratio{a} = mean(allSGpwr{a}./allFGpwr{a})';
            allSGpwr{a} = allSGpwr{a}(:);
            allFGpwr{a} = allFGpwr{a}(:);
            allSGFGratio{a} = allSGpwr{a}./allFGpwr{a};
        end
        
        
    end
    
elseif figopt==3 %relationship between speed and power in frequency bins, raw and Z-scored
    for b = 1:size(bands,1)
        bandstart = bands(b,1);
        bandend = bands(b,2);
        avgpwr = NaN(length(nanimals),length(velbins)-1);
        avgnormpwr = NaN(length(nanimals),length(velbins)-1);
        
        for v = 1:length(velbins)-1
            for a = 1:nanimals
                pwrtemp = [];
                normpwrtemp = [];
                for c = 1:length(pwr{a}{v}) %for each channel
                    pwrtemp = [pwrtemp mean(squeeze(mean(pwr{a}{v}{c}(1,bandstart:bandend,:),2)))];
                    normpwrtemp = [normpwrtemp mean(squeeze(mean(normpwr{a}{v}{c}(1,bandstart:bandend,:),2)))];
                end
                
                if ~isempty(pwr{a}{v})
                    avgpwr(a,v) = mean(pwrtemp);
                    avgnormpwr(a,v) = mean(normpwrtemp);
                end
            end
        end
        for a = 1:nanimals
            [R(a),p(a)] = corr(velbins(1:end-1)',avgpwr(a,:)');
            [nR(v),np(v)] = corr(velbins(1:end-1),avgnormpwr(a,:));
        end
        bandnames{b}
        R'
        p'
        nR'
        np'
    end
    
elseif figopt==4 %baseline & SD
    avgbase = NaN(length(nanimals),length(bands)-1);
    avgsd = NaN(length(nanimals),length(bands)-1);
    for b = 1:size(bands,1)
        bandstart = bands(b,1);
        bandend = bands(b,2);
        
        for a = 1:nanimals
            basetemp = [];
            sdtemp = [];
            for c = 1:length(base{a}) %for each channel
                basetemp = [basetemp mean(squeeze(mean(base{a}{c}(1,bandstart:bandend,:),2)))];
                sdtemp = [sdtemp mean(squeeze(mean(sd{a}{c}(1,bandstart:bandend,:),2)))];
            end
            
            if ~isempty(base{a})
                avgbase(a,b) = mean(basetemp);
                avgsd(a,b) = mean(sdtemp);
            end
        end
    end
    avgbase
    avgsd
    
elseif figopt==5 %correlation across regions
    %     location1 = '(isequal($area,''ca1'') && contains($layer,''pyr''))';
    %     location1 = '(isequal($area,''ca1'') && contains($layer,''sr''))';
    location1 = '(isequal($area,''ca1'') && contains($layer,''slm''))';
    %     location2 = '(isequal($area,''ca3'') && (contains($layer,''pyr'') | contains($layer,''sr'')))';
    % location = '( isequal($area,''dg'')&& contains($layer,''mol''))';
    location2 = '( isequal($area,''dg'') && (contains($layer,''gc'')||contains($layer,''hil'')))';
    
    g = f{1};
    for b = 1:size(bands,1)
        bandstart = bands(b,1);
        bandend = bands(b,2);
        for a = 1:nanimals
            nanimals = length(g);
            results = g(a).output.continuous_specgram.results;
            freqs = results{1}{1}.freqs;
            
            %load chinfo file to lookup site locations
            animinfo = animaldef(g(a).animal{1});
            infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
            load(infofile)
            
            %location 1
            %identify sites
            pwrtemp = evaluatefilter(chinfo,location1);
            sites{a} = unique(pwrtemp(:,3));
            [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
            
            %collect power across all time windows for each site
            loc1{a} = [];
            for c = 1:length(siteinds)
                loc1{a}{c} = [];
                for e = 1:length(results{1})
                    [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                    if length(results{1}{e}.freqs)>1 && c<=length(siteinds)
                        loc1{a}{c} = [loc1{a}{c}; mean(results{1}{e}.Pprobe{siteinds(c)}(:,bandstart:bandend),2)];
                    end
                end
            end
            loc1avg{a} = mean(cell2mat(loc1{a}),2);
            
            %repeat for location 2
            %identify sites
            pwrtemp = evaluatefilter(chinfo,location2);
            sites{a} = unique(pwrtemp(:,3));
            [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
            
            %collect power across all time windows for each site
            loc2{a} = [];
            for c = 1:length(siteinds)
                loc2{a}{c} = [];
                for e = 1:length(results{1})
                    [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                    if length(results{1}{e}.freqs)>1 && c<=length(siteinds)
                        loc2{a}{c} = [loc2{a}{c}; mean(results{1}{e}.Pprobe{siteinds(c)}(:,bandstart:bandend),2)];
                    end
                end
            end
            loc2avg{a} = mean(cell2mat(loc2{a}),2);
            
            [R(a),p(a)] = corr(loc1avg{a},loc2avg{a});
        end
        bandnames{b}
        R'
        p'
    end
    
elseif figopt==6 %distribution of power over time windows
    clearvars avgpwr histpwr histratio
    pwrbins = 0:100:10000; %0:20:1000; %
    pwrcenters = pwrbins(1:end-1)+50; %+10
    ratiobins = 0:1:30;
    ratiocenters = ratiobins(1:end-1)+0.5;
    
    for v = 1:length(velbins)-1
        fprintf('Velocity > %d and < %d cm/s\n',velbins(v),velbins(v+1))
        for b = 1:size(bands,1)
            bandstart = bands(b,1);
            bandend = bands(b,2);
            
            for a = 1:nanimals
                avgpwr{a}{b} = [];
                for c = 1:length(distpwr{a}{v}) %avg over freqs in bin
                    avgpwr{a}{b} = [avgpwr{a}{b} mean(distpwr{a}{v}{c}(:,bandstart:bandend),2)];
                end
                avgpwr{a}{b} = mean(avgpwr{a}{b},2); %avg over channels
                histpwr(a,:) = histcounts(avgpwr{a}{b},pwrbins)/length(avgpwr{a}{b});
            end
            
            figure %avg over animals
            plot(pwrcenters,nanmean(histpwr),'Color',[0.543 0 0],'Linewidth',1); %,[0.543 0 0]
            hold on
            sem = nanstd(histpwr,0,1)/sqrt(9);
            h = fill([pwrcenters fliplr(pwrcenters)], [(nanmean(histpwr)-sem) fliplr(nanmean(histpwr)+sem)],[0.543 0 0],'FaceAlpha',.5); %[.543 0 0]
            set(h,'EdgeColor','none')
            title(location)
            box off
            hold on
            %             plot(pwrcenters,histpwr')
        end
        
        %repeat for slow/fast gamma ratio
        for a = 1:nanimals
            histratio(a,:) = histcounts(avgpwr{a}{2}./avgpwr{a}{3},ratiobins)/length(avgpwr{a}{2});
            sdratio(a) = std(avgpwr{a}{2}./avgpwr{a}{3});
        end
        
        location
        sdratio'
        
        figure %avg over animals
        histratio = histratio(9:15,:);
        plot(ratiocenters,nanmean(histratio),'Color',[0.543 0 0],'Linewidth',1); %,[0.543 0 0]
        hold on
        sem = nanstd(histratio,0,1)/sqrt(8);
        h = fill([ratiocenters fliplr(ratiocenters)], [(nanmean(histratio)-sem) fliplr(nanmean(histratio)+sem)],[0.543 0 0],'FaceAlpha',.5); %[.543 0 0]
        set(h,'EdgeColor','none')
        title(location)
        box off
    end
end