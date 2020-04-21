clearvars -except f figopt
figopt = 2;

% location = '(isequal($area,''dg'') && contains($layer,''*mua*''))';
location = '(isequal($area,''ca1'') && contains($layer,''sr''))';
% location = '( isequal($area,''ca3'')&& (contains($layer,''pyr'') || contains($layer,''sr'')))';
% location = '( isequal($area,''dg'') && contains($layer,''val'') && (contains($layer,''gc'') || contains($layer,''hil'')))';

% Figure and Font Sizes
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'DefaultAxesFontName','Arial')

%figopts:
%2 = average power at 1 timepoint during event (eg SG power immediately
%after ripple onset)
%3 = curve of power, baseline, & SD over all freqs
%4 = baseline & SD
%5 = same as 2 except with matched baselines & SD between 2 epochs
%6 = same as 9 except with matched baselines & SD between 2 epochs
%7 = bootstrapped timecourse
%8 = correlation between different subregions
%9 = distribution of power z-score over events (eg inc/dec fraction)
%10 = [baseline power] to [during event power] in uV
%11 = plot power in every site for each mouse
%12 = distribution of SG baseline power in time bins

%% calculate  power in various freq bands at one timepoint for specific areas
if figopt==2
    set(0,'DefaultAxesFontName','Arial')
    animnames = [];
    animgeno = [];
    cols = [];
    clearvars pwr
    
    for a = 1:length(f)
        results = f(a).output.RTspecgram3.results;
        %load chinfo file to lookup site locations
        animinfo = animaldef(f(a).animal{1});
        infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
        load(infofile)
        %identify sites in location of interest
        temp = evaluatefilter(chinfo,location);
        sites{a} = unique(temp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
        
        probe{a} = [];
        pwr{a} = [];
        for c = 1:length(siteinds)  %iterate thru selected sites
            %combine data over epochs
            probe{a}{c} = [];
            pwr{a}{c} = [];
%             lastind = 0;
            times = results{1}{1}.times;
            freqs = results{1}{1}.freqs;
            for e = 1:length(results{1})
                [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                if length(results{1}{e}.freqs)>1 && c<=length(siteinds)
                    probe{a}{c} = cat(3,probe{a}{c},results{1}{e}.Sprobe{1}{siteinds(c)});
                    pwr{a}{c} = [pwr{a}{c}; squeeze(nanmean(results{1}{e}.Sprobe{1}{siteinds(c)}(5,4:7,:)))];
                    eprobe{a}{c}{e} = results{1}{e}.Sprobe{1}{siteinds(c)};
                end
                %end
            end
        end
        
        pwr{a} = cell2mat(pwr{a});
        avgpwr{a} = nanmean(pwr{a},2);
        allpwr{a} = pwr{a}(:);
    end
end

%% calculate z-scored power across all bands at one timepoint for specific areas
if figopt==3
    set(0,'DefaultAxesFontName','Arial')
    animnames = [];
    animgeno = [];
    cols = [];
    freqs = f(1).output.RTspecgram3.results{1}{1}.freqs;
    for a = 1:length(f)
        results = f(a).output.RTspecgram3.results;
        %load chinfo file to lookup site locations
        infofile = sprintf('%s%schinfo.mat',f(a).animal{2},f(a).animal{3});
        load(infofile)
        %identify sites in location of interest
        temp = evaluatefilter(chinfo,location);
        sites{a} = unique(temp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
        
        for c = 1:length(siteinds)  %iterate thru selected sites
            %combine data over epochs
            probe{a}{c} = [];
            baseline{a}{c} = [];
            stds{a}{c} = [];
            for e = 1:length(results{1})
                if results{1}{e}.riptimes(1)>0
                    probe{a}{c} = cat(3,probe{a}{c},results{1}{e}.Sprobe{1}{siteinds(c)});
                    baseline{a}{c} = cat(3,baseline{a}{c},results{1}{e}.meanPprobe{1}{siteinds(c)});
                    stds{a}{c} = cat(3,stds{a}{c},results{1}{e}.stdPprobe{1}{siteinds(c)});
                end
            end
        end
        
        animnames = [animnames; f(a).animal(3) ];
        cols = [cols; f(a).animal{5} ];
        animgeno = [animgeno; f(a).animal(4)];
    end
    
    %get rid of animals with no sites in requested area
    validanims = ~cellfun(@isempty, sites);
    animnames = animnames(validanims);
    cols = cols(validanims);
    animgeno = animgeno(validanims);
    probe = probe(validanims);
    baseline = baseline(validanims);
    stds = stds(validanims);
    
    temp = [];
    for g = 1:length(gnames)
        match = regexp(animgeno,gnames(g));
        temp=[temp g*(~cellfun(@isempty,match))];
    end
    groupnum = sum(temp,2);
    
    timeind = 5;  %time ms after thresh crossing
    avgpwr = zeros(length(animnames), length(freqs));
    avgbase = zeros(length(animnames), length(freqs));
    avgstd = zeros(length(animnames), length(freqs));
    
    for a = 1:length(animnames)
        for freq = 1:length(freqs)
            temppwr = [];
            tempbase = [];
            tempstd = [];
            for c = 1:length(probe{a}) %for each channel
                temppwr = [temppwr mean(squeeze(mean(probe{a}{c}(timeind,freq,:),2)))];
                tempbase = [tempbase mean(baseline{a}{c}(1,freq,:))];
                tempstd = [tempstd mean(stds{a}{c}(1,freq,:))];
            end
            avgpwr(a,freq) = nanmean(temppwr);
            avgbase(a,freq) = nanmean(tempbase);
            avgstd(a,freq) = nanmean(tempstd);
        end
    end
    
    %Z-score
    figure
    bins = freqs;
    for g = 1:length(gnames)
        aggdist{g} = [];
        for i = find(groupnum==g)'
            aggdist{g} = [aggdist{g}; avgpwr(i,:)];
        end
        hold on
        sem = std(aggdist{g})/sqrt(size(aggdist{g},1));
        h = fill([bins fliplr(bins)], [(mean(aggdist{g},1)-sem) fliplr(mean(aggdist{g},1)+sem)],sem_colors(gindex(g),:),'FaceAlpha',.5);
        set(h,'EdgeColor','none')
        lines(g) = plot(bins,nanmean(aggdist{g},1),'Color',colors(gindex(g),:),'Linewidth',3,'DisplayName',gnames{g});
    end
    xlabel('Frequency')
    ylabel('Power (z-score)')
    legend(lines)
    legend('boxoff')
    
    %Baseline
    figure
    for g = 1:length(gnames)
        aggdist{g} = [];
        for i = find(groupnum==g)'
            aggdist{g} = [aggdist{g}; avgbase(i,:)];
        end
        hold on
        sem = std(aggdist{g})/sqrt(size(aggdist{g},1));
        h = fill([bins fliplr(bins)], [(mean(aggdist{g},1)-sem) fliplr(mean(aggdist{g},1)+sem)],sem_colors(gindex(g),:),'FaceAlpha',.5);
        set(h,'EdgeColor','none')
        lines(g) = plot(bins,nanmean(aggdist{g},1),'Color',colors(gindex(g),:),'Linewidth',3,'DisplayName',gnames{g});
    end
    xlabel('Frequency')
    ylabel('Baseline (uV)')
    legend(lines)
    legend('boxoff')
    
    %SD
    figure
    for g = 1:length(gnames)
        aggdist{g} = [];
        for i = find(groupnum==g)'
            aggdist{g} = [aggdist{g}; avgstd(i,:)];
        end
        hold on
        sem = std(aggdist{g})/sqrt(size(aggdist{g},1));
        h = fill([bins fliplr(bins)], [(mean(aggdist{g},1)-sem) fliplr(mean(aggdist{g},1)+sem)],sem_colors(gindex(g),:),'FaceAlpha',.5);
        set(h,'EdgeColor','none')
        lines(g) = plot(bins,nanmean(aggdist{g},1),'Color',colors(gindex(g),:),'Linewidth',3,'DisplayName',gnames{g});
    end
    xlabel('Frequency')
    ylabel('SD (uV)')
    legend(lines)
    legend('boxoff')
end

%% plot zscore baseline and stdev
%only works with one condition

if figopt ==4
    animnames = [];
    cols = [];
    animgeno = [];
    basechan = [];
    sdchan = [];
    baseline = [];
    sd = [];
    
    for a = 1:length(f)
        results = f(a).output.RTspecgram3.results;
        %load chinfo file to lookup site locations
        animinfo = animaldef(f(a).animal{1});%f(a).animal;
        infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
        load(infofile)
        %identify sites in location of interest
        temp = evaluatefilter(chinfo,location);
        sites{a} = unique(temp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
        clear temp
        
        basechan = [];
        sdchan = [];
        if(~isempty(siteinds))
            for g = 1:length(results)  %iterate thru conditions
                for c = 1:length(siteinds)  %iterate thru selected sites
                    baseepochs = [];
                    sdepochs = [];
                    freqs = results{g}{1}.freqs;
                    for e = 1:length(results{g}) %combine data over epochs
                        [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                        baseepochs = cat(3, baseepochs, results{g}{e}.meanPprobe{1}{siteinds(c)});
                        sdepochs = cat(3, sdepochs, results{g}{e}.stdPprobe{1}{siteinds(c)});
                    end
                    basechan = [basechan; mean(baseepochs,3)];  %{a}(c,:) combine data over epochs
                    sdchan = [sdchan; mean(sdepochs,3)];
                end
                baseline = [baseline; mean(basechan,1)];%(a,:) combine data over channels
                sd = [sd; mean(sdchan,1)];
            end
        else
            baseline = [baseline; nan(1,44)];
            sd = [sd; nan(1,44)];
        end
        animnames = [animnames; f(a).animal(3) ];
        cols = [cols; f(a).animal{5} ];
        animgeno = [animgeno; f(a).animal(4)];
        
    end
    bands = [4 7];
    bandnames = {'slow gamma'};
    
    for b = 1:size(bands,1)
        bandstart = bands(b,1);
        bandend = bands(b,2);
        
        for a = 1:length(animnames)
            meanbaseline(a) = mean2(baseline(a,bandstart:bandend)*double(12500/65536));
            meansd(a) = mean2(sd(a,bandstart:bandend)*double(12500/65536));
        end
        disp(bandnames{b})
        meanbaseline'
        meansd'
        
    end
end
%%
if figopt==5
    for a = 1:length(f)
        results = f(a).output.RTspecgram3.results;
        vehresults = f1(a).output.RTspecgram3.results;
        
        animinfo = animaldef(f(a).animal{1});
        infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
        load(infofile)
        
        temp = evaluatefilter(chinfo,location);
        sites{a} = unique(temp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),sites{a});
        probe{a} = [];
        pwr{a} = repmat(-5000,5000,length(siteinds));
        
        for c = 1:length(siteinds)  %iterate thru selected sites
            %combine data over epochs
            probe{a}{c} = [];
            for e = 1:min(length(results{1}),length(vehresults{1}))
                [sites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),sites{a});
                if ~isempty(results{1}{e}.riptimes) && c<=length(siteinds)
                    cnobase = results{1}{e}.meanPprobe{1}{siteinds(c)};
                    vehbase = vehresults{1}{e}.meanPprobe{1}{siteinds(c)};
                    cnostd = results{1}{e}.stdPprobe{1}{siteinds(c)};
                    vehstd = vehresults{1}{e}.stdPprobe{1}{siteinds(c)};
                    z = results{1}{e}.Sprobe{1}{siteinds(c)}; %z-scored by cno
                    newz = ((z.*cnostd+cnobase)-vehbase)./vehstd; %re z-scored by veh
                    probe{a}{c} = cat(3,probe{a}{c},newz);
                    pwr{a}(lastind+1:length(results{1}{e}.riptimes)+lastind,c) = squeeze(nanmean(newz(5,4:7,:)));
                    
                end
            end
        end
        pwr{a} = nanmean(pwr{a},2);
        pwr{a} = pwr{a}(pwr{a}>-4000);
    end
end