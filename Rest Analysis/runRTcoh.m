% clearvars -except f figopt

animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''*FAIL*'') && $dur==120)'};

datafilter = [];  % load all channels...use embedded filter in run function to get specific ones
datafilter{1} = {'chinfo','(~isequal($area,''cc''))'};

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);

% set varargin for RTcoh (channel filtering criteria)
trigcrit = '(isequal($area,''ca1'') && contains($layer,''*pyr 1*''))';
cohcrit = '( isequal($area,''ca1'')&& contains($layer,''**sr**''))';
probecrit = '( isequal($area,''ca3'')&& (contains($layer,''*pyr*'') || contains($layer,''*sr*'')))';
f = setfilterfunction(f, 'RTcoh', {'eeg','ripples','chinfo'},'trig',trigcrit,'coh',cohcrit,'probe',probecrit,'minstd',3, 'cwin',[.1 .1]);
f = runfilter(f);

%% plot results!
figopt = 4;
% clearvars -except f figopt gnames animals location

if figopt==1  %plot individual channels
    for a = 1:length(f)
        results = f(a).output.RTcoh.results;
        for g = 1:length(results)  %iterate thru conditions
            for t = 1:size(results{1}{1}.cohindex,1)  %iterate thru trigchans
                %combine data over epochs
                
                times = results{g}{1}.times;
                freqs = results{g}{1}.freqs;
                
                for p = 1:size(results{1}{1}.probeindex,1) %iterate through probechans
                    %iterate thru probe chans, combine eps, plot
                    
                    coh{t}{p} = [];
                    for e = 1:length(results{g}) %combine events from all eps
                        if length(results{g}{e}.freqs)>1
                            times = results{g}{e}.times;
                            freqs = results{g}{e}.freqs;
                            coh{t}{p} = cat(3,coh{t}{p},results{g}{e}.coherence{t}{p});
                        end
                    end
                    
                    meancoh{t}{p} = mean(coh{t}{p},3);
                    numevents = size(coh{t}{p},3);
                    figure
                    colormap hot
                    imagesc(times, freqs,meancoh{t}{p}');
                    colorbar;
                    set(gca, 'Clim', [.5 .9]);
                    set(gca,'YDir','normal');
                    details = sprintf('%s coh btwn %d and %d, %d 5SD rips on ch%d', ...
                        f(a).animal{3},results{1}{1}.cohindex(t,3),results{1}{1}.probeindex(p,3),numevents,results{1}{1}.cohindex(t,3));
                    title(details, 'FontSize',20,'Fontweight','normal');
                    ylabel('Freq','FontSize',20,'Fontweight','normal');
                    xlabel('Time(s)','FontSize',20,'Fontweight','normal');
                    set(gca,'XLim',[min(times) max(times)]);
                    set(gca,'YLim',[min(freqs) max(freqs)]);
                    % Plot Line at 0 ms - Start of ripple
                    hold on
                    ypts = freqs;
                    xpts = 0*ones(size(ypts));
                    plot(xpts , ypts, 'w--','Linewidth',2);
                end
            end
        end
    end
    
    
elseif figopt == 2  % plot all channels in correct configuration
    for a = 1:length(f)
        figure
        hold on
        colormap hot
        
        infofile = sprintf('%s%schinfo.mat',f(a).animal{2},f(a).animal{3});
        load(infofile)
        results = f(a).output.RTcoh.results;
        for g = 1:length(results)  %iterate thru conditions
            for t = 1:size(results{1}{1}.trigindex,1)  %iterate thru trigchans, should only be 1
                %combine data over epochs
                times = results{1}{1}.times;
                freqs = results{1}{1}.freqs;
                
                for p = 1:size(results{1}{1}.probeindex,1) %iterate through probechans
                    coh{t}{p} = [];
                    for e = 1:length(results{g}) %epochs
                        %combine events from all eps
                        if results{g}{e}.ripsizes(1)>0
                            coh{t}{p} = cat(3,coh{t}{p},results{g}{e}.coherence{t}{p});
                            %                             ecoh{t}{p}{e} = mean(squeeze(mean(results{g}{e}.coherence{t}{p}(5,4:7,:))));
                        end
                    end
                    numevents = size(coh{t}{p},3);
                    meancoh{t}{p} = mean(coh{t}{p},3);
                    
                    figure
                    plot(cell2mat(ecoh{t}{p}))
                    title(sprintf('%s ch%d %s %s', f(a).animal{1}, results{1}{1}.probeindex(p,3),...
                        chinfo{1}{1}{results{1}{1}.probeindex(p,3)}.area, chinfo{1}{1}{results{1}{1}.probeindex(p,3)}.layer));
                    fprintf('%s ch%d %s %s\t%f\t%f\n', f(a).animal{1}, results{1}{1}.probeindex(p,3),...
                        chinfo{1}{1}{results{1}{1}.probeindex(p,3)}.area,...
                        chinfo{1}{1}{results{1}{1}.probeindex(p,3)}.layer,...
                        mean(cell2mat(ecoh{t}{p}(vehicle))), mean(cell2mat(ecoh{t}{p}(highsess))));
                end
                
                shank = size(results{1}{1}.probeindex,1);
                
                if(~isempty(meancoh{t}{1}))
                    %four shank
                    set(gcf,'Position',[1 1 1595 964]);
                    disp('four sh')
                    for s = 1:shank
                        chanpos = results{1}{1}.probeindex(s,3);
                        
                        subplot(4,8,chanpos)
                        
                        imagesc(times, freqs(1:32),meancoh{t}{s}(:,1:32)');
                        set(gca, 'Clim', [.4 .7]);%
                        set(gca,'YDir','normal');
                        details = sprintf('%s %s', chinfo{1}{1}{chanpos}.area, chinfo{1}{1}{chanpos}.layer) ;
                        title(details, 'FontSize',10,'Fontweight','normal');
                        %ylabel('Freq','FontSize',10,'Fontweight','normal');
                        %xlabel('Time(s)','FontSize',10,'Fontweight','normal');
                        set(gca,'XLim',[min(times) max(times)]);
                        set(gca,'YLim',[min(freqs) max(freqs(1:32))]);
                        % Plot Line at 0 ms - Start of ripple
                        hold on
                        ypts = freqs(1:32);
                        xpts = 0*ones(size(ypts));
                        plot(xpts , ypts, 'k--','Linewidth',2);
                        colorbar;
                        
                    end
                    suptitle = sprintf('%s coh %drips on ch%d, 5SD in longhc immob>30',f(a).animal{3},numevents,results{g}{t}.trigindex(t,3));
                    [ax,h]=suplabel(suptitle ,'t');
                    tightfig;
                    
                    %end
                    
                    saveas(gcf, sprintf('%s CA1SR1 RTCoh.pdf',f(a).animal{3}))
                end
                
                
            end
        end
        
        
    end
    
    
    %plot quantification of low gamma coherence over timecourse
    
elseif figopt==3
    
    location = '(isequal($area,''ca1'')&& contains($layer,''pyr''))';%
    
    for a = 1:length(f)
        figure;
        hold on;
        results = f(a).output.RTcoh.results;
        %load chinfo struct to label site locations
        infofile = sprintf('%s%schinfo.mat',f(a).animal{2},f(a).animal{3});
        load(infofile)
        %identify sites in location of interest
        temp = evaluatefilter(chinfo,location);
        psites{a} = unique(temp(:,3));
        [psites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),psites{a});
        
        for g = 1:length(results)  %iterate thru conditions
            for t = 1:size(results{1}{1}.trigindex,1)  %iterate thru trigchans, should only be 1
                %combine data over epochs
                
                times = results{1}{1}.times;
                freqs = results{1}{1}.freqs;
                
                for p = 1:size(results{1}{1}.probeindex,1) %iterate through probechans
                    
                    coh{t}{p} = [];
                    
                    for e = 1:length(results{g}) %epochs
                        %combine events from all eps
                        if results{g}{e}.ripsizes(1)>0
                            coh{t}{p} = cat(3,coh{t}{p},results{g}{e}.coherence{t}{p});
                        end
                    end
                    numevents = size(coh{t}{p},3);
                    meancoh{t}{p} = mean(coh{t}{p},3);
                    
                end
                
                %plot all RTS
                shank = size(results{1}{1}.probeindex,1);
                
                if shank <31
                    %two shank
                    set(gcf,'Position',[1 1 3095 364]);
                    disp('twosh')
                    for s = 1:shank
                        chanpos = results{1}{1}.probeindex(s,3);
                        
                        if chanpos<3
                            subplot(2,16,chanpos)
                        else
                            subplot(2,16,chanpos+1)
                        end
                        plot(times, mean(meancoh{t}{s}(:,3:7),2));
                        %set(gca, 'Clim', [.4 .7]);
                        %set(gca,'YDir','normal');
                        details = sprintf('%s %s', chinfo{1}{1}{chanpos}.area , chinfo{1}{1}{chanpos}.layer) ;
                        title(details, 'FontSize',10,'Fontweight','normal');
                        %ylabel('Freq','FontSize',10,'Fontweight','normal');
                        %xlabel('Time(s)','FontSize',10,'Fontweight','normal');
                        set(gca,'XLim',[min(times) max(times)]);
                        set(gca,'YLim',[0 1]);
                        % Plot Line at 0 ms - Start of ripple
                        %hold on
                        %ypts = freqs(1:26);
                        %xpts = 0*ones(size(ypts));
                        %plot(xpts , ypts, 'k--','Linewidth',2);
                        %colorbar;
                        
                    end
                    
                    %
                    
                    suptitle = sprintf('%s coh %drips on ch%d, 5SD during longhc immob>30',f(a).animal{3},numevents,results{g}{t}.trigindex(t,3));
                    [ax,h]=suplabel(suptitle ,'t');
                    tightfig;
                    
                else
                    %four shank
                    set(gcf,'Position',[1 1 1595 964]);
                    disp('four sh')
                    for s = 1:shank
                        chanpos = results{1}{1}.probeindex(s,3);
                        
                        subplot(4,8,chanpos)
                        
                        plot(times,mean(meancoh{t}{s}(:,3:7),2));
                        %set(gca, 'Clim', [.4 .7]);
                        %set(gca,'YDir','normal');
                        details = sprintf('%s %s', chinfo{1}{1}{chanpos}.area, chinfo{1}{1}{chanpos}.layer) ;
                        title(details, 'FontSize',10,'Fontweight','normal');
                        %ylabel('Freq','FontSize',10,'Fontweight','normal');
                        %xlabel('Time(s)','FontSize',10,'Fontweight','normal');
                        set(gca,'XLim',[min(times) max(times)]);
                        set(gca,'YLim',[0 1]);
                        % Plot Line at 0 ms - Start of ripple
                        hold on
                        %ypts = freqs(1:26);
                        %xpts = 0*ones(size(ypts));
                        %plot(xpts , ypts, 'k--','Linewidth',2);
                        %colorbar;
                        
                    end
                    suptitle = sprintf('%s coh %drips on ch%d, 5SD in longhc immob>30',f(a).animal{3},numevents,results{g}{t}.trigindex(t,3));
                    [ax,h]=suplabel(suptitle ,'t');
                    tightfig;
                    
                end
                
                %saveas(gcf, sprintf('Coh longhc %s 5SD immobover30 nooverlap',f(a).animal{3}),'tif')
                
                
            end
        end
        
        
    end
    
    %% calculate  coh in various freq bands at one timepoint for specific areas
    % average of each individual mouse
    
elseif figopt==4
    set(0,'DefaultAxesFontName','Arial')
%         location = '( isequal($area,''dg'') && contains($layer,''val'') && (contains($layer,''*gc*'') || contains($layer,''*hil*'')))';
    location = '( isequal($area,''ca3'')&& (contains($layer,''pyr'') || contains($layer,''sr'')))';
    %     location = '( isequal($area,''ca3'')&& contains($layer,''**mua**''))';
    
    
    %plot power at time in bands:
    bands = [4 7];
    bandnames = {'slow gamma'};
    baseind = 1;  %first time bin
    peakind = 3;%  %time bin immediately after thresh crossing
    
    animnames = [];
    animgeno = [];
    for b = 1:size(bands,1)
        bandstart = bands(b,1);
        bandend = bands(b,2);
        animnames = [];
        animgeno = [];
        indices = [];
        for a = 1:length(f)
            allcoh{a} = [];
            results = f(a).output.RTcoh.results;
            %load chinfo file to lookup site locations
            animinfo = animaldef(f(a).animal{3});
            infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
            load(infofile)
            %identify sites in location of interest
            temp = evaluatefilter(chinfo,location);
            psites{a} = unique(temp(:,3));
            [psites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),psites{a});
            if (~isempty(siteinds) && ~isempty(results{1}{1}.cohindex))
                times = results{1}{1}.times;
                freqs = results{1}{1}.freqs;
                %                 base{a} = repmat(-5000,5000,length(siteinds));
                %                 peak{a} = repmat(-5000,5000,length(siteinds));
                for c = 1:size(results{1}{1}.cohindex,1)
                    for p = 1:length(siteinds)  %iterate thru selected sites
                        lastind=0;
                        %combine data over epochs
                        basecoh{a}{c}{p} = [];
                        peakcoh{a}{c}{p} = [];
                        delta{a}{c}{p} = [];
                        meancoh{a}{c}{p} = [];
                        normcoh{a}{c}{p} = [];
                        deltamean{a}{c}{p} = [];
                        
                        for e = 1:length(results{1})
                            [psites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),psites{a});
                            if length(results{1}{e}.freqs)>1 && p<=length(siteinds)
                                [psites{a},siteinds] = intersect(results{1}{e}.probeindex(:,3),psites{a});
                                basec = squeeze(nanmean(results{1}{e}.coherence{c}{siteinds(p)}(baseind,bandstart:bandend,:)));
                                peakc = squeeze(nanmean(results{1}{e}.coherence{c}{siteinds(p)}(peakind,bandstart:bandend,:)));
                                basecoh{a}{c}{p} = [basecoh{a}{c}{p}; basec];
                                peakcoh{a}{c}{p} = [peakcoh{a}{c}{p}; peakc];
                                delta{a}{c}{p} = [delta{a}{c}{p}; peakc-basec];
                                meanc = nanmean(results{1}{e}.meanCoh{c}{siteinds(p)}(bandstart:bandend));
                                stdc = nanmean(results{1}{e}.stdCoh{c}{siteinds(p)}(bandstart:bandend));
                                meancoh{a}{c}{p} = [meancoh{a}{c}{p}; meanc];
                                deltamean{a}{c}{p} = [deltamean{a}{c}{p}; peakc-meanc];
                                normcoh{a}{c}{p} = [normcoh{a}{c}{p}; (peakc-meanc)/stdc];
                            end
                        end
                        chanbasecoh{a}(c,p) = mean(basecoh{a}{c}{p});
                        chanpeakcoh{a}(c,p) = mean(peakcoh{a}{c}{p});
                        chandelta{a}(c,p) = mean(delta{a}{c}{p});
                        chanmeancoh{a}(c,p) = mean(meancoh{a}{c}{p});
                        channormcoh{a}(c,p) = mean(normcoh{a}{c}{p});
                        chandeltamean{a}(c,p) = mean(deltamean{a}{c}{p});
                        
                    end
                    
                    trigcoh{a}{c} = cell2mat(delta{a}{c});
                    allcoh{a} = [allcoh{a} trigcoh{a}{c}];
                    
                end
                avgcoh{a} = mean(allcoh{a},2);
                allcoh{a} = allcoh{a}(:);
                
                
            end
        end
    end
    
    
elseif figopt==5  % sort based on slow gamma power to normalize
    
    set(0,'DefaultAxesFontName','Arial')
    clearvars -except f figopt
    
    %location = '(isequal($area,''ca3''))';
    location = '( isequal($area,''dg'')&& ~contains($layer,''low''))';
    %location = '( isequal($area,''ca1'')&& contains($layer,''sr 1''))'; %
    animnames = [];
    animgeno = [];
    cols = [];
    for a = 1:length(f)
        results = f(a).output.RTcoh.results;
        %load chinfo file to lookup site locations
        infofile = sprintf('%s%schinfo.mat',f(a).animal{2},f(a).animal{3});
        load(infofile)
        %identify sites in location of interest
        temp = evaluatefilter(chinfo,location);
        psites{a} = unique(temp(:,3));
        [psites{a},siteinds] = intersect(results{1}{1}.probeindex(:,3),psites{a});
        
        for p = 1:length(siteinds)  %iterate thru selected sites
            %combine data over epochs
            probe{a}{p} = [];
            times = results{1}{1}.times;
            freqs = results{1}{1}.freqs;
            for e = 1:length(results{1})
                if results{1}{e}.riptimes(1)>0
                    probe{a}{p} = cat(3,probe{a}{p},results{1}{e}.coherence{1}{siteinds(p)});
                end
            end
        end
        try
            numev(a) = length(probe{a}{p});
        catch
            numev(a) = 0;
        end
        animnames = [animnames; f(a).animal(3) ];
        cols = [cols; f(a).animal{5} ];
        animgeno = [animgeno; f(a).animal(4)];
    end
    
    
    
end