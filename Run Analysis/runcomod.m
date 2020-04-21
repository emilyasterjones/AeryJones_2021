clearvars -except f location

% pre-processing parameters
% version 2: low freq combination
PhaseFreqVector=0:.5:12;
AmpFreqVector=0:5:150;
PhaseFreq_BandWidth=0.5;
AmpFreq_BandWidth=5;

% version 3: single bands
% PhaseFreqVector=5;
% AmpFreqVector=20;
% PhaseFreq_BandWidth=6;
% AmpFreq_BandWidth=30;

% Animal
animals = {'Adobo_LT', 'Caraway_LT', 'Chives_LT', ...
    'Cinnamon_LT', 'Coriander_LT', 'Fenugreek_LT', 'GaramMasala_LT', 'Salt_LT', ...
    'Baharat_LT', 'Cardamom_LT', 'Jerk_LT', 'Mace_LT', 'Mustard_LT', ...
    'Tarragon_LT', 'Vanilla_LT',...
    'Basil_LT', 'Cumin_LT', 'Dill_LT', 'Nutmeg_LT', 'Paprika_LT', 'Parsley_LT', 'Sumac_LT',...
    'Pepper_LT', 'Sage_LT', 'Anise_LT', 'Thyme_LT', 'OldBay_LT', 'Rosemary_LT',...
    'Provence_LT', 'Saffron_LT'};

epochfilter = [];
epochfilter{1} = {'task','(contains($env,''LT'') && ~contains($descript,''FAIL'') && isequal($treatment,''veh''))'};

datafilter = [];  % load all channels...use embedded filter in run function to get specific ones
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
timefilter{1} =  {'pos', '($vel > 1)'};

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
%Specify function and run
f = setfilterfunction(f, 'comod', {'eeg'},'PhaseFreqVector',PhaseFreqVector,'AmpFreqVector',AmpFreqVector,'PhaseFreq_BandWidth',PhaseFreq_BandWidth,'AmpFreq_BandWidth',AmpFreq_BandWidth);
f = runfilter(f);

figopt = 3;

%% Plot
%figopt1 = full spectrum
%figopt2 = low gamma zoom
%figopt3 = quantification
%figopt4 = peak freq of modulation
%figopt5 = single site heatmap

%all chans by location, averaged over eps.  Full spectrum
if figopt==1
    
    for a = 1:length(f)
        
        results = f(a).output.comod.results{1}; %condition 1 only
        
        %combine data over epochs for each channel
        
        for c = 1:length(results{1}.comodulogram) %iterate through channels
            %iterate thru probe chans, combine eps, plot
            totaltime = 0;   %same for all channels
            comodgram{c} = [];
            for e = 1:length(results) %epochs
                %combine events from all eps
                comodgram{c} = cat(3,comodgram{c},results{e}.comodulogram{c});
                totaltime = totaltime+results{e}.totalsecofeeg{c};
            end
            meancomod{c} = mean(comodgram{c},3);
        end
        
        %plot all RTS
        figure;
        phaseaxis = results{1}.frequency_phase + results{1}.phasefreq_bandwidth/2;
        ampaxis = results{1}.frequency_amplitude + results{1}.ampfreq_bandwidth/2;
        shank = results{1}.index.chinfo(end);
        
        if shank <32 %two shank
            
            set(gcf,'Position',[1 1 3095 364]);
            
            for s = 1:length(results{1}.index.chinfo)
                chan = results{1}.index.chinfo(s);
                
                if chan<3
                    subplot(2,16,chan)
                else
                    subplot(2,16,chan+1)
                end
                
                contourf(phaseaxis, ampaxis,meancomod{s}',30,'lines','none')
                caxis([0 .002]);
                title(chan)
                if chan==shank
                    colorbar
                end
                
            end
            suptitle = sprintf('%s comod avg over %d sec, vel>3',f(a).animal{3},int16(totaltime));
            [ax,h]=suplabel(suptitle ,'t');
            tightfig;
            
        else %four shank
            
            set(gcf,'Position',[1 1 1595 964]);
            
            for s = 1:length(results{1}.index.chinfo)
                chan = results{1}.index.chinfo(s);
                
                subplot(4,8,chan)
                contourf(phaseaxis, ampaxis,meancomod{s}',30,'lines','none')
                %                 caxis([0 .002]);
                title(chan)
                if chan==shank
                    colorbar
                end
                
            end
            suptitle = sprintf('%s comod avged over %d seconds, vel>3',f(a).animal{3},int16(totaltime));
            [ax,h]=suplabel(suptitle ,'t');
            tightfig;
        end
    end
end
%%

if figopt==2  %low gamma zoom
    
    for a = 1:length(f)
        
        results = f(a).output.comod.results{1}; %condition 1 only
        
        %combine data over epochs for each channel
        
        for c = 1:length(results{1}.comodulogram) %iterate through channels
            %iterate thru probe chans, combine eps, plot
            totaltime = 0;   %same for all channels
            comodgram{c} = [];
            for e = 1:length(results) %epochs
                %combine events from all eps
                comodgram{c} = cat(3,comodgram{c},results{e}.comodulogram{c});
                totaltime = totaltime+results{e}.totalsecofeeg{c};
            end
            meancomod{c} = mean(comodgram{c},3);
        end
        
        
        %plot all RTS
        figure;
        phaseaxis = results{1}.frequency_phase + results{1}.phasefreq_bandwidth/2;
        ampaxis = results{1}.frequency_amplitude + results{1}.ampfreq_bandwidth/2;
        shank = results{1}.index.chinfo(end);
        
        if shank <32 %two shank
            
            set(gcf,'Position',[1 1 3095 364]);
            
            for s = 1:length(results{1}.index.chinfo)
                chan = results{1}.index.chinfo(s);
                
                if chan<3
                    subplot(2,16,chan)
                else
                    subplot(2,16,chan+1)
                end
                
                contourf(phaseaxis, ampaxis(1:10),meancomod{s}(:,1:10)',30,'lines','none')
                
                caxis([0 .0005]);
                %set(gca,'YDir','normal');
                %colormap('jet')
                title(chan)
                if chan==shank
                    colorbar
                end
                
            end
            suptitle = sprintf('%s comod avg over %d sec, vel>3 LOW GAMMA',f(a).animal{3},int16(totaltime));
            [ax,h]=suplabel(suptitle ,'t');
            tightfig;
            
        else %four shank
            
            set(gcf,'Position',[1 1 1595 964]);
            
            for s = 1:length(results{1}.index.chinfo)
                chan = results{1}.index.chinfo(s);
                
                subplot(4,8,chan)
                contourf(phaseaxis, ampaxis(1:10),meancomod{s}(:,1:10)',30,'lines','none')
%                 caxis([0 .0005]);
                %set(gca,'YDir','normal');
                %colormap('jet')
                title(chan)
                if chan==shank
                    colorbar
                end
                
            end
            suptitle = sprintf('%s comod avged over %d seconds, vel>3 LOW GAMMA',f(a).animal{3},int16(totaltime));
            [ax,h]=suplabel(suptitle ,'t');
            tightfig;
        end
    end
end

% calculate comod index & phase of peak modulation
% only works for a single theta & gamma frequency band
if figopt == 3
    
    location = '(isequal($area,''ca1'') && contains($layer,''pyr''))';
    % location = '(isequal($area,''ca1'') && contains($layer,''sr''))';
    % location = '(isequal($area,''ca1'') && contains($layer,''slm''))';
    % location = '(isequal($area,''ca3'') && (contains($layer,''pyr'') | contains($layer,''sr'')))';
    % location = '( isequal($area,''dg'')&& contains($layer,''mol''))';
    % location = '( isequal($area,''dg'') && (contains($layer,''gc'')||contains($layer,''hil'')))';
    
    format compact
    format shortG
    clearvars coh modovercycle maxmod maxphase
    
    nbin = 18;
    center = (-pi+pi/(2*nbin)):pi/nbin:(pi-pi/(2*nbin));
    
    % calculate comod index
    for a = 1:length(f)
        coh{a} = [];
        allcomod{a} = [];
        results = f(a).output.comod.results;
        %load chinfo file to lookup site locations
        animinfo = animaldef(f(a).animal{1});
        infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
        load(infofile)
        
        %identify sites in location of interest
        pwrtemp = evaluatefilter(chinfo,location);
        sites{a} = unique(pwrtemp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.index.chinfo,sites{a}); %probeindex(:,3)
        for c = 1:length(siteinds)
            allcomod{a}{c} = [];
            for e = 1:length(results{1})
                [sites{a},siteinds] = intersect(results{1}{e}.index.chinfo,sites{a});
                if length(results{1}{e}.comodulogram)>1 && c<=length(siteinds)
                    coh{a}(c,e) = results{1}{e}.comodulogram{siteinds(c)};
                    allcomod{a}{c} = [allcomod{a}{c}; results{1}{e}.comodulogram_bin{siteinds(c)}];
                    modovercycle{a}(c,e,:) = results{1}{e}.comodbyphase{siteinds(c)};
                    [mx, bin] = max(results{1}{e}.comodbyphase{siteinds(c)});
                    maxmod{a}(c,e) = mx;
                    maxphase{a}(c,e) = center(bin+1);
                end
            end
        end
        meancoh(a) = mean(mean(coh{a}));
        meanmodovercycle(a,:) = mean(mean(modovercycle{a}));
        meanmaxmod(a) = mean(mean(maxmod{a}));
        meanmaxphase(a) = mean(mean(maxphase{a}));
        
        allcomod{a} = cell2mat(allcomod{a});
        avgcomod{a} = nanmean(allcomod{a},2);
        allcomod{a} = allcomod{a}(:);
    end
    
    % meancoh'
    % meanmaxmod'
    % meanmaxphase'
    %
    % figure
    % plotdata = meanmodovercycle(9:15,:);
    % plot(center(1:end-1),nanmean(plotdata),'Color',[0.543 0 0],'Linewidth',1) %,[0.543 0 0]
    % hold on
    % sem = nanstd(plotdata)/sqrt(size(plotdata,2));
    % h = fill([center(1:end-1) fliplr(center(1:end-1))], [(nanmean(plotdata)-sem) fliplr(nanmean(plotdata)+sem)],[0.543 0 0],'FaceAlpha',.5); %[.543 0 0]
    % set(h,'EdgeColor','none')
    % title(location)
    % box off
    % xlim([-pi pi])
    % xticks([-pi -pi/2 0 pi/2 pi])
    % xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    
end

% over theta/gamma freq bins, find freq of max mod
% either phase or amp freq bins (set manually) must be 1 bin
if figopt==4
    
    %     freqs = 5.25:0.5:10.75;
    
    %     location = '(isequal($area,''ca1'') && contains($layer,''*pyr*''))';
    % location = '(isequal($area,''ca1'') && contains($layer,''*sr*''))';
    % location = '(isequal($area,''ca1'') && contains($layer,''*slm*''))';
    % location = '(isequal($area,''ca3'') && (contains($layer,''*pyr*'') | contains($layer,''*sr*'')))';
    % location = '( isequal($area,''dg'')&& contains($layer,''mol*''))';
    % location = '( isequal($area,''dg'') && (contains($layer,''*gc*'')||contains($layer,''*hil*'')))';
    
    format compact
    format shortG
    clearvars maxfreq
    
    % calculate comod index
    for a = 1:length(f)
        coh{a} = [];
        results = f(a).output.comod.results;
        %load chinfo file to lookup site locations
        animinfo = animaldef(f(a).animal{1});
        infofile = sprintf('%s%schinfo.mat',animinfo{2},animinfo{3});
        load(infofile)
        
        %identify sites in location of interest
        pwrtemp = evaluatefilter(chinfo,location);
        sites{a} = unique(pwrtemp(:,3));
        [sites{a},siteinds] = intersect(results{1}{1}.index.chinfo,sites{a}); %probeindex(:,3)
        for c = 1:length(siteinds)
            for e = 1:length(results{1})
                [sites{a},siteinds] = intersect(results{1}{e}.index.chinfo,sites{a});
                if length(results{1}{e}.comodulogram)>1 && c<=length(siteinds)
                    [mx, ind] = max(results{1}{e}.comodulogram{siteinds(c)});
                    maxfreq{a}(c,e) = freqs(ind);
                end
            end
        end
        meanmaxfreq(a) = mean(mean(maxfreq{a}));
    end
    
    %     meanmaxfreq'
    
    %
%     figure
%     plotdata = meanmodovercycle(9:15,:);
%     plot(center(1:end-1),nanmean(plotdata),'Color','k','Linewidth',1) %,[0.543 0 0]
%     hold on
%     sem = nanstd(plotdata)/sqrt(size(plotdata,2));
%     h = fill([center(1:end-1) fliplr(center(1:end-1))], [(nanmean(plotdata)-sem) fliplr(nanmean(plotdata)+sem)],'k','FaceAlpha',.5); %[.543 0 0]
%     set(h,'EdgeColor','none')
%     title(location)
%     box off
    
    
end

if figopt==5
    
    for a = 1:length(f)
        
        results = f(a).output.comod.results{1}; %condition 1 only
        
        %combine data over epochs for each channel
        
        for c = 1:length(results{1}.comodulogram) %iterate through channels
            %iterate thru probe chans, combine eps, plot
            totaltime = 0;   %same for all channels
            comodgram{c} = [];
            for e = 1:length(results) %epochs
                %combine events from all eps
                comodgram{c} = cat(3,comodgram{c},results{e}.comodulogram{c});
                totaltime = totaltime+results{e}.totalsecofeeg{c};
            end
            meancomod{c} = mean(comodgram{c},3);
        end
        
        phaseaxis = results{1}.frequency_phase + results{1}.phasefreq_bandwidth/2;
        ampaxis = results{1}.frequency_amplitude + results{1}.ampfreq_bandwidth/2;
        
        for s = 1:length(results{1}.index.chinfo)
            figure
            %         s=15;
            chan = results{1}.index.chinfo(s);
            imagesc(phaseaxis, ampaxis,meancomod{s}');
            %         caxis([0 .002]);
            colormap hot %jet
            set(gca,'YDir','normal');
            title(chan)
            colorbar
        end
    end
end