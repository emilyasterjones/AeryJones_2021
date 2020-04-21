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
epochfilter{1} = {'task','(contains($env,''LT'') && isequal($treatment, ''veh'') && ~contains($descript,''FAIL''))'};

datafilter = [];  % load all channels...use embedded filter in run function to get specific ones
datafilter{1} = {'chinfo','(~isequal($area,''dead''))'};

timefilter = [];
timefilter{1} = {'pos', '($vel>1)'};

% set channel filtering criteria
cohcrit = '(isequal($area,''ca1'') && contains($layer,''pyr''))';
probecrit = '( isequal($area,''ca3'')&& (contains($layer,''pyr'') || contains($layer,''sr'')))';
% probecrit = '( isequal($area,''dg'') && (contains($layer,''gc'')||contains($layer,''hil'')))';

% run filters
f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);
f = setfilterfunction(f, 'calccoherence', {'eeg','chinfo'}, 'coh',cohcrit,'probe',probecrit);
f = runfilter(f);

figopt=1;

%% Quantify
% figopts:
% 1 = coherence & freq of max coherence in specific freq bands
% 2 = coherence across all freqs

if figopt==1
    %frequency bands
    results = f(1).output.calccoherence.results;
    theta = [5 11];
    slowgamma = [20 50];
    fastgamma = [50 110];
    bandnames = {'theta', 'slow gamma', 'fast gamma'};
    
    meancoh = NaN(length(f),length(bandnames));
    stdcoh = NaN(length(f),length(bandnames));
    maxcohfreq = NaN(length(f),length(bandnames));
    clearvars emeancoh estdcoh emaxcoh
    
    for b = 1:length(bandnames)
        for a = 1:length(f)
            results = f(a).output.calccoherence.results;
            for e = 1:length(results{1})
                emeancoh{e} = [];
                estdcoh{e} = [];
                emaxcohfreq{e} = [];
                if ~isempty(results{1}{e}.cohindex) && ~isempty(results{1}{e}.probeindex)
                    freqs = results{1}{e}.freqs;
                    bands = [lookup(theta(1),freqs) lookup(theta(2),freqs); ...
                        lookup(slowgamma(1),freqs) lookup(slowgamma(2),freqs); ...
                        lookup(fastgamma(1),freqs) lookup(fastgamma(2),freqs)];
                    bandstart = bands(b,1);
                    bandend = bands(b,2);
                    
                    bin_freqs = results{1}{e}.freqs_bin{1}{1};
                    bin_bands = [lookup(theta(1),bin_freqs) lookup(theta(2),bin_freqs); ...
                        lookup(slowgamma(1),bin_freqs) lookup(slowgamma(2),bin_freqs); ...
                        lookup(fastgamma(1),bin_freqs) lookup(fastgamma(2),bin_freqs)];
                    bin_bandstart = bin_bands(b,1);
                    bin_bandend = bin_bands(b,2);
                    
                    for c = 1:length(results{1}{e}.cohindex(:,3))
                        for p = 1:length(results{1}{e}.probeindex(:,3))
                            if length(results{1}{e}.freqs)>1 && p<=length(results{1}{e}.probeindex(:,3)) && c<=length(results{1}{e}.cohindex(:,3))
                                emeancoh{e}(c,p) = nanmean(results{1}{e}.coherence{c}{p}(bandstart:bandend));
                                estdcoh{e}(c,p) = nanstd(results{1}{e}.coherence{c}{p}(bandstart:bandend));
                                [mx, bin] = max(results{1}{e}.coherence{c}{p}(bandstart:bandend));
                                bandfreqs = freqs(bandstart:bandend);
                                emaxcohfreq{e}(c,p) = bandfreqs(bin);
                                
                                nestedcoh{a}{c}{p}{e} = nanmean(results{1}{e}.coherence_bin{c}{p}(:,bin_bandstart:bin_bandend),2);
                                [~, bin] = max(results{1}{e}.coherence_bin{c}{p}(:,bin_bandstart:bin_bandend),[],2);
                                bandfreqs = bin_freqs(bin_bandstart:bin_bandend);
                                nestedfreqs{a}{c}{p}{e} = bandfreqs(bin)';
                            end
                        end
                    end
                    emeancoh{e} = mean(mean(emeancoh{e})); %average over channel pairs
                    estdcoh{e} = mean(mean(estdcoh{e}));
                    emaxcohfreq{e} = mean(mean(emaxcohfreq{e}));
                end
            end
            
            for c = 1:length(nestedcoh{a})
                nprobe = length(nestedcoh{a}{c});
                for p = 1:length(nestedcoh{a}{c})
                    epoch = [];
                    epoch_freqs = [];
                    for e = 1:length(nestedcoh{a}{c}{p})
                        epoch = [epoch; nestedcoh{a}{c}{p}{e}];
                        epoch_freqs = [epoch_freqs; nestedfreqs{a}{c}{p}{e}];
                    end
                    flatcoh{a}(:,nprobe*(c-1)+p) = epoch;
                    flatfreqs{a}(:,nprobe*(c-1)+p) = epoch_freqs;
                end
            end
                
            avgcoh{a} = nanmean(flatcoh{a},2);
            allcoh{a} = flatcoh{a}(:);
            avgcohfreqs{a} = nanmean(flatfreqs{a},2);
            allcohfreqs{a} = flatfreqs{a}(:);
            
            meancoh(a,b) = mean(cell2mat(emeancoh));
            stdcoh(a,b) = mean(cell2mat(estdcoh));
            maxcohfreq(a,b) = mean(cell2mat(emaxcohfreq));
        end
        if b==1
            avgthetacohfreqs = avgcohfreqs;
            allthetacohfreqs = allcohfreqs;
        elseif b==2
            avgSGcohfreqs = avgcohfreqs;
            allSGcohfreqs = allcohfreqs;
        else
            avgFGcohfreqs = avgcohfreqs;
            allFGcohfreqs = allcohfreqs;
        end
            
    end
    
elseif figopt==2
    clearvars ecoh echancoh eprobecoh ebincoh
    bins = 5.5:2:109.5;
    meancoh = NaN(length(f),length(bins));
    
    for a = 1:length(f)
        results = f(a).output.calccoherence.results;
        for e = 1:length(results{1})
            if ~isempty(results{1}{e}.cohindex) && ~isempty(results{1}{e}.probeindex)
                
                %set freq bins
                %sometimes freq list is not the same between epochs
                freqs = results{1}{e}.freqs;
                binfreqs = cell(length(bins),1);
                for freqind = 1:length(freqs)
                    [~, idx] = min(abs(bins - freqs(freqind)));
                    binfreqs{idx} = [binfreqs{idx}, freqind];
                end
                
                for c = 1:length(results{1}{e}.cohindex(:,3))
                    for p = 1:length(results{1}{e}.probeindex(:,3))
                        if p<=length(results{1}{e}.probeindex(:,3)) && c<length(results{1}{e}.cohindex(:,3))
                            eprobecoh{p} = results{1}{e}.coherence{c}{p};
                        end
                    end
                    echancoh{c} = mean(eprobecoh{p},2);
                end
                ecoh{e} = mean(echancoh{c},2); %average over channel pairs
                
                ebincoh{e} = NaN(length(bins),1);
                for binind = 1:length(bins)
                    ebincoh{e}(binind) = mean(ecoh{e}(binfreqs{binind}));
                end
            end
        end
        meancoh(a,:) = mean(cell2mat(ebincoh),2);
    end
    
    figure
    plot(bins,nanmean(meancoh),'Color',[0.543 0 0],'Linewidth',1); %,[0.543 0 0]
    hold on
    sem = nanstd(meancoh,0,1)/sqrt(7);
    h = fill([bins fliplr(bins)], [(nanmean(meancoh)-sem) fliplr(nanmean(meancoh)+sem)],[0.543 0 0],'FaceAlpha',.5); %[.543 0 0]
    set(h,'EdgeColor','none')
end