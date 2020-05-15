load('All_Figures.mat');
% requires data structure with the following:
% structs pv, sst, pvsst, and ev, with fields comprised of (n_animalsx1)
% cells of all data, 1 field for each feature plotted
% strings feature (y-axis labels) and subfig (figure panel identifier)
% doubles limits (y limits for most plots) and delta_limits
% (y limits for Figures 6 and S5)
% structs fm (rest features) and fm2 (run features) containing a logical
% for each genotype (n_animalsx1, female=0, males=1)

figind = 1:15; %1:15 for REST features, 16:length(subfig) for RUN features
logtransform = 1;

%% Figures 2, 3, and S2 OR 3, 4, S3, and S4, plus Table S2
p = zeros(length(figind),3);
df = zeros(length(figind),3);
F = zeros(length(figind),3);
for i = figind
    % treatment effects within genotypes
    eval(['lme = calcLMM(pv.veh',char(subfig(i)),',pv.cno',char(subfig(i)),',logtransform);'])
    [p(i,1),F(i,1),~,df(i,1)] = coefTest(lme);
    eval(['lme = calcLMM(sst.veh',char(subfig(i)),',sst.cno',char(subfig(i)),',logtransform);'])
    [p(i,2),F(i,2),~,df(i,2)] = coefTest(lme);
    
    % treatment-genotype interaction effects
    vehpv = eval(['pv.veh',char(subfig(i))]);
    cnopv = eval(['pv.cno',char(subfig(i))]);
    vehsst = eval(['sst.veh',char(subfig(i))]);
    cnosst = eval(['sst.cno',char(subfig(i))]);
    [p(i,3),F(i,3),df(i,3)] = calcInteraction(vehpv, cnopv, vehsst, cnosst, logtransform);
end

% print legends
for i = figind
    fprintf('PV: p = %.2g; SST: p = %.2g; PV vs SST: p = %.2g\n',p(i,1),p(i,2),p(i,3))
end

% correct for multiple comparisons
stars = HBcorrect(p);

% print statistical details for Table S2
for i = figind
    fprintf('%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%s\n',df(i,1),F(i,1),df(i,2),F(i,2),df(i,3),F(i,3),stars(i,3))
end

% Plot
figure
for i = figind
    subplot(4,5,mod(i-1,10)*2+1)
    eval(['plotsubfig(pv.veh',char(subfig(i)),',pv.cno',char(subfig(i)),...
        ',[55/255 116/255 184/255],feature(i),subfig(i),limits(i,:),stars(i,1));'])
    subplot(4,5,mod(i-1,10)*2+2)
    eval(['plotsubfig(sst.veh',char(subfig(i)),',sst.cno',char(subfig(i)),...
        ',[231/255 172/255 80/255],feature(i),subfig(i),limits(i,:),stars(i,2));'])
    if ~mod(i,10) || i == figind(end)
        set(gcf, 'PaperUnits', 'inches', 'PaperSize', [8.5 11],...
            'PaperPositionMode', 'manual', 'PaperPosition', [0.5 0.5 7.5 10]);
        saveas(gcf,sprintf('Plots%d.pdf',ceil(i/10)))
        figure
    end
end

%% Figures 6, S5, & Table S4
p = zeros(length(figind),4);
df = zeros(length(figind),4);
F = zeros(length(figind),4);

for i = figind
    % treatment effects within genotypes
    eval(['lme = calcLMM(pvsst.veh',char(subfig(i)),',pvsst.cno',char(subfig(i)),',logtransform);'])
    [p(i,3),F(i,3),~,df(i,3)] = coefTest(lme);
    eval(['lme = calcLMM(ev.veh',char(subfig(i)),',ev.cno',char(subfig(i)),',logtransform);'])
    [p(i,4),F(i,4),~,df(i,4)] = coefTest(lme);
    
    % treatment-genotype interaction effects
    vehpv = eval(['pv.veh',char(subfig(i))]);
    cnopv = eval(['pv.cno',char(subfig(i))]);
    vehsst = eval(['sst.veh',char(subfig(i))]);
    cnosst = eval(['sst.cno',char(subfig(i))]);
    vehpvsst = eval(['pvsst.veh',char(subfig(i))]);
    cnopvsst = eval(['pvsst.cno',char(subfig(i))]);
    [p(i,1),F(i,1),df(i,1)] = calcInteraction(vehpv, cnopv, vehpvsst, cnopvsst, logtransform);
    [p(i,2),F(i,2),df(i,2)] = calcInteraction(vehsst, cnosst, vehpvsst, cnopvsst, logtransform);
end

% print legends
for i = figind
    fprintf('PV/SST: p = %.2g; PV vs PV/SST: p = %.2g; SST vs PV/SST: p = %.2g; EV: p = %.2g\n',...
        p(i,3),p(i,1),p(i,2),p(i,4))
end

% correct for multiple comparisons
stars = HBcorrect(p);

% print statistical details for Table S2
for i = figind
    fprintf('%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n',df(i,3),F(i,3),df(i,1),F(i,1),...
        df(i,2),F(i,2),df(i,4),F(i,4))
end

% Plot
figure
groups = {'pv','sst','pvsst','ev'};
colors = [40 116 200; 231 172 80; 80 175 87; 128 128 128];
colors = colors/256;
for i = figind
    subplot(5,3,mod(i-1,15)+1)
    line([0 5],[0 0],'Color','k','LineStyle',':','LineWidth',1)
    for g = 1:length(groups)
        eval(['[lme, vehcoeff, cnocoeff] = calcLMM(',groups{g},'.veh',char(subfig(i)),','...
            ,groups{g},'.cno',char(subfig(i)),',logtransform);'])
        
        if logtransform
            beta = exp(sum(fixedEf)) - exp(fixedEf(1));
            fixedCI = coefCI(lme);
            betaCI = beta - (exp(sum(fixedCI)) - exp(fixedCI(1,:)));
            betaCI2 = exp(sum(fixedCI)) - exp(fixedCI(1,:));
        else
            fixedEf = fixedEffects(lme);
            beta = fixedEf(2);
            fixedCI = coefCI(lme);
            betaCI2 = fixedCI(2,:);
        end
        
        hold on
        patch([g-.2 g+.2 g+.2 g-.2], [betaCI2(1) betaCI2(1) betaCI2(2) betaCI2(2)], ...
            colors(g,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        points = plotSpread(cnocoeff-vehcoeff,'spreadWidth',2);
        scatter(points(:,1)+(g-1), points(:,2), 10, colors(g,:), 'filled')
        line([g-0.2 g+0.2], [beta beta], 'Color', 'k', 'LineWidth', 1)
        yt = ylim;
        text(g,yt(1),groups(g),'FontSize',8,'HorizontalAlignment','center')
        text(g,yt(1),stars(i,g),'FontWeight','bold','HorizontalAlignment','center')
        
    end
    xlim([0.5 4.5])
    ylim([delta_limits(i,1) delta_limits(i,2)])
    ylabel(feature(i),'FontSize',8)
    text(0,delta_limits(i,1),subfig(i),'FontSize',11,'FontWeight','bold','VerticalAlignment','bottom')
    set(gca,'TickDir','out','Color','none','LineWidth',1,'FontSize',7,...
        'XColor','none','YColor','k','TickLength',[0.02 0])
    line([0 5],[delta_limits(i,1) delta_limits(i,1)],'Color','k','LineWidth',1)
    
    if ~mod(i,15) || i == figind(end)
        set(gcf, 'PaperUnits', 'inches', 'PaperSize', [8.5 11],...
            'PaperPositionMode', 'manual', 'PaperPosition', [0.5 0.5 7.5 10]);
        saveas(gcf,sprintf('Deltaplots%d.pdf',ceil(i/15)))
        figure
    end
end


%% Table S3
p = zeros(length(figind),1);
df = zeros(length(figind),1);
F = zeros(length(figind),1);
pvmean = zeros(length(figind),1);
sstmean = zeros(length(figind),1);
for i = figind
    eval(['lme = calcunpairedLMM(pv.veh',char(subfig(i)),', sst.veh',char(subfig(i)),', logtransform);'])
    [p(i,1),F(i,1),~,df(i,1)] = coefTest(lme);
    eval(['pvmean(i,1) = mean(cellfun(@mean, pv.veh',char(subfig(i)),'));'])
    eval(['sstmean(i,1) = mean(cellfun(@mean, sst.veh',char(subfig(i)),'));'])
end
stars = HBcorrect(p);
for i = figind
    fprintf('%s\t%.2g\t%.2g\t%d\t%.2g\t%.2g\t%s\n',feature(i),...
        pvmean(i,1),sstmean(i,1),df(i,1),F(i,1),p(i,1),stars(i))
end


%% Table S5
% use fm for RUN features, fm2 for REST features
p = zeros(length(figind),1);
df = zeros(length(figind),1);
F = zeros(length(figind),1);
femalesmean = zeros(length(figind),1);
malesmean = zeros(length(figind),1);
for i = figind
    males = eval(['[pv.veh',char(subfig(i)),'(fm.pv) sst.veh',...
        char(subfig(i)),'(fm.sst) pvsst.veh',char(subfig(i)),...
        '(fm.pvsst) ev.veh',char(subfig(i)),'(fm.ev)];']);
    females = eval(['[pv.veh',char(subfig(i)),'(~fm.pv) sst.veh',...
        char(subfig(i)),'(~fm.sst) pvsst.veh',char(subfig(i)),...
        '(~fm.pvsst) ev.veh',char(subfig(i)),'(~fm.ev)];']);
    lme = calcunpairedLMM(males, females, logtransform);
    [p(i,1),F(i,1),~,df(i,1)] = coefTest(lme);
    malesmean(i,1) = mean(cellfun(@mean, males));
    femalesmean(i,1) = mean(cellfun(@mean, females));
end
stars = HBcorrect(p);
for i = figind
    fprintf('%s\t%.2g\t%.2g\t%d\t%.2g\t%.2g\t%s\n',feature(i),...
        femalesmean(i,1),malesmean(i,1),df(i,1),F(i,1),p(i,1),stars(i))
end


%% Table S6
% use fm for RUN features, fm2 for REST features
p = zeros(length(figind),3);
df = zeros(length(figind),3);
F = zeros(length(figind),3);
femalesbeta = zeros(length(figind),3);
malesbeta = zeros(length(figind),3);
groups = {'pv','sst','pvsst'};
for g = 1:length(groups)
    groups{g}
    for i = figind
        vehmales = eval([groups{g},'.veh',char(subfig(i)),'(fm.',groups{g},')']);
        cnomales = eval([groups{g},'.cno',char(subfig(i)),'(fm.',groups{g},')']);
        vehfemales = eval([groups{g},'.veh',char(subfig(i)),'(~fm.',groups{g},')']);
        cnofemales = eval([groups{g},'.cno',char(subfig(i)),'(~fm.',groups{g},')']);
        
        lme = calcLMM(vehmales, cnomales, logtransform);
        fixedEf = fixedEffects(lme);
        if logtransform
            malesbeta(i,g) = exp(sum(fixedEf)) - exp(fixedEf(1));
        else
            malesbeta(i,g) = fixedEf(2);
        end
        
        lme = calcLMM(vehfemales, cnofemales, logtransform);
        fixedEf = fixedEffects(lme);
        if logtransform
            femalesbeta(i,g) = exp(sum(fixedEf)) - exp(fixedEf(1));
        else
            femalesbeta(i,g) = fixedEf(2);
        end
        
        [p(i,g),F(i,g),df(i,g)] = calcInteraction(vehmales, cnomales, vehfemales, cnofemales, logtransform);
    end
    stars = HBcorrect(p(:,g));
    for i = figind
        fprintf('%s\t%.2g\t%.2g\t%d\t%.2g\t%.2g\t%s\n',feature(i),...
            femalesbeta(i,g),malesbeta(i,g),df(i,g),F(i,g),p(i,g),stars(i,1))
    end
end