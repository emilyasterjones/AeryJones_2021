% plot group comparisons, used for Figures 2, S2, 3, and S3
% inputs are same data structures as in run_compare_centers

char = 'A':'Z';

figure
s = 1;
for i = 1:length(features)
    for g = 1:length(grpidx)
        stats = compare_centers(vehresults(grpidx{g},i),cnoresults(grpidx{g},i),'paired',1);
        fprintf('%s %s: %s\t%s\n',features(i),groups(g),get_stars(stats.p),stats.resultsstr);
        subplot(4,5,s)
        bar(-1,0,0.65)
        hold on
        bar(1,nanmean(vehresults(grpidx{g},i)),0.65,'FaceAlpha',0,'LineWidth',1.25)
        bar(2,nanmean(cnoresults(grpidx{g},i)),0.65,'FaceAlpha',0,'LineWidth',1.25,'EdgeColor',[0.543 0 0])
        for a = 1:length(grpidx{g})
            idx = grpidx{g}(a);
            plot([1 2],[vehresults(idx,i),cnoresults(idx,i)],'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        end
        ylim(limits(i,:))
        yt = yticks;
        text(1,yt(1),'Veh','FontSize',8,'HorizontalAlignment','center') %labelpos
        text(2,yt(1),'CNO','FontSize',8,'HorizontalAlignment','center') %labelpos
        text(1.5,yt(1),get_stars(stats.p),'FontWeight','bold','HorizontalAlignment','center') %yt(end)-interval/2
        set(gca,'TickDir','out','Color','none','LineWidth',1,'FontSize',7,...
            'XColor','k','YColor','k','TickLength',[0.02 0])
        xlim([0.2 2.8])
        ylabel(features(i),'FontSize',8)
        text(-1.25,yt(1),char(i),'FontSize',11,'FontWeight','bold','VerticalAlignment','bottom')
        box off
        s = s+1;
        
    end
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [8.5 11],...
        'PaperPositionMode', 'manual', 'PaperPosition', [0.5 0.5 7.5 10]);
    saveas(gcf,'Fig2.pdf')
end