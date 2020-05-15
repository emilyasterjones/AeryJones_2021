function plotsubfig(vehdata,cnodata,color,feature,subfig,limits,stars)
bar(-1,0,0.65)
hold on
vehpoints = cellfun(@nanmean,vehdata);
cnopoints = cellfun(@nanmean,cnodata);
fprintf('%s: %g %g\n',feature,min([vehpoints cnopoints]),max([vehpoints cnopoints]));
bar(1,nanmean(vehpoints),0.65,'FaceAlpha',0,'LineWidth',1.25,'EdgeColor',color)
bar(2,nanmean(cnopoints),0.65,'FaceColor',color,'LineWidth',1.25,'EdgeColor',color)
plot([repmat(1,length(vehdata),1)'; repmat(2,length(vehdata),1)'],[vehpoints; cnopoints],'Color',[0.5 0.5 0.5],'LineWidth',0.75)
ylim(limits)
yt = ylim;
text(1.5,yt(1),stars,'FontWeight','bold','HorizontalAlignment','center')
set(gca,'TickDir','out','Color','none','LineWidth',1,'FontSize',7,...
    'XColor','k','YColor','k','TickLength',[0.02 0])
xlim([0.2 2.8])
ylabel(feature,'FontSize',8)
text(-1.25,yt(1),subfig,'FontSize',11,'FontWeight','bold','VerticalAlignment','bottom')
box off
end