% Compare 2 groups using appropriate statistic

% genotypes: N(animals) matrix
% sexes: N(animals) matrix
% vehresults & cnoresults: N(animals) x N(features) matrix
% features: N(features) matrix with list of feature names
% limits: N(features)x2 matrix of x axis limits

groups = unique(genotypes);
for g = 1:length(groups)
    grpidx{g} = [];
    for a = 1:length(genotypes)
        if isequal(genotypes(a),groups(g))
            grpidx{g} = [grpidx{g} a];
        end
    end
end

disp('Genotypes paired t tests');
for i = 1:length(features)
    for g = 1:length(grpidx)
        stats = compare_centers(vehresults(grpidx{g},i),cnoresults(grpidx{g},i),'paired',1);
        fprintf('%s %s: %s\t%s\n',features(i),groups(g),get_stars(stats.p),stats.resultsstr);
    end
end

sgroups = unique(sexes);
for g = 1:length(sgroups)
    sgrpidx{g} = [];
    for a = 1:length(sexes)
        if isequal(sexes(a),sgroups(g))
            sgrpidx{g} = [sgrpidx{g} a];
        end
    end
end

disp('Sexes vehicle t tests')
for i = 1:length(features)
    stats = compare_centers(vehresults(sgrpidx{1},i),vehresults(sgrpidx{2},i));
    fprintf('%s\t%f\t%f\t%s\t%f\t%f\t%f\n',features(i),nanmean(vehresults(sgrpidx{1},i)),...
        nanmean(vehresults(sgrpidx{2},i)),stats.test,stats.df,stats.testval,stats.p);
end

disp('Sexes deltas t tests')
for g = 1:length(groups)
    fidx = intersect(sgrpidx{1}, grpidx{g});
    midx = intersect(sgrpidx{2}, grpidx{g});
    for i = 1:length(features)
        fdelta = cnoresults(fidx,i) - vehresults(fidx,i);
        mdelta = cnoresults(midx,i) - vehresults(midx,i);
        stats = compare_centers(fdelta,mdelta);
        fprintf('%s\t%s\t%f\t%f\t%s\t%f\t%f\t%f\n',groups(g),features(i),nanmean(fdelta),...
            nanmean(mdelta),stats.test,stats.df,stats.testval,stats.p);
    end
end

disp('Vehicle Baselines')
for i = 1:length(features)
    for g = 1:length(grpidx)-1
        for h = g+1:length(grpidx)
            [~,p,~,~] = vartest2(vehresults(grpidx{g},i), vehresults(grpidx{h},i));
            if p<0.01
                fprintf('%s: %s and %s var different\n',features(i),groups(g),groups(h));
            end
        end
    end
    [p,~,mcstats] = anova1(vehresults(:,i),genotypes,'off');
    fprintf('%s veh genotypes 1-way ANOVA p = %s\n',features(i),get_stars(p));
    mcp = multcompare(mcstats,'Display','off');
    for m = 1:length(mcp(:,6))
        if mcp(m,6)<0.05
            fprintf('%s veh genotypes post-hoc t test: %s vs %s p = %s\n',...
                features(i),groups(mcp(m,1)),groups(mcp(m,2)),get_stars(mcp(m,6)));
        end
    end
end

disp('Deltas')
deltas = cnoresults - vehresults;
ncomp = length(grpidx)*(length(grpidx)-1)/2;
for i = 1:length(features)
    for g = 1:length(grpidx)-1
        for h = g+1:length(grpidx)
            stats = compare_centers(deltas(grpidx{g},i),deltas(grpidx{h},i),'multcompare',ncomp);
            fprintf('%s Deltas %s vs %s: %s\n',features(i),groups(g),groups(h),stats.resultsstr);
        end
    end
    [p,~,mcstats] = anova1(deltas(:,i),genotypes,'off');
    fprintf('%s Deltas 1-way ANOVA p = %s\n',features(i),get_stars(p));
end