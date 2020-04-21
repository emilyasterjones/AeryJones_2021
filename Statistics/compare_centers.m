function stats = compare_centers(x1, x2, varargin)

%Compares 2 populations using appropriate tests for measures of centrality
%Performs Shapiro-Wilk test for normality & F test for equal variance
%Uses appropriate statistic to compare populations
%Sidak correction for multiple tests if needed
%Input: 2 vectors to be compared
%Varargin: whether the comparison is paired, whether it should be corrected
%for multiple comparisons
%Output: stats stucture containing name of test used, p value, test statistic,
%degrees of freedom (if applicable), and effect size

paired = 0;
multcompare = 1;

for option = 1:2:length(varargin)-1
    if ischar(varargin{option})
        switch(varargin{option})
            case 'paired'
                paired = varargin{option+1};
            case 'multcompare'
                multcompare = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

[~, p1] = swtest(x1); %Shapiro-Wilk test for normality
[~, p2] = swtest(x2);

if (p1<0.01 || p2<0.01)
	stats.center = 'median';
	stats.x1m = nanmedian(x1);
	stats.x2m = nanmedian(x2);
	stats.teststat = 'Z';
	if (paired)
		stats.test = 'Wilcoxon matched pairs signed rank test';
		[p, ~, sr] = signrank(x1,x2,'method','approximate');
	else
		stats.test = 'Wilcoxon rank sum test';
		[p, ~, sr] = ranksum(x1,x2,'method','approximate');
	end
	stats.testval = sr.zval;
    stats.df = NaN;
else	
	stats.center = 'mean';
	stats.x1m = nanmean(x1);
	stats.x2m = nanmean(x2);
	stats.teststat = 't';
	[~,p,~,~] = vartest2(x1, x2); %F test for equal variance
	if (p<0.01)
		if (paired)
            stats.test = 'Paired t test';
			[~,p,~,tdf] = ttest(x1, x2);
		else
			stats.test = 'Unpaired Welch''s t test';
			[~,p,~,tdf] = ttest2(x1, x2, 'Vartype', 'unequal');
		end
	else
		if (paired)
			stats.test = 'Paired t test';
			[~,p,~,tdf] = ttest(x1, x2);
		else
			stats.test = 'Unpaired t test';
			[~,p,~,tdf] = ttest2(x1, x2);
		end
	end
	stats.df = tdf.df;
	stats.testval =tdf.tstat;
end

stats.p = 1-(1-p)^multcompare; %Sidak correction for multiple comparisons

stats.resultsstr = sprintf('%s, %s(%d) = %.2f, p = %.4f, %s: x1 = %.3f, x2 = %.3f',...
    stats.test, stats.teststat, stats.df, stats.testval, stats.p, stats.center,...
    stats.x1m, stats.x2m);