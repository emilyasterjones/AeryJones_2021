group = 'EV'
% 
clearvars -except group
load(['E:\LMM Data\',group,'.mat'])
cno = who('cno*');
veh = who('veh*');

fprintf('Metric\tVeh Mean\tCNO Mean\tEffect size\tdf\tF\tp\tF log\tp\n')
for i = 1:length(cno)   
    eval(['lme = calcLMM(',veh{i},',',cno{i},',0);'])
    fixedEf = fixedEffects(lme);
    [p,F,~,df2] = coefTest(lme);
    eval(['vehmean = mean(cellfun(@nanmean,',veh{i},'));'])
    eval(['cnomean = mean(cellfun(@nanmean,',cno{i},'));'])
    fprintf('%s\t%g\t%g\t%g\t%g\t%g\t%g\t%s\n',veh{i},vehmean,cnomean,fixedEf(2),df2,F,p,get_stars(p))
        
end