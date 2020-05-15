function [p, F, df] = calcInteraction(vehdata1,cnodata1,vehdata2,cnodata2,logtransform)
lme = calcLMM([vehdata1 vehdata2], [cnodata1 cnodata2], logtransform);
jointlme = calcjointLMM(vehdata1, cnodata1, vehdata2, cnodata2, logtransform);
table = compare(lme, jointlme);
df = table.deltaDF(2);
F = table.LRStat(2);
p = table.pValue(2);
end