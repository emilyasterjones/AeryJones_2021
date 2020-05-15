function lme = calcjointLMM(vehdata1,cnodata1,vehdata2,cnodata2,logtransform)

if logtransform
    vehmins = cellfun(@min,vehdata1);
    cnomins = cellfun(@min,cnodata1);
    datamin = min([vehmins cnomins]);
    if datamin<=0
        bump = -datamin+0.001;
        vehdata1 = cellfun(@(x) x+bump,vehdata1,'UniformOutput',false);
        cnodata1 = cellfun(@(x) x+bump,cnodata1,'UniformOutput',false);
    end
    vehdata1 = cellfun(@log,vehdata1,'UniformOutput',false);
    cnodata1 = cellfun(@log,cnodata1,'UniformOutput',false);
    
    vehmins = cellfun(@min,vehdata2);
    cnomins = cellfun(@min,cnodata2);
    datamin = min([vehmins cnomins]);
    if datamin<=0
        bump = -datamin+0.001;
        vehdata2 = cellfun(@(x) x+bump,vehdata2,'UniformOutput',false);
        cnodata2 = cellfun(@(x) x+bump,cnodata2,'UniformOutput',false);
    end
    vehdata2 = cellfun(@log,vehdata2,'UniformOutput',false);
    cnodata2 = cellfun(@log,cnodata2,'UniformOutput',false);
end

%initialize
animal = zeros(length(vehdata1{1})+length(cnodata1{1}),1);
len = [vehdata1{1}; cnodata1{1}];
x = NaN(1,4);

for a = 1:length(vehdata1)
    % Collect data for table format
    len = [vehdata1{a}; cnodata1{a}];
    %categorical variables are stored as 0, 1, 2, ...
    animal = repmat(a-1,length(vehdata1{a})+length(cnodata1{a}),1);
    condition = [zeros(length(vehdata1{a}),1); ones(length(cnodata1{a}),1)];
    genotype = zeros(length(vehdata1{a})+length(cnodata1{a}),1);
    x = [x; len animal condition genotype];
end

for a = 1:length(vehdata2)
    % Collect data for table format
    len = [vehdata2{a}; cnodata2{a}];
    %categorical variables are stored as 0, 1, 2, ...
    animal = repmat(a-1+length(vehdata1),length(vehdata2{a})+length(cnodata2{a}),1);
    condition = [zeros(length(vehdata2{a}),1); ones(length(cnodata2{a}),1)];
    genotype = ones(length(vehdata2{a})+length(cnodata2{a}),1);
    x = [x; len animal condition genotype];
end
x = x(2:end,:);

%fit a linear mixed effects model to the data
tbl = table(x(:,1),x(:,2),x(:,3),x(:,4),'VariableNames',{'Measure','Animal','Treatment','Genotype'});

lme = fitlme(tbl,'Measure~Genotype*Treatment+((Genotype*Treatment)|Animal)');
% fprintf('%d\t%f\n',1,lme.ModelCriterion.AIC)
% double(lme.Coefficients(end,6))
% lmefit = fitted(lme);
end