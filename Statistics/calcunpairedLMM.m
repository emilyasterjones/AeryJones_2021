function lme = calcunpairedLMM(vehdata1,vehdata2,logtransform)

if logtransform
    datamin = min(cellfun(@min,vehdata1));
    if datamin<=0
        bump = -datamin+0.001;
        vehdata1 = cellfun(@(x) x+bump,vehdata1,'UniformOutput',false);
    end
    vehdata1 = cellfun(@log,vehdata1,'UniformOutput',false);
    
    datamin = min(cellfun(@min,vehdata2));
    if datamin<=0
        bump = -datamin+0.001;
        vehdata2 = cellfun(@(x) x+bump,vehdata2,'UniformOutput',false);
    end
    vehdata2 = cellfun(@log,vehdata2,'UniformOutput',false);
end

%initialize
animal = zeros(length(vehdata1{1}),1);
len = [vehdata1{1}];
x = NaN(1,3);

for a = 1:length(vehdata1)
    % Collect data for table format
    len = vehdata1{a};
    %categorical variables are stored as 0, 1, 2, ...
    animal = repmat(a-1,length(vehdata1{a}),1);
    sex = zeros(length(vehdata1{a}),1);
    x = [x; len animal sex];
end

for a = 1:length(vehdata2)
    % Collect data for table format
    len = vehdata2{a};
    %categorical variables are stored as 0, 1, 2, ...
    animal = repmat(a-1+length(vehdata1),length(vehdata2{a}),1);
    sex = ones(length(vehdata2{a}),1);
    x = [x; len animal sex];
end
x = x(2:end,:);

%fit a linear mixed effects model to the data
tbl = table(x(:,1),x(:,2),x(:,3),'VariableNames',{'Measure','Animal','Sex'});

lme = fitlme(tbl,'Measure~Sex+(Sex|Animal)');
% fprintf('%d\t%f\n',1,lme.ModelCriterion.AIC)
% double(lme.Coefficients(end,6))
% lmefit = fitted(lme);
end