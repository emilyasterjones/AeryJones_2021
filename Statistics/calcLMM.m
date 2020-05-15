function [lme, vehcoeff, cnocoeff] = calcLMM(vehdata,cnodata,logtransform)

bump = 0;
if logtransform
    vehmins = cellfun(@min,vehdata);
    cnomins = cellfun(@min,cnodata);
    datamin = min([vehmins cnomins]);
    if datamin<=0
        bump = -datamin+0.001;
        vehdata = cellfun(@(x) x+bump,vehdata,'UniformOutput',false);
        cnodata = cellfun(@(x) x+bump,cnodata,'UniformOutput',false);
    end
    vehdata = cellfun(@log,vehdata,'UniformOutput',false);
    cnodata = cellfun(@log,cnodata,'UniformOutput',false);
end

%initialize
animal = zeros(length(vehdata{1})+length(cnodata{1}),1);
len = [vehdata{1}; cnodata{1}];
x = NaN(1,3);

for a = 1:length(vehdata)
    % Collect data for table format
    len = [vehdata{a}; cnodata{a}];
    %categorical variables are stored as 0, 1, 2, ...
    animal = repmat(a-1,length(vehdata{a})+length(cnodata{a}),1);
    condition = [zeros(length(vehdata{a}),1); ones(length(cnodata{a}),1)];
    x = [x; len animal condition];
end
x = x(2:end,:);

%fit a linear mixed effects model to the data
tbl = table(x(:,1),x(:,2),x(:,3),'VariableNames',{'Measure','Animal','Treatment'});
lme = fitlme(tbl,'Measure~Treatment+(Treatment|Animal)');
lmefit = fitted(lme);
vehcoeff = zeros(length(vehdata),1);
cnocoeff = zeros(length(vehdata),1);
for a = 1:length(vehdata)
    vehcoeff(a) = unique(lmefit(tbl.Animal == a-1 & tbl.Treatment == 0));
    cnocoeff(a) = unique(lmefit(tbl.Animal == a-1 & tbl.Treatment == 1));
end

if logtransform
    vehcoeff = exp(vehcoeff);
    cnocoeff = exp(cnocoeff);
    if datamin<=0
        vehcoeff = vehcoeff - bump;
        cnocoeff = cnocoeff - bump;
    end 
end


% %plot the residuals & manually inspect them; should be Gaussian & uncorrelated
% figure
% subplot(1,3,1)
% plotResiduals(lme)
% title('Res. distr.')
% subplot(1,3,2)
% qqplot(residuals(lme))
% title('Res. QQplot')
% subplot(1,3,3)
% scatter(tbl.Animal,residuals(lme))
% title('Res. related to animal')

% disp(lme)
end