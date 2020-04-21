% %% Construct filter
% clearvars -except animals
animals = {'Adobo_LT', 'Caraway_LT', 'Chives_LT', ...
    'Cinnamon_LT', 'Coriander_LT', 'Fenugreek_LT', 'GaramMasala_LT', 'Salt_LT', ...
    'Baharat_LT', 'Cardamom_LT', 'Jerk_LT', 'Mace_LT', 'Mustard_LT', ...
    'Tarragon_LT', 'Vanilla_LT',...
    'Basil_LT', 'Cumin_LT', 'Dill_LT', 'Nutmeg_LT', 'Paprika_LT', 'Parsley_LT', 'Sumac_LT',... 
    'Pepper_LT', 'Sage_LT', 'Anise_LT', 'Thyme_LT', 'OldBay_LT', 'Rosemary_LT',...
    'Provence_LT', 'Saffron_LT'};

epochfilter{1} = {'task','(isequal($treatment,''Veh'') && contains($env,''LT'') && ~contains($descript,''FAIL''))'};

datafilter = [];  %doesnt matter what channel you pick
datafilter{1} = {'chinfo','(isequal($area,''ca1'') && contains($layer,''pyr''))'};

timefilter = [];  %doesnt matter since going to ignore excludetimes
timefilter{1} =  {'pos', '($vel < 50)'};

binrange = [0 1 3 5 50];
% exponents = -2:5;
% binrange = 2.^exponents;
% binrange = 0:0.5:20;

f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);

%Specify function and run
f = setfilterfunction(f, 'velquant', {'pos'},'appendindex',1,'bins',binrange);
f = runfilter(f);

velresults = f;

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

f = createfilter('animal', animals, 'epochs', epochfilter, 'excludetime', timefilter);

%% aggregate data over animals
binnum = length(binrange)-1;
animnames = [];

for a = 1:length(f)
    results = velresults(a).output.velquant.results;
    totaltime{a} = 0;
    avgspeed{a} = 0;
    maxspeed{a} = 0;
    rawcounts{a} = zeros(1,binnum);
    reltemp{a} = zeros(1,binnum);
    sleepcounts{a} = [];
    sleepday{a} = [];
    totaltimeday{a} = [];
    sleep{a} = 0;
    for g = 1:length(results)  %iterate thru conditions
        for e = 1:length(results{g})  %iterate through epochs, combine counts
            totaltimeday{a}{e} = results{g}{e}.totaltime;
            totaltime{a} = totaltime{a} + results{g}{e}.totaltime;
            avgspeed{a} = avgspeed{a} + results{g}{e}.avgspeed;
            maxspeed{a} = maxspeed{a} + results{g}{e}.maxspeed;
            rawcounts{a} = rawcounts{a} + results{g}{e}.rawcounts';
            reltemp{a} = [reltemp{a}; results{g}{e}.relcounts'];
            
            sleepcounts{a} = f(a).excludetime{g}{e}(:,2)-f(a).excludetime{g}{e}(:,1);
            sleepday{a}{e} = sum(sleepcounts{a});
            sleep{a} = sleep{a} + sum(sleepcounts{a});
        end
        avgspeed{a} = avgspeed{a}/length(results{g});
        maxspeed{a} = maxspeed{a}/length(results{g});
    end
    avgspeed{a} = avgspeed{a}/length(results);
    maxspeed{a} = maxspeed{a}/length(results);
    animnames = [animnames; f(a).animal(3) ];
    rawcounts{a} = rawcounts{a}/totaltime{a};
    reltemp{a} = mean(reltemp{a}(2:end,:),1);
end

%% plot speed per animal
figure
hold on
for a = 1:length(f)
    subplot(5,6,a)
    bar(rawcounts{a},'BarWidth',1)
end

figure
hold on
for a = 1:length(f)
    bar(a, sum(rawcounts{a}(2:end))*totaltime{a})
end
title('Average time/session >1cm/s')
set(gca, 'XTickLabel', '')
ypos = -max(ylim)/50;
text(1:a,repmat(ypos,a,1), animnames','horizontalalignment','right','Rotation',90,'FontSize',15)
ylabel('mins')

figure
hold on
for a = 1:length(f)
    bar(a, avgspeed{a})
end
title('Average velocity')
set(gca, 'XTickLabel', '')
ypos = -max(ylim)/50;
text(1:a,repmat(ypos,a,1), animnames','horizontalalignment','right','Rotation',90,'FontSize',15)
ylabel('cm/s')

%% plot total time
figure
hold on
for a = 1:length(f)
    h = bar(a, totaltime{a}/60);
    col = char(f(a).animal(5));
    set(h,'FaceColor',col)
end
title('Total time recorded')
set(gca, 'XTickLabel', '')
ypos = -max(ylim)/50;
text(1:a,repmat(ypos,a,1), animnames','horizontalalignment','right','Rotation',90,'FontSize',15)
ylabel('mins')

%% plot percentage
figure
for s = 1:binnum
    subplot(5,1,s)
    hold on
    for a = 1:length(f)
        h = bar(a, 100*rawcounts{a}(s)/totaltime{a});
        col = bitget(find('krgybmcw'==char(f(a).animal(5)))-1,1:3); %convert char color to RGB
        set(h,'FaceColor',col)
    end
    title(sprintf('Percent time spent at %d-%d cm/s',binrange(s),binrange(s+1)))
    set(gca, 'XTickLabel', '')
    ypos = -max(ylim)/50;
    text(1:a,repmat(ypos,a,1), animnames','horizontalalignment','right','Rotation',90,'FontSize',15)

end

%% plot time totals
figure
for s = 1:binnum
    subplot(5,1,s)
    hold on
    for a = 1:length(f)
        h = bar(a, rawcounts{a}(s)/60);
        col = bitget(find('krgybmcw'==char(f(a).animal(5)))-1,1:3); %convert char color to RGB
        set(h,'FaceColor',col)
    end
    title(sprintf('Total time spent at %d-%d cm/s',binrange(s),binrange(s+1)))
    set(gca, 'XTickLabel', '')
    ypos = -max(ylim)/50;
    text(1:a,repmat(ypos,a,1), animnames','horizontalalignment','right','Rotation',90,'FontSize',15)
    ylabel('mins')
end

%% plot time asleep
figure
hold on
for a = 1:length(f)
    h = bar(a, (totaltime{a}-sleep{a})/60);
    %disp((totaltime{a}-sleep{a})/60); %%%%%%%%%%%%%
    col = char(f(a).animal(5));
    set(h,'FaceColor',col)
end
title('Total time asleep (0-1cm/s over >30s)')
set(gca, 'XTickLabel', '')
ypos = -max(ylim)/50;
text(1:a,repmat(ypos,a,1), animnames','horizontalalignment','right','Rotation',90,'FontSize',15)

%% plot time asleep each session
figure
hold on
for a = 1:length(f)
    results = velresults(a).output.velquant.results;
    for e = 1:length(results{g})
        h = bar(e+(a-1)*length(results{g}), (totaltimeday{a}{e}-sleepday{a}{e})/60); %a+length(f)*(e-1)
        disp((totaltimeday{a}{e}-sleepday{a}{e})/60);
        col = char(f(a).animal(5));
        set(h,'FaceColor',col)
    end
end
title('Total time asleep each day')
set(gca, 'XTickLabel', '')
ypos = -max(ylim)/50;
text([1:e:a*e],repmat(ypos,a,1), animnames','horizontalalignment','left','Rotation',0,'FontSize',15)