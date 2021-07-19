% Construct filter
figopt=1;
clearvars -except figopt f

animals = {'Bones', 'Odo', 'Sulu', 'Worf', 'Beverly', 'Chakotay', 'Nerys',...
    'SevenOfNine', 'Rain', 'OBrien', 'Picard', 'Riker', 'Kes', 'Neelix', 'Quark',...
    'Sato', 'Garrett', 'Guinan', 'Keeler', 'Dax', 'TPol', 'Doctor', 'Tuvok', 'Bashir',...
    'Scotty'};

epochfilter = [];
epochfilter{1} = {'task','(contains(lower($env), ''veh'') && ~contains($descript,''FAIL'') && $dur==120)'};

datafilter = [];  %only works with single channel specified
datafilter{1} = {'chinfo','(isequal($area,''ca1'') && contains($layer,''pyr 1''))'};

timefilter = [];
timefilter{1} = {'<function> get2dstate <argname> immobilecutoff <argval> 1','($immobilitytime > 30)'};

binrange = [5 50];
f = createfilter('animal', animals, 'epochs', epochfilter, 'data', datafilter, 'excludetime', timefilter);

%Specify function and run
f = setfilterfunction(f, 'calcriprates2', {'ripples'},'appendindex',2,'bins',binrange);
f = runfilter(f);

% save('E:\LFP Data\Hyperactivity\EJ Analysis\E3 E4 CA1 ripple rates.mat', 'f');

%% Plot results

%ONLY WORKS WITH ONE CHANNEL PER MOUSE
% ------------------------------
% Figure and Font Sizes
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
set(0,'DefaultAxesFontName','Arial')
tfont = 20; % title font
xfont = 20;
yfont = 20;

if figopt==1   % PLOT RIPRATES by indiv and group
    %aggregate data over animals
    binnum = length(binrange)-1;
    animnames = [];
    animgeno = [];
    cols = [];
    ratemat = nan(length(f),5);
    for a = 1:length(f)
        if ~isempty(f(a).output)
            results = f(a).output.calcriprates2.results;
            for g = 1:length(results)  %iterate thru conditions
                totaldur = 0;
                durs{a} = [];
                rates{a} = [];
                sizes{a} = [];
                counts{a} = [];
                timecounts{a} = [];
                for e = 1:length(results{g})  %iterate through sessions & epochs, combine counts
                    totaldur = totaldur + results{g}{e}.validdur;
                    durs{a} = [durs{a}; results{g}{e}.validdur];
                    counts{a} = [counts{a}; results{g}{e}.counts(1:binnum,:)];
                    rates{a} = [rates{a} results{g}{e}.rates(1:binnum,:)];  %get rid of extra row (all zeros)
                    ratemat(a,e) = results{g}{e}.rates(1:binnum,:);
                    sizes{a} = [sizes{a}; results{g}{e}.sizes];
                    timecounts{a} = [timecounts{a}; results{g}{e}.timecounts'];
                end
                
            end
        else
            totaldur = 0;
            rates{a} = 0;
            sizes{a} = 0;
        end
        animnames = [animnames; f(a).animal(3) ];
        animgeno = [animgeno; f(a).animal(4)];
    end
    
    for s = 1:binnum
        for a = 1:length(f)
            avgrate(a) = mean(rates{a}(s,:));
        end
        
    end
    avgrate'
    
    %% plot baseline and SD for rip detection
elseif figopt==2
    
    animnames = [];
    animgeno = [];
    cols = [];
    
    for a = 1:length(f)
        if ~isempty(f(a).output)
            results = f(a).output.calcriprates2.results;
            for g = 1:length(results)  %iterate thru conditions
                baseline{a} = [];
                sd{a} = [];
                for e = 1:length(results{g})  %iterate through epochs, combine counts
                    baseline{a} = [baseline{a} results{g}{e}.baseline];
                    sd{a} = [sd{a} results{g}{e}.std];
                end
            end
        else
            baseline{a} = [];
            sd{a} = [];
        end
        animnames = [animnames; f(a).animal(3) ];
        cols = [cols; f(a).animal{5} ];
        animgeno = [animgeno; f(a).animal(4)];
    end
    
    b = [];
    for a = 1:length(f)
        b(a) = mean(baseline{a});
    end
    b'
    for a = 1:length(f)
        b(a) = mean(sd{a});
    end
    b'
    
    for a = 1:length(f)
        subplot(1,2,1)
        hold on
        h = bar(a, mean(baseline{a}));
        errorbar(a,mean(baseline{a}),std(baseline{a}),'.')
        fprintf('%f\t%f\n',mean(baseline{a}),mean(sd{a}));
        col = cols(a);
        set(h,'FaceColor',col)
        title('baseline')
        
        subplot(1,2,2)
        hold on
        h = bar(a, mean(sd{a}));
        errorbar(a,mean(sd{a}),std(sd{a}),'.')
        set(h,'FaceColor',col)
        title('standard deviation')
    end
end