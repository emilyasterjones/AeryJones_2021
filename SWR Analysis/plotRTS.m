figopt = 2;
clear title

%% plot a single channel

% ------------------------------
% Figure and Font Sizes
set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',16);
tfont = 20; % title font
xfont = 20;
yfont = 20;
% ---------------------------------------
    % plot ALL RTS for a mouse in correct layout
if figopt ==2
    % ------------------------------
    % Figure and Font Sizes
    set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
    set(0,'defaultaxesfontsize',10);
    tfont = 10; % title font
    xfont = 10;
    yfont = 10;
    clr = {'b','r','g','c','m','y','k','r'};
    % ---------------------------------------
    
    
    for a = 1:length(f)  %should always be 1
        
        results = f(a).output.RTspecgram3.results;
        %         vehresults = fveh(a).output.RTspecgram3.results;
        for g = 1:length(results)  %iterate thru conditions
            for t = 1:size(results{1}{1}.trigindex,1)  %iterate thru trigchans
                %combine data over epochs for trig and probe
                
                numevents = 0;
                
                for p = 1:size(results{1}{1}.probeindex,1) %iterate through probechans
                    %iterate thru probe chans, combine eps, plot
                    
                    probe{t}{p} = [];
                    for e = 1:length(results{g}) %epochs
                        if length(results{g}{e}.riptimes)>1
                            %                         valids = find(results{1}{e}.stds{1} >= 5); %ensure that rips are over 5SD threshold
                            %                         probe{t}{p} = cat(3,probe{t}{p},results{g}{e}.Sprobe{t}{p}(:,1:38,valids));  %38 full spectrum, 1:11 SG
                            probe{t}{p} = cat(3,probe{t}{p},results{g}{e}.Sprobe{t}{p}(:,:,:)); %%%1:11
                            Sprobe = results{g}{e}.Sprobe{t}{p};
                            meanPprobe = results{g}{e}.meanPprobe{t}{p};
                            stdPprobe = results{g}{e}.stdPprobe{t}{p};
                            % VmeanPprobe = vehresults{g}{e}.meanPprobe{t}{p};
                            % VstdPprobe = vehresults{g}{e}.stdPprobe{t}{p};
                            
                            % for i = 1:size(Sprobe,1)
                            %     for j = 1:size(Sprobe,3)
                            %         Sprobe(i,:,j) = Sprobe(i,:,j).*stdPprobe + meanPprobe;
                            %         Sprobe(i,:,j) = (Sprobe(i,:,j) - VmeanPprobe)./VstdPprobe;
                            %     end
                            % end
                            probe{t}{p} = cat(3,probe{t}{p},Sprobe);
                            
                            times = results{g}{e}.times;
                            freqs = results{g}{e}.freqs(:); %38 full spectrum, 1:11 SG  %%%1:11
                        end
                    end
                    probemean{t}{p} = median(probe{t}{p},3);%mean
                end
                numevents = size(probe{t}{1},3);
                
                
                %plot all RTS
                figure;
                hold on;
                
                shank = size(results{1}{1}.probeindex,1);
                
                if ~isempty(probe{t}{1})
                    if shank <20
                        %two shank
                        set(gcf,'Position',[1 1 3095 364]);
                        disp('twosh')
                        for s = 1:shank
                            chanpos = results{1}{1}.probeindex(s,3);
                            %                         disp(sprintf('s %d chanpos %d', s, chanpos))
                            if chanpos<3
                                subplot(2,16,chanpos)
                            else
                                subplot(2,16,chanpos+1)
                            end
                            imagesc(times, freqs,probemean{t}{s}')
                            set(gca, 'Clim', [0 3]);
                            set(gca,'YDir','normal');
                            details = sprintf('%d', chanpos);
                            title(details, 'FontSize',10,'Fontweight','normal');
                            %ylabel('Freq','FontSize',10,'Fontweight','normal');
                            %xlabel('Time(s)','FontSize',10,'Fontweight','normal');
                            set(gca,'XLim',[min(times) max(times)]);
                            set(gca,'YLim',[min(freqs) max(freqs)]);
                            % Plot Line at 0 ms - Start of ripple
                            hold on
                            ypts = freqs;
                            xpts = 0*ones(size(ypts));
                            plot(xpts , ypts, 'k--','Linewidth',2);
                            %colorbar;
                        end
                        suptitle = sprintf('%s %drips trigggered off ch%d',f(a).animal{3},numevents,results{g}{t}.trigindex(t,3));
                        [ax,h]=suplabel(suptitle ,'t');
                        tightfig;
                        
                    else
                        %four shank
                        set(gcf,'Position',[1 1 1595 964]);
                        disp('four sh')
                        for s = 1:shank
                            chanpos = results{1}{1}.probeindex(s,3);
                            
                            subplot(4,8,chanpos)
                            imagesc(times, freqs,probemean{t}{s}');
                            colormap hot %jet
                            set(gca, 'Clim', [0 4]);
                            set(gca,'YDir','normal');
                            details = sprintf('%d', chanpos);
                            title(details, 'FontSize',10,'Fontweight','normal');
                            %ylabel('Freq','FontSize',10,'Fontweight','normal');
                            %xlabel('Time(s)','FontSize',10,'Fontweight','normal');
                            set(gca,'XLim',[min(times) max(times)]);
                            set(gca,'YLim',[min(freqs) max(freqs)]);
                            % Plot Line at 0 ms - Start of ripple
                            hold on
                            ypts = freqs;
                            xpts = 0*ones(size(ypts));
                            plot(xpts , ypts, 'w--','Linewidth',2);
                            colorbar;
                            
                        end
                        suptitle = sprintf('%s %drips trigggered off ch%d, 5SD',f(a).animal{3},numevents,results{g}{t}.trigindex(t,3));
                        [ax,h]=suplabel(suptitle ,'t');
                        tightfig;
                    end
                end
            end
        end
    end

elseif figopt==3  % plot all chans JUST low freqs
    
    % ------------------------------
    % Figure and Font Sizes
    set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
    set(0,'defaultaxesfontsize',10);
    tfont = 10; % title font
    xfont = 10;
    yfont = 10;
    clr = {'b','r','g','c','m','y','k','r'};
    % ---------------------------------------
    
    
    for a = 1:length(f)  %should always be 1
        
        results = f(a).output.RTspecgram3.results;
        
        infofile = sprintf('%s%schinfo.mat',f(a).animal{2},f(a).animal{3});
        load(infofile)
        
        for g = 1:length(results)  %iterate thru conditions
            for t = 1:size(results{1}{1}.trigindex,1)  %iterate thru trigchans
                %combine data over epochs for trig and probe
                
                numevents = 0;
                times = results{1}{1}.times(56:135);
                freqs = results{1}{1}.freqs(1:11);
                
                for p = 1:size(results{1}{1}.probeindex,1) %iterate through probechans
                    %iterate thru probe chans, combine eps, plot
                    
                    probe{t}{p} = [];
                    for e = 1:length(results{g}) %epochs
                        %combine events from all eps
                        valids = find(results{1}{e}.stds{1} >= 5); %ensure that rips are over 5SD threshold
                        probe{t}{p} = cat(3,probe{t}{p},results{g}{e}.Sprobe{t}{p}(56:135,1:11,valids));
                    end
                    numevents = size(probe{t}{p},3);
                    probemean{t}{p} = mean(probe{t}{p},3);
                end
                
                
                %plot all RTS
                figure;
                hold on;
                %four shank
                set(gcf,'Position',[1 1 1595 964]);
                disp('four sh')
                for s = 1:shank
                    chanpos = results{1}{1}.probeindex(s,3);
                    
                    subplot(4,8,chanpos)
                    
                    imagesc(times, freqs,probemean{t}{s}');
                    set(gca, 'Clim', [0 3]);
                    set(gca,'YDir','normal');
                    details = sprintf('%d %s %s', chanpos,chinfo{1}{1}{chanpos}.area,chinfo{1}{1}{chanpos}.layer);
                    title(details, 'FontSize',10,'Fontweight','normal');
                    %ylabel('Freq','FontSize',10,'Fontweight','normal');
                    %xlabel('Time(s)','FontSize',10,'Fontweight','normal');
                    set(gca,'XLim',[min(times) max(times)]);
                    set(gca,'YLim',[min(freqs) max(freqs)]);
                    % Plot Line at 0 ms - Start of ripple
                    hold on
                    ypts = freqs;
                    xpts = 0*ones(size(ypts));
                    plot(xpts , ypts, 'k--','Linewidth',2);
                    if s==shank
                        colorbar;
                    end
                    
                end
                suptitle = sprintf('%s %drips trigggered off ch%d, 3SD',f(a).animal{3},numevents,results{g}{t}.trigindex(t,3));
                [ax,h]=suplabel(suptitle ,'t');
                tightfig;
            end
        end
    end
    
elseif figopt==4  %plot a single designated channel
    
    % Figure and Font Sizes
    set(0,'defaultaxesfontweight','normal'); set(0,'defaultaxeslinewidth',2);
    set(0,'defaultaxesfontsize',18);
    tfont = 16; % title font
    xfont = 16;
    yfont = 16;
    sesscount = 0;
    animcount = 0;
    probemean = 0;
    rips = NaN(length(f),8);
    for a = 1:length(f)
        results = f(a).output.RTspecgram4.results;
        %         vehresults = f1(a).output.RTspecgram3.results;
        for g = 1:length(results)  %iterate thru conditions
            for t = 1:size(results{1}{1}.trigindex,1)  %iterate thru trigchans
                %combine data over epochs for probe
                
                for p = 1:size(results{1}{1}.probeindex,1) %iterate through probechans
                    %iterate thru probe chans, combine eps, plot
                    
                    probe = [];
                    for e = 1:length(results{g}) %epochs
                        rips(a,e) = length(results{g}{e}.riptimes);
                        if length(results{g}{e}.freqs)>1
                            %combine events from all eps
                            peak = max(max(results{g}{e}.Sprobe{t}{p}));
                            results{g}{e}.Sprobe{t}{p} = results{g}{e}.Sprobe{t}{p};
                            %                             results{g}{e}.Sprobe{t}{p} = (results{g}{e}.Sprobe{t}{p}-vehresults{g}{e}.meanPprobe{t}{p})./vehresults{g}{e}.stdPprobe{t}{p};
                            probe = cat(3,probe,results{g}{e}.Sprobe{t}{p});
                            times = results{g}{e}.times;
                            freqs = results{g}{e}.freqs;
                        end
                    end
                    probemean = probemean + nanmean(probe,3);
                    animcount = animcount+1;
                    
                end
            end
        end
    end
    probemean = probemean./animcount;
    figure;
    hold on;
    
    imagesc(times, freqs, probemean'); %(16:56,1:10)
    set(gca, 'Clim', [0 3]);
    colorbar;
    colormap hot
    set(gca,'YDir','normal');
    details = sprintf('%s ch%d', f(a).animal{3},results{1}{1}.probeindex(p,3));
    title(details, 'FontSize',20,'Fontweight','normal');
    ylabel('Freq','FontSize',20,'Fontweight','normal');
    xlabel('Time(s)','FontSize',20,'Fontweight','normal');
    set(gca,'XLim',[min(times) max(times)]);
    set(gca,'YLim',[min(freqs) max(freqs)]);
    % Plot Line at 0 ms - Start of ripple
    hold on
    ypts = freqs;
    xpts = 0*ones(size(ypts));
    plot(xpts , ypts, 'k--','Linewidth',2);
    
end