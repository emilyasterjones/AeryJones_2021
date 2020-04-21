function out = comod(index, excludetimes, eeg, varargin)
% how phase of theta modulates power of gamma
% based on dfakk_MI from kenny

%Defaults
srate = 1000;

% Define vector of frequencies whose amplitudes are modulated
% version 1: original
PhaseFreqVector=0:2:50;
AmpFreqVector=10:5:300;
PhaseFreq_BandWidth=4;
AmpFreq_BandWidth=10;

nbin = 18;

%Set options
for option = 1:2:length(varargin)-1
    if isstr(varargin{option})
        switch(varargin{option})
            case 'nbin'
                nbin = varargin{option+1};
            case 'PhaseFreqVector'
                PhaseFreqVector = varargin{option+1};
            case 'AmpFreqVector'
                AmpFreqVector = varargin{option+1};
            case 'PhaseFreq_BandWidth'
                PhaseFreq_BandWidth = varargin{option+1};
            case 'AmpFreq_BandWidth'
                AmpFreq_BandWidth = varargin{option+1};
            otherwise
                error(['Option ',varargin{option},' unknown.']);
        end
    else
        error('Options must be strings, followed by the variable');
    end
end

center = (-pi+pi/(2*nbin)):pi/nbin:(pi-pi/(2*nbin));

for c = 1:length(index.chinfo)
    
    e = double(eeg{index.epochs(1)}{index.epochs(2)}{index.chinfo(c)}.data);
    timevec = geteegtimes(eeg{index.epochs(1)}{index.epochs(2)}{index.chinfo(c)});
    data_length = length(e);
    
    %% Filter and Hilbert transform
    Comodulogram = zeros(length(PhaseFreqVector),length(AmpFreqVector));
    Comodbyphase = zeros(length(PhaseFreqVector),length(AmpFreqVector),2*nbin-1);
    AmpFreqTransformed = zeros(length(AmpFreqVector), data_length);
    PhaseFreqTransformed = zeros(length(PhaseFreqVector), data_length);
    
    for ii=1:length(AmpFreqVector)
        Af1 = AmpFreqVector(ii);
        Af2 = Af1+AmpFreq_BandWidth;
        coeff = designeegfilt(srate,Af1,Af2); % just filtering
        d = 1;
        AmpFreq = int16(filtfilt(coeff,d,e));
        AmpFreqTransformed(ii, :) = abs(hilbert(AmpFreq)); % getting the amplitude envelope
    end
    
    for jj=1:length(PhaseFreqVector)
        Pf1 = PhaseFreqVector(jj);
        Pf2 = Pf1 + PhaseFreq_BandWidth;
        coeff = designeegfilt(srate,Pf1,Pf2); % just filtering
        d = 1;
        PhaseFreq = int16(filtfilt(coeff,d,e));
        PhaseFreqTransformed(jj, :) = angle(hilbert(PhaseFreq)); % this is getting the phase time series
    end
    
    
    %% Apply excludetimes to filtered eeg % pat kenny on back for this block of code :D
    includedindices = ~isExcluded(timevec,excludetimes);
    totalsecofeeg = sum(includedindices)/srate;
    AmpFreqTransformed = AmpFreqTransformed(:,includedindices);
    PhaseFreqTransformed = PhaseFreqTransformed(:,includedindices);
    
    
    %% Do comodulation calculation
    for ii=1:length(PhaseFreqVector)
        comodindex = zeros(1,length(AmpFreqVector));
        comodphase = zeros(length(AmpFreqVector),2*nbin-1);
        for jj=1:length(AmpFreqVector)
            %Compute mean of amplitudes at each phase bin
            temp = accumarray(lookup(PhaseFreqTransformed(ii, :),center(2:end)),...
                AmpFreqTransformed(jj, :), [2*nbin-1 1], @(x) mean(x));
            comodphase(jj,:) = temp/sum(temp);
            comodindex(jj) = (log(2*nbin-1) + sum(comodphase(jj,:).*log(comodphase(jj,:))))/log(2*nbin-1);
        end
        Comodbyphase(ii,:,:) = comodphase;
        Comodulogram(ii,:) = comodindex;
    end
    
    %% Repeat with 1s windows
    timebin = 1:1000:length(AmpFreqTransformed);
    
    for b = 1:length(timebin)-1
        AmpFreqTransformed_bin = AmpFreqTransformed(:,timebin(b):timebin(b+1));
        PhaseFreqTransformed_bin = PhaseFreqTransformed(:,timebin(b):timebin(b+1));
        
        for ii=1:length(PhaseFreqVector)
            comodindex = zeros(1,length(AmpFreqVector));
            comodphase = zeros(length(AmpFreqVector),2*nbin-1);
            for jj=1:length(AmpFreqVector)
                temp = accumarray(lookup(PhaseFreqTransformed_bin(ii, :),center(2:end)),...
                    AmpFreqTransformed_bin(jj, :), [2*nbin-1 1], @(x) mean(x));
                comodphase(jj,:) = temp/sum(temp);
                comodindex(jj) = (log(2*nbin-1) + sum(comodphase(jj,:).*log(comodphase(jj,:))))/log(2*nbin-1);
            end
            Comodulogram_bin(b,ii,:) = comodindex;
        end
    end

    out.comodulogram{c} = Comodulogram;
    out.comodbyphase{c} = Comodbyphase;
    out.totalsecofeeg{c} = totalsecofeeg;
    out.comodulogram_bin{c} = Comodulogram_bin;
    
end
out.frequency_amplitude = AmpFreqVector;
out.frequency_phase = PhaseFreqVector;
out.phasefreq_bandwidth = PhaseFreq_BandWidth;
out.ampfreq_bandwidth = AmpFreq_BandWidth;
out.index = index;
end