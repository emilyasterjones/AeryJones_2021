
function [filtwts] = designeegfilt(srate,locutoff,hicutoff)

% eegfilt() -  (high|low|band)-iass filter data using two-way least-squares 
%              FIR filtering. Multiple data channels and epochs supported.
%              Requires the MATLAB Signal Processing Toolbox.
% Usage:
%  >> [filtwts] = eegfilt(data,srate,locutoff,hicutoff, ...
%                                             epochframes,filtorder);
% Inputs:
%   srate       = data sampling rate (Hz)
%   locutoff    = low-edge frequency in pass band (Hz)  {0 -> lowpass}
%   hicutoff    = high-edge frequency in pass band (Hz) {0 -> highpass}
%   epochframes = frames per epoch (filter each epoch separately {def/0: data is 1 epoch}
%   filtorder   = length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt     = [0|1] reverse filter (i.e. bandpass filter to notch filter). {0}
%
% Outputs:
%    filtwts    = filter coefficients [smoothdata <- filtfilt(filtwts,1,data)]
%


if ~exist('firls')
   error('*** eegfilt() requires the signal processing toolbox. ***');
end

nyq            = srate*0.5;  % Nyquist frequency
%MINFREQ = 0.1/nyq;
MINFREQ = 0;

minfac         = 3;    % this many (lo)cutoff-freq cycles in filter 
min_filtorder  = 15;   % minimum filter length
trans          = 0.15; % fractional width of transition zones

if locutoff>0 & hicutoff > 0 & locutoff > hicutoff,
    error('locutoff > hicutoff ???\n');
end
if locutoff < 0 | hicutoff < 0,
   error('locutoff | hicutoff < 0 ???\n');
end

if locutoff>nyq,
    error('Low cutoff frequency cannot be > srate/2');
end

if hicutoff>nyq
    error('High cutoff frequency cannot be > srate/2');
end

if nargin<5
   filtorder = 0;
end
if nargin<6
   revfilt = 0;
end

if isempty(filtorder) | filtorder==0,
   if locutoff>0,
     filtorder = minfac*fix(srate/locutoff);
   elseif hicutoff>0,
     filtorder = minfac*fix(srate/hicutoff);
   end
     
   if filtorder < min_filtorder
        filtorder = min_filtorder;
    end
end

if nargin<4
	epochframes = 0;
end

if (1+trans)*hicutoff/nyq > 1
	error('high cutoff frequency too close to Nyquist frequency');
end;

if locutoff > 0 & hicutoff > 0,    % bandpass filter
    f=[MINFREQ (1-trans)*locutoff/nyq locutoff/nyq hicutoff/nyq (1+trans)*hicutoff/nyq 1]; 
    m=[0       0                      1            1            0                      0]; 
elseif locutoff > 0                % highpass filter
 if locutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - highpass cutoff freq must be > %g Hz\n\n',MINFREQ*nyq));
 end
 f=[MINFREQ locutoff*(1-trans)/nyq locutoff/nyq 1]; 
 m=[   0             0                   1      1];
elseif hicutoff > 0                %  lowpass filter
 if hicutoff/nyq < MINFREQ
    error(sprintf('eegfilt() - lowpass cutoff freq must be > %g Hz',MINFREQ*nyq));
 end
 f=[MINFREQ hicutoff/nyq hicutoff*(1+trans)/nyq 1]; 
 m=[     1           1              0                 0];
else
    error('You must provide a non-0 low or high cut-off frequency');
end
if revfilt
    m = ~m;
end;

filtwts = firls(filtorder,f,m); % get FIR filter coefficients

end