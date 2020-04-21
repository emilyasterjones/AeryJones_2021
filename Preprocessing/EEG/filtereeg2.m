function filteredeegstruct = filtereeg2(raweegstruct, filterstruct, varargin)
% use as reference for how to filter
%   filtereeg2(raweegstruct, filterstruct, options)
%
% Applies the specified zero-phase finite impulse response FILTER to the
% the specified RAWEEG structure. Returns a FILTEREDEEG structure whose
% data field contains the filtered amplitude, instantaneous phase
% (extracted by Hilbert transform), and envelope magnitude (extracted by
% Hilbert transform).
%
% Note that the sampling rate of the FILTER must match the sampling rate
% of the eeg data. 
% filterstruct is one of two types of structure:
%    for FIR filters, it has the following fields:
%	'kernel' 	The FIR filter kernel
%	'samprate'	The sampling rate used to create the kernel
%	'descript'	A description of the kernel
%
%    for IIR filters, the structure is the output of the sptool gui with fields:
%	'tf'		The numerator and denominator of the IIR filter kernel
%	'Fs'		The sampling rate used to create the kernel
%	'descript'	A description of the filter
%
% Finally, be aware that the very beginning and end time segments of the
% filtered output are not to be trusted because the required samples are
% missing from the input.
%
% differs from filtereeg in how deals with NaNs.  Filtereeg cuts out NaN,
% performs filter, then replaces NaN.  This code applies filter to data
% with NaNs, may enlarge NaN regions.

% written by asinger, from filtereeg

% options:
%	'int16', 1	Specifies the data should be in int16 format.  Phases
%			are multiplied by 10000 in this case.
%			Default 0 (double precision)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ERROR CHECKING


% Check that the supplied raweeg data is a valid column vector
if isempty(raweegstruct.data)
    error('eeg data is empty')
elseif length(raweegstruct.data) ~= numel(raweegstruct.data)
    error('eeg data must be a column vector')
else
    [numrows, numcols] = size(raweegstruct.data);
    if numrows==1
        error('eeg data must be a column vector')
    end
end


% Figure out what sort of filter we have
firfilter = 0;
if (isfield(filterstruct, 'kernel'))
    firfilter = 1;
    if isempty(filterstruct.kernel)
	error('filter kernel is empty')
    elseif length(filterstruct.kernel) ~= numel(filterstruct.kernel)
	error('filter kernel must be a vector')
    end
    samprate = filterstruct.samprate;
elseif (~isfield(filterstruct, 'tf'))
    error('the filter is not in the right format');
else
    samprate = filterstruct.Fs;
end




% Error checking on the samprates
samprate_ratio = raweegstruct.samprate/samprate;

% Check that sampling rates of the filter and the data match to within some
% very small epsilon; warn the user if this is not the case
if abs(samprate_ratio - 1) > 1e-3
    warning('filterstruct.samprate does not match eegstruct.samprate');
end

% Check that filterstruct.samprate does not exceed eegstruct.samprate
if (samprate_ratio - 1) < -1e-3
    disp(samprate_ratio);
    error('filterstruct.samprate must not exceed eegstruct.samprate');
end

% Check that samprate_ratio is close to an integer
if abs(round(samprate_ratio) - samprate_ratio) > 1e-3
    error('filterstruct.samprate must be very close to an integer fraction of eegstruct.samprate');
end

% assign the options:
useint16 = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'int16'
            useint16 = varargin{option+1};
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOW FOR THE ACTUAL BUSINESS

%filtdata = filtfilt(filterstruct.kernel, 1 , raweegstruct);
%clear raweegstruct
%hdata = hilbert(filtdata);
%env = abs(hdata);
%phase = angle(hdata);

% edit the descript field of this cell to reflect that this is 
% *filtered* data, not the original raw eeg trace
if round(samprate_ratio)==1
    filteredeegstruct.descript = ...
    strvcat(raweegstruct.descript,'filtered through: ', ...
    [blanks(size(filterstruct.descript,1))' filterstruct.descript]);
else
    filteredeegstruct.descript = ...
    strvcat(raweegstruct.descript,'decimated, filtered through: ', ...
    [blanks(size(filterstruct.descript,1))' filterstruct.descript], ...
    'and interpolated back to original samprate');
end

%filtereegstruct
filteredeegstruct.samprate = raweegstruct.samprate;
filteredeegstruct.starttime = raweegstruct.starttime;
if (firfilter) 
    filtdata = filtfilt(filterstruct.kernel, 1 , double(raweegstruct.data));
else
    filtdata = filtfilt(filterstruct.tf.num, filterstruct.tf.den, ...
    			double(raweegstruct.data));
end
clear raweegstruct
% apply Hilbert transform
hdata = hilbert(filtdata);
env = abs(hdata);
phase = angle(hdata);
clear hdata
if (~useint16)
    filteredeegstruct.data = intNaN(length(filtdata),3);
    filteredeegstruct.data(:,1) = filtdata;
    clear filtdata
    filteredeegstruct.data(:,2) = phase;
    clear phase
    filteredeegstruct.data(:,3) = env;
    clear env
    filteredeegstruct.fields = ...
    'filtered_amplitude1 instantaneous_phase2 envelope_magnitude3';
else
    filteredeegstruct.data = zeros(length(filtdata),3, 'int16');
    filteredeegstruct.data(:,1) = int16(filtdata);
    clear filtdata
    filteredeegstruct.data(:,2) = int16(phase*10000);
    clear phase
    filteredeegstruct.data(:,3) = int16(env);
    clear env
    filteredeegstruct.fields = ...
    'filtered_amplitude instantaneous_phase*10000 envelope_magnitude';
end

% add in a copy of the filter kernel
if (firfilter)
    filteredeegstruct.filterkernel = filterstruct.kernel;
else
    filteredeegstruct.filterkernel = filterstruct.tf;
end
filteredeegstruct.filtersamprate = samprate;

