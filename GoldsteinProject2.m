function GoldsteinProject2

% This function is a script for event-related potential 
% analysis using the SPM12 EEG utilities. Notably, it does not make use of
% SPM12 jobman function, but directly calls the respective spm_eeg_xxx
% functions as was required in previous versions of SPM. 
%
%          Input - None
%
%         Output - 4 Figures comparing different analysis parameters
%                  Figure 1 - compares different rereferencing paramters
%                  Figure 2 - compares different high-pass filter cut-offs
%                  Figure 3 - compares different low-pass filter cut-offs
%                  Figure 4 - compares peak2peak artefact rejection to none
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------
clc
close all

% data source directory
sdir   = 'C:\Users\Ursula\Documents\MATLAB\Project2\EEG_Example_Data_and_Script';

% subject stubs
SJ      = '04_LIHOV';

% create input structure for spm_eeg_erp
S       = [];
S.sdir  = [sdir filesep];
S.SJ    = SJ;

% perform analysis
% Reference Graph
                Average(S)

%  Figure 1 - compare re-referencing
       TP10Rereferencing(S)
          NoRereferencing(S)
        
% Figure 2 - compare high contrasts
                    HC11(S)
                    HC22(S)

% Figure 3 - compare low contrasts
                    LC11(S)
                    LC22(S)
       
% Figure 4 - compare artefact correction to none
                Artefact(S)

% visualize results for 2 occipital/parietal electrodes 
% -------------------------------------------------------------------------
      Ave = fullfile(sdir, SJ,  ['meflfhr_' SJ '.mat'])      ; % concatenate filename
  TP10Ref = fullfile(sdir, SJ,  ['2m2e2fl2fh2r_' SJ '.mat']) ; % concatenate filename
  NoReref = fullfile(sdir, SJ,  ['3m3e3fl3fh_' SJ '.mat'])   ; % concatenate filename
      HC1 = fullfile(sdir, SJ,  ['4m4e4fl4fhr_' SJ '.mat'])  ; % concatenate filename
      HC2 = fullfile(sdir, SJ,  ['5m5e5fl5fhr_' SJ '.mat'])  ; % concatenate filename
      LC1 = fullfile(sdir, SJ,  ['6m6e6flfhr_' SJ '.mat'])   ; % concatenate filename
      LC2 = fullfile(sdir, SJ,  ['7m7e7flfhr_' SJ '.mat'])   ; % concatenate filename
Artefacts = fullfile(sdir, SJ,  ['maeflfhr_' SJ '.mat'])     ; % concatenate filename

R       = spm_eeg_load(Ave)                                  ; % load the MEEG object/data
N       = spm_eeg_load(NoReref)                              ; % load the MEEG object/data
T       = spm_eeg_load(TP10Ref)                              ; % load the MEEG object/data
H       = spm_eeg_load(HC1)                                  ; % load the MEEG object/data
I       = spm_eeg_load(HC2)                                  ; % load the MEEG object/data
L       = spm_eeg_load(LC1)                                  ; % load the MEEG object/data
M       = spm_eeg_load(LC2)                                  ; % load the MEEG object/data
A       = spm_eeg_load(Artefacts)                            ; % load the MEEG object/data


r       = time(R).*1000                                      ; % peristimulus time
n       = time(N).*1000                                      ; % peristimulus time
t       = time(T).*1000                                      ; % peristimulus time
b       = time(H).*1000                                      ; % peristimulus time
c       = time(I).*1000                                      ; % peristimulus time
l       = time(L).*1000                                      ; % peristimulus time
m       = time(M).*1000                                      ; % peristimulus time
a       = time(A).*1000                                      ; % peristimulus time

eoi     = {'PO8', 'O2'}                                      ; % electrodes of interest
 
 eoi_idx = indchannel(R, eoi)                           ; % indices of electrodes of interest
eoi_idx2 = indchannel(T, eoi)                           ; % indices of electrodes of interest
eoi_idx3 = indchannel(N, eoi)                           ; % indices of electrodes of interest
eoi_idx4 = indchannel(H, eoi)                           ; % indices of electrodes of interest
eoi_idx5 = indchannel(I, eoi)                           ; % indices of electrodes of interest
eoi_idx6 = indchannel(L, eoi)                           ; % indices of electrodes of interest
eoi_idx7 = indchannel(M, eoi)                           ; % indices of electrodes of interest
eoi_idx8 = indchannel(A, eoi)                           ; % indices of electrodes of interest


h = figure                                              ; % new figure
set(h, 'Color', [1 1 1])                                ; % figure formatting

% plot data for selected channels for rereferencing conditions
for i = 1:length(eoi_idx)
    subplot(1,2,i)     
    hold on
    plot(n, N(eoi_idx2(i),:,1)', 'b'  , 'LineWidth', 2)
    plot(n, N(eoi_idx2(i),:,2)',  '--b' ,'LineWidth', 2) 
    plot(r, R(eoi_idx(i),:,1)', 'r'  , 'LineWidth', 2)
    plot(r, R(eoi_idx(i),:,2)','--r', 'LineWidth', 2)
    plot(t, T(eoi_idx3(i),:,1)', 'k'  , 'LineWidth', 2)
    plot(t, T(eoi_idx3(i),:,2)',  '--k' ,'LineWidth', 2)
    title([ R.chanlabels{eoi_idx(i)}], 'FontWeight', 'Normal')
    xlabel('ms')
    ylabel('\muV', 'Rotation', 0)
    set(gca, 'FontSize', 10, 'LineWidth', 2, 'FontName', 'Times New Roman')  
    xlim([r(1) r(end)])
    ylim([-25 25])
    if i == 1
        legend('Recording HC', 'Recording LC', 'Average HC', 'Average LC', 'TP10 HC', 'TP10 LC')
    end

end 

h = figure                                              ; % new figure
set(h, 'Color', [1 1 1])                                ; % figure formatting

% plot data for selected channels for high-pass filter cut-offs conditions
for i = 1:length(eoi_idx)
    subplot(1,2,i)     
    hold on
    plot(b, H(eoi_idx4(i),:,1)', 'b'  , 'LineWidth', 2)
    plot(b, H(eoi_idx4(i),:,2)',  '--b' ,'LineWidth', 2) 
    plot(r, R(eoi_idx(i),:,1)', 'r'  , 'LineWidth', 2)
    plot(r, R(eoi_idx(i),:,2)','--r', 'LineWidth', 2)
    plot(c, I(eoi_idx5(i),:,1)', 'k'  , 'LineWidth', 2)
    plot(c, I(eoi_idx5(i),:,2)',  '--k' ,'LineWidth', 2)
    title([ R.chanlabels{eoi_idx(i)}], 'FontWeight', 'Normal')
    xlabel('ms')
    ylabel('\muV', 'Rotation', 0)
    set(gca, 'FontSize', 10, 'LineWidth', 2, 'FontName', 'Times New Roman')  
    xlim([r(1) r(end)])
    ylim([-25 25])
    if i == 1
        legend('HPF 1 HC', 'HPF 1 LC', 'Average HC', 'Average LC', 'HPF 2 HC', 'HPF 2 LC')
    end

end 

h = figure                                              ; % new figure
set(h, 'Color', [1 1 1])                                ; % figure formatting

% plot data for selected channels for low-pass filter cut-offs conditions
for i = 1:length(eoi_idx)
    subplot(1,2,i)     
    hold on
    plot(l, L(eoi_idx6(i),:,1)', 'b'  , 'LineWidth', 2)
    plot(l, L(eoi_idx6(i),:,2)',  '--b' ,'LineWidth', 2) 
    plot(r, R(eoi_idx(i),:,1)', 'r'  , 'LineWidth', 2)
    plot(r, R(eoi_idx(i),:,2)','--r', 'LineWidth', 2)
    plot(m, M(eoi_idx7(i),:,1)', 'k'  , 'LineWidth', 2)
    plot(m, M(eoi_idx7(i),:,2)',  '--k' ,'LineWidth', 2)
    title([ R.chanlabels{eoi_idx(i)}], 'FontWeight', 'Normal')
    xlabel('ms')
    ylabel('\muV', 'Rotation', 0)
    set(gca, 'FontSize', 10, 'LineWidth', 2, 'FontName', 'Times New Roman')  
    xlim([r(1) r(end)])
    ylim([-25 25])
    if i == 1
        legend('LPF 20 HC', 'LPF 20 LC', 'Average HC', 'Average LC', 'LPF 10 HC', 'LPF 10 LC')
    end

end 

h = figure                                              ; % new figure
set(h, 'Color', [1 1 1])                                ; % figure formatting

% plot data for selected channels for artefact rejection conditions
for i = 1:length(eoi_idx)
    subplot(1,2,i)     
    hold on
    plot(a, A(eoi_idx8(i),:,1)', 'b'  , 'LineWidth', 2)
    plot(a, A(eoi_idx8(i),:,2)',  '--b' ,'LineWidth', 2) 
    plot(r, R(eoi_idx(i),:,1)', 'r'  , 'LineWidth', 2)
    plot(r, R(eoi_idx(i),:,2)','--r', 'LineWidth', 2)
    title([ R.chanlabels{eoi_idx(i)}], 'FontWeight', 'Normal')
    xlabel('ms')
    ylabel('\muV', 'Rotation', 0)
    set(gca, 'FontSize', 10, 'LineWidth', 2, 'FontName', 'Times New Roman')  
    xlim([r(1) r(end)])
    ylim([-25 25])
    if i == 1
        legend('Peak2Peak HC', 'Peak2Peak LC', 'Average HC', 'Average LC')
    end
end 

end

function Average(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% Data conversion
% -------------------------------------------------------------------------
S.dataset               = [sdir SJ '.eeg' ]             ; % source file
S.outfile               = [sdir SJ]                     ; % target file
S.mode                  = 'continuous'                  ; % convert data as continuous 
S.checkboundary         = 1                             ; % 1 = check if there are breaks in the file and do not read across those breaks [default], 0 = ignore breaks (not recommended).
S.saveorigheader        = 0                             ; % 1 = save original data header with the dataset, 0 = do not keep the original header [default]

% convert data set
spm_eeg_convert(S)            

% Rereference to average reference
% -------------------------------------------------------------------------
% prepare montaging
S                       = []                            ; % spm_eeg_montage structure initilization
S.D                     = [sdir SJ '.mat']              ; % source filename
D                       = spm_eeg_load(S.D)             ; % load MEEG object to obtain channel labels
S.montage.labelorg      = D.chanlabels                  ; % N x 1 cell-array: original labels
S.montage.labelnew      = D.chanlabels                  ; % N x 1 cell-array: new labels

% create re-reference matrix 
N                       = numel(S.montage.labelorg)     ; % number of channels
ref_tra                 = NaN(N,N)                      ; % reference matrix initialization

% create new montage matrix for  average reference
for i = 1:N
    for j = 1:N        
        if i == j 
            ref_tra(i,j) = (N-1)/N;
        else
            ref_tra(i,j) = -1/N;
        end
    end
end

S.montage.tra           = ref_tra                       ; % M x N matrix
S.keepothers            = 'no'                          ; % 'yes'/'no': keep or discard the channels not involved in the montage [default: 'yes']
S.prefix                = 'r_'                          ; % prefix for the output file (default - 'M')

% create new data montage
spm_eeg_montage(S);

% High-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir 'r_' SJ '.mat']         ; % source filename
S.band                  = 'high'                        ; % filterband [low|high|bandpass|stop]
S.freq                  = 0.1                           ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = 'fh'                          ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Low-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir 'fhr_' SJ '.mat']       ; % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 40                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = 'fl'                          ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir 'flfhr_' SJ '.mat']         ; % source filename
S.timewin               = [-100 500]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = 'e'                               ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
    S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
    S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
    S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

% perform epoching
spm_eeg_epochs(S); 

% Artefact correction
% -------------------------------------------------------------------------
% none

% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir 'eflfhr_' SJ '.mat']    ; % source filename
S.robust                = 1                             ; % use standard averaging
S.prefix                = 'm'                           ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

function TP10Rereferencing(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set with rereferencing to the electrode TP10.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% Data conversion
% -------------------------------------------------------------------------
S.dataset               = [sdir SJ '.eeg' ]             ; % source file
S.outfile               = [sdir SJ]                     ; % target file
S.mode                  = 'continuous'                  ; % convert data as continuous 
S.checkboundary         = 1                             ; % 1 = check if there are breaks in the file and do not read across those breaks [default], 0 = ignore breaks (not recommended).
S.saveorigheader        = 0                             ; % 1 = save original data header with the dataset, 0 = do not keep the original header [default]

% convert data set
spm_eeg_convert(S)            

% Rereference to average reference
% -------------------------------------------------------------------------
% prepare montaging
S                       = []                            ; % spm_eeg_montage structure initilization
S.D                     = [sdir SJ '.mat']              ; % source filename
D                       = spm_eeg_load(S.D)             ; % load MEEG object to obtain channel labels
S.montage.labelorg      = D.chanlabels                  ; % N x 1 cell-array: original labels
S.montage.labelnew      = D.chanlabels                  ; % N x 1 cell-array: new labels

% create re-reference matrix 
N                       = numel(S.montage.labelorg)     ; % number of channels
ref_tra                 = NaN(N,N)                      ; % reference matrix initialization

% create new montage matrix for TP10 reference
for i = 1:N
    for j = 1:N        
        if i == j 
            if i == 22
                ref_tra(i,j) = 0;
            else
                ref_tra(i,j) = 1;
            end
        else
            if j == 22
                ref_tra(i,j) = -1;
            else
                ref_tra(i,j) = 0;
            end
        end
    end
end

S.montage.tra           = ref_tra                       ; % M x N matrix
S.keepothers            = 'no'                          ; % 'yes'/'no': keep or discard the channels not involved in the montage [default: 'yes']
S.prefix                = '2r_'                         ; % prefix for the output file (default - 'M')

% create new data montage
spm_eeg_montage(S);

% High-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                             ; % spm_eeg_filter structure initilization
S.D                     = [sdir '2r_' SJ '.mat']         ; % source filename
S.band                  = 'high'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 0.1                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                  ; % optional filter type, IIR (default) or FIR
S.order                 = 5                              ; % filter order (default - 5 for Butterworth)
S.prefix                = '2fh'                          ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Low-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir '2fh2r_' SJ '.mat']     ; % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 40                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '2fl'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '2fl2fh2r_' SJ '.mat']      ; % source filename
S.timewin               = [-100 500]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = '2e'                              ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
    S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
    S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
    S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

% perform epoching
spm_eeg_epochs(S); 

% Artefact correction
% -------------------------------------------------------------------------
% None                                       


% Averaging
% -------------------------------------------------------------------------
S                       = []                                 ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '2e2fl2fh2r_' SJ '.mat']   ; % source filename
S.robust                = 1                                  ; % use standard averaging
S.prefix                = '2m'                               ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

function NoRereferencing(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set with no rereferecing.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% Data conversion
% -------------------------------------------------------------------------
S.dataset               = [sdir SJ '.eeg' ]             ; % source file
S.outfile               = [sdir SJ]                     ; % target file
S.mode                  = 'continuous'                  ; % convert data as continuous 
S.checkboundary         = 1                             ; % 1 = check if there are breaks in the file and do not read across those breaks [default], 0 = ignore breaks (not recommended).
S.saveorigheader        = 0                             ; % 1 = save original data header with the dataset, 0 = do not keep the original header [default]

% convert data set
spm_eeg_convert(S)            

% Rereference to average reference
% -------------------------------------------------------------------------
% None

% High-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir  SJ '.mat']             ; % source filename
S.band                  = 'high'                        ; % filterband [low|high|bandpass|stop]
S.freq                  = 0.1                           ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '3fh_'                        ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Low-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir '3fh_' SJ '.mat']       ; % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 40                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '3fl'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '3fl3fh_' SJ '.mat']        ; % source filename
S.timewin               = [-100 500]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = '3e'                              ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
    S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
    S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
    S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

% perform epoching
spm_eeg_epochs(S); 

% Artefact correction
% -------------------------------------------------------------------------
% none

% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '3e3fl3fh_' SJ '.mat']  ; % source filename
S.robust                = 1                             ; % use standard averaging
S.prefix                = '3m'                          ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

function HC11(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set with a high-pass filter cut-off of 1 Hz.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% High-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir 'r_' SJ '.mat']         ; % source filename
S.band                  = 'high'                        ; % filterband [low|high|bandpass|stop]
S.freq                  = 1                             ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '4fh'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Low-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir '4fhr_' SJ '.mat']      ; % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 40                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '4fl'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '4fl4fhr_' SJ '.mat']       ; % source filename
S.timewin               = [-100 500]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = '4e'                              ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
    S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
    S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
    S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

% perform epoching
spm_eeg_epochs(S); 

% Artefact correction
% -------------------------------------------------------------------------
% None

% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '4e4fl4fhr_' SJ '.mat'] ; % source filename
S.robust                = 1                             ; % use standard averaging
S.prefix                = '4m'                          ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

function HC22(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set with a high-pass filter cut-off of 2 Hz.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% High-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir 'r_' SJ '.mat']         ; % source filename
S.band                  = 'high'                        ; % filterband [low|high|bandpass|stop]
S.freq                  = 2                             ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '5fh'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Low-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir '5fhr_' SJ '.mat']      ; % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 40                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '5fl'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '5fl5fhr_' SJ '.mat']       ; % source filename
S.timewin               = [-100 500]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = '5e'                              ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
    S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
    S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
    S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

% perform epoching
spm_eeg_epochs(S); 

% Artefact correction
% -------------------------------------------------------------------------
% None

% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '5e5fl5fhr_' SJ '.mat'] ; % source filename
S.robust                = 1                             ; % use standard averaging
S.prefix                = '5m'                          ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

function LC11(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set with a low-pass filter cut-off of 20 Hz.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% Low-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir 'fhr_' SJ '.mat']       ; % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 20                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '6fl'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '6flfhr_' SJ '.mat']        ; % source filename
S.timewin               = [-100 500]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = '6e'                              ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
    S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
    S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
    S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

% perform epoching
spm_eeg_epochs(S); 

% Artefact correction
% -------------------------------------------------------------------------
% None

% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '6e6flfhr_' SJ '.mat']  ; % source filename
S.robust                = 1                             ; % use standard averaging
S.prefix                = '6m'                          ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

function LC22(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set with a low-pass filter cut-off of 10 Hz.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% Low-Pass Filtering
% -------------------------------------------------------------------------
S                       = []                            ; % spm_eeg_filter structure initilization
S.D                     = [sdir 'fhr_' SJ '.mat']       ; % source filename
S.band                  = 'low'                         ; % filterband [low|high|bandpass|stop]
S.freq                  = 10                            ; % cutoff frequency [Hz]
S.type                  = 'butterworth'                 ; % optional filter type, IIR (default) or FIR
S.order                 = 5                             ; % filter order (default - 5 for Butterworth)
S.prefix                = '7fl'                         ; % prefix for the output file (default - 'f')

% filter the data
spm_eeg_filter(S); 

% Epoching
% --------------------------------------------------------------------------
S                       = []                                ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '7flfhr_' SJ '.mat']        ; % source filename
S.timewin               = [-100 500]                        ; % time window in PST ms
S.bc                    = 1                                 ; % baseline-correct the data (1 - yes, 0 - no).
S.prefix                = '7e'                              ; % prefix for the output file (default - 'e')
S.trialdef              = []                                ; % initialize structure array for trial definition

% the current markers of interest are S 1 and S 2
for i = 1:2
    S.trialdef(1,i).conditionlabel  = ['S  ', num2str(i)]    ; % string label for the condition in the output file
    S.trialdef(1,i).eventtype       = 'Stimulus'             ; % string label of event type
    S.trialdef(1,i).eventvalue      = ['S  ', num2str(i)]    ; % string, numeric, or empty value of the condition
end

% perform epoching
spm_eeg_epochs(S); 

% Artefact correction
% -------------------------------------------------------------------------
% None

% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir '7e7flfhr_' SJ '.mat']  ; % source filename
S.robust                = 1                             ; % use standard averaging
S.prefix                = '7m'                          ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

function Artefact(S)

% This function performs a basic event-related potential analysis for a
% single participant's data set with a peak 2 peak artefact correction.
%
%   Inputs
%       S   : data location structure with fields
%               .spm_dir    : string, source directory
%               .SJ         : string, subject stub
%   Outputs
%           None, saves data to disc
%
% Copyright (C) Dirk Ostwald & Ursula Goldstein
% -------------------------------------------------------------------------

% Initialization
% ------------------------------------------------------------------------
spm('defaults', 'eeg')                                  ; % load the spm eeg defaults
sdir                    = S.sdir                        ; % specify source directory
SJ                      = S.SJ                          ; % specify subject stub

% Artefact correction
% -------------------------------------------------------------------------

S = []                                                      ; % spm_eeg_artefact structure initilization 
S.D = [sdir 'eflfhr_' SJ '.mat']                            ; % source filename
S.mode = 'reject'                                           ; % mode to reject artefacts
S.badchanthresh = 0.2                                       ; % when is a channel bad
S.methods.channels = 'all'                                  ; % specify number of channels
S.methods.fun = 'peak2peak'                                 ; % specify method used
S.methods.settings.threshold = 50                           ; % specify threshhold for method
S.append = true                                             ; %
S.prefix = 'a'                                              ; % prefix for the output file (default - 'a') 

% perform artefact correction and rejection
spm_eeg_artefact(S)

% Averaging
% -------------------------------------------------------------------------
S                       = []                            ; % initialize spm_eeg_epochs input structure
S.D                     = [sdir 'aeflfhr_' SJ '.mat']   ; % source filename
S.robust                = 1                             ; % use standard averaging
S.prefix                = 'm'                           ; % prefix for the output file (default - 'm')

% perform robust averaging
spm_eeg_average(S); 

end

