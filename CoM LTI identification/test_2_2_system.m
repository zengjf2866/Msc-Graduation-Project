%%
% test a 2x2 system identification
% continous-time system
clear;
clc;

A = [-0.003 0.039 0 -0.322; -0.065 -0.319 7.74 0; 0.020 -0.101 -0.429 0; 0 0 1 0];
B = [0.010 1; -0.18 -0.04; -1.16 0.598; 0 0];
C = [1 0 0 0; 0 -1 0 7.74];
D = zeros(2,2);

sys = ss(A,B,C,D);
bodemag(sys)

nu = size(B,2);
ny = size(C,1);

%%
% identification
% excitation signal generation
% MIMO experiment with only odd frequency excited
% generate odd frequency first, identify first 10000 Hz
% ExcitedHarm is the excited harmonics required for MIMO identification

DefFreq.fs      =   5000; % sampling frequency of the generator
DefFreq.fres    =   0.1;% frequency spacing in Hz between the (odd) harmonics 
DefFreq.fmin    =   0.01;% lowest excited frequency in Hz (for lowpass design of odd multisines fmin = fres/2)
DefFreq.fmax    =   500;% largest excited frequency in Hz
% DefFreq.frat    =   % ratio between consecutive (odd) harmonics for a logarithmic frequency spacing

Nblock = 5;                 % one out of three consecutive harmonics is randomly eliminated
Spacing = 'lin';            % linear frequency spacing
MultiType = 'odd';          % no even excited harmonics

total_periods = 15; % defination of total periods

[ExcitedHarm, N_odd_lin, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);

fs = NewDefFreq.fs;

% --------------------------------------------------------------------------------------------------
% generate orthogonal excitation

N = 2*DefFreq.fs;                                  % number of time domain samples , each period lasts for 2 seconds and have 40000 data points in total
nh = length(ExcitedHarm);

% definition amplitude spectra inputs
AmplitudeExcitedHarm = ones(nu, nh);
% AmplitudeExcitedHarm(4:6, :) = 0.225; % relative amplitudes

% defition rms values inputs
RmsValues = [4; 4];

TheSignal = {};
for i = 1:1:total_periods
    TheSignal{i} = Calc_MIMO_Multisine(ExcitedHarm, N, AmplitudeExcitedHarm, RmsValues);
end

% these are the orthogonal exciation signals d used in 6 independent experiments
multisine_exp_1 = zeros(nu,nh*total_periods);
multisine_exp_2 = zeros(nu,nh*total_periods);
% multisine_exp_3 = zeros(6,nh*total_periods);
% multisine_exp_4 = zeros(6,nh*total_periods);
% multisine_exp_5 = zeros(6,nh*total_periods);
% multisine_exp_6 = zeros(6,nh*total_periods);
for i = 0:1:(total_periods-1)
    for j = 1:1:N
        multisine_exp_1(:,(i*N+j)) = TheSignal{i+1}(:,1,j); 
        multisine_exp_2(:,(i*N+j)) = TheSignal{i+1}(:,2,j);
%         multisine_exp_3(:,(i*N+j)) = TheSignal{i+1}(:,3,j);
%         multisine_exp_4(:,(i*N+j)) = TheSignal{i+1}(:,4,j);
%         multisine_exp_5(:,(i*N+j)) = TheSignal{i+1}(:,5,j);
%         multisine_exp_6(:,(i*N+j)) = TheSignal{i+1}(:,6,j);
    end
end

sampling_frequency = fs;
L = length(multisine_exp_1(1,:));
total_time = (0:L-1)*(N_odd_lin*total_periods/sampling_frequency)/L; % total time is 30 seconds
% f = 0.5:1:1999.5;

%%
% Interact with simulink model, 6 independent experiments and obtain d,u
% and e for closed-loop identification
clear data;
Ts = 1/fs;
u0 = multisine_exp_1;
u0_2 = multisine_exp_2;
StopTime = length(u0)/fs - Ts;

sys_dt = c2d(sys,Ts);
Ad = sys_dt.A; Bd = sys_dt.B; Cd = sys_dt.C; Dd = sys_dt.D;
for i = 1:2
   if i == 1
       EXCSIG = u0;
   else
       EXCSIG = u0_2;
   end    
   data{i} = sim('boeing747.slx');
end

%%
% Extract the data
y_1 = data{1,1}.y.signals.values'; 
d_1 = data{1,1}.d.signals.values'; 

y_2 = data{1,2}.y.signals.values'; 
d_2 = data{1,2}.d.signals.values'; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling the ArbLocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data 
data_1.u = d_1;                            % row index is the input number; the column index the time instant (in samples) 
data_1.y = y_1;                             % row index is the output number; the column index the time instant (in samples) 
data_1.Ts = 1/fs;                         % sampling period

% method
method.dof = 10;                        % degrees of freedom of the variance estimate
method.order = 2;                       % order local polynomial approximation
method.startfreq = 0.1;              % defines the start frequency of the analysis 
method.stopfreq = 500;                % defines the stop frequency of the analysis

% local polynomial estimate FRF and its variance
[CZ, Z, freq, G, CvecG, dof, CL] = ArbLocalPolyAnal(data_1, method);

% estimated output noise power spectrum
S = CZ.n(1:ny,1:ny,:); 

% estimated variances FRM entries: keep the diagonal elements CvecG only
F = length(freq);
varG = zeros(ny, nu, F);
for kk=1:F
    varG(:, :, kk) = reshape(diag(CvecG(:, :, kk)), [ny, nu]);
end % kk

%%
% comparison estimated and true FRM
figure(1)
% G_2 = frd(G_2,freq,FrequencyUnit = 'Hz');
sys_resp = freqresp(sys,freq,'Hz');

mm = 0;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		plot(freq, db(squeeze(G(jj,ii,:))), 'r', freq, db(squeeze(sys_resp(jj,ii,:))), 'k--');
	end % ii
end % jj
subplot(ny,nu,1);
title('identified G: red; G_true: black;');
zoom on; shg











