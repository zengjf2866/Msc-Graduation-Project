%%
% test example of a 6x6 system (6 inputs, 6 outputs)
% open-loop stable
clear;
clc;

A = [-0.01 2.2 3.3 4.4 5.5 6.6;
      0 -0.02 -2.5 3.1 -5.8 3.6;
      0 0 -0.03 7.1 -1.8 2.6;
      0 0 0 -0.04 6 1.2;
      0 0 0 0 -0.05 -3;  
      0 0 0 0 0 -0.06];
  
B = rand(6,6);
C = rand(6,6);
D = zeros(6,6);

sys = ss(A,B,C,D);
% bodemag(sys)

nu = 6;
ny = 6;

%%
% identification
% excitation signal generation
% MIMO experiment with only odd frequency excited
% generate odd frequency first, identify first 10000 Hz
% ExcitedHarm is the excited harmonics required for MIMO identification

DefFreq.fs      =   20000; % sampling frequency of the generator
DefFreq.fres    =   0.5;% frequency spacing in Hz between the (odd) harmonics 
DefFreq.fmin    =   0.5;% lowest excited frequency in Hz (for lowpass design of odd multisines fmin = fres/2)
DefFreq.fmax    =   2000;% largest excited frequency in Hz
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
RmsValues = [4; 4; 4; 4; 4; 4];

TheSignal = {};
for i = 1:1:total_periods
    TheSignal{i} = Calc_MIMO_Multisine(ExcitedHarm, N, AmplitudeExcitedHarm, RmsValues);
end

% these are the orthogonal exciation signals d used in 6 independent experiments
multisine_exp_1 = zeros(nu,nh*total_periods);
% multisine_exp_2 = zeros(nu,nh*total_periods);
% multisine_exp_3 = zeros(6,nh*total_periods);
% multisine_exp_4 = zeros(6,nh*total_periods);
% multisine_exp_5 = zeros(6,nh*total_periods);
% multisine_exp_6 = zeros(6,nh*total_periods);
for i = 0:1:(total_periods-1)
    for j = 1:1:N
        multisine_exp_1(:,(i*N+j)) = TheSignal{i+1}(:,1,j); 
%         multisine_exp_2(:,(i*N+j)) = TheSignal{i+1}(:,2,j);
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
StopTime = length(u0)/fs - Ts;

sys_dt = c2d(sys,Ts);
Ad = sys_dt.A; Bd = sys_dt.B; Cd = sys_dt.C; Dd = sys_dt.D;
for i = 1:1
   if i == 1
       EXCSIG = u0;
   end    
   data{i} = sim('system66.slx');
end

%%
% Extract the data
y_1 = data{1,1}.y.signals.values'; 
d_1 = data{1,1}.d.signals.values'; 

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
method.startfreq = 0.5;              % defines the start frequency of the analysis 
method.stopfreq = 1999,5;                % defines the stop frequency of the analysis

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
		semilogx(freq, db(squeeze(G(jj,ii,:))), 'r', freq, db(squeeze(sys_resp(jj,ii,:))), 'k--');
	end % ii
end % jj
subplot(ny,nu,1);
title('identified G: red; G_true: black;');
zoom on; shg











