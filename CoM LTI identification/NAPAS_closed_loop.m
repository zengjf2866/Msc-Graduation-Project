%%
clear;
clc;

CFB = load('PD_Control_NAPAS.mat');
CFB = CFB.C;
CFB = ss(CFB);

%%
% rigid-body controller construction
% load rigid-body controller
% CFB = load('system_identification_MIMO - test\CFB.mat');
% CFB = load('Control.mat'); % use low bandwidth controller for better signal-to-noise ratio
% CFB = CFB.C;

% CFB = load('CFB.mat'); % use low bandwidth controller for better signal-to-noise ratio
% CFB = CFB.CFB;
Ac = CFB.A; Bc = CFB.B; Cc = CFB.C; Dc = CFB.D;
% controller discretization
Controller_ct = ss(Ac,Bc,Cc,Dc);

% MMPA construction
% A, B and C matrices of MMPA are fixed at specific time instant and position
Ap = load('experiment_A.mat'); Ap = Ap.A;
Bp = load('experiment_B.mat'); Bp = Bp.B;
Cp = load('experiment_C.mat'); Cp = Cp.C;
% plant discretization
Plant_ct = ss(Ap,Bp,Cp,zeros(6,6));

%%
loops_1 = loopsens(Plant_ct,Controller_ct); 

%%
% bodemag(loops_1.So,'r',{1,2e3*2*pi})
% legend('Sensitivity')
% grid on

%%
bodemag(CFB,'r',{1,2e3*2*pi})
legend('Two controllers')
grid on

%%
% excitation signal generation
% MIMO experiment with only odd frequency excited
% generate odd frequency first
% ExcitedHarm is the excited harmonics required for MIMO identification

DefFreq.fs      =   40000; % sampling frequency of the generator
DefFreq.fres    =   0.1;% frequency spacing in Hz between the (odd) harmonics 
DefFreq.fmin    =   0.5;% lowest excited frequency in Hz (for lowpass design of odd multisines fmin = fres/2)
DefFreq.fmax    =   4000;% largest excited frequency in Hz
% DefFreq.frat    =   % ratio between consecutive (odd) harmonics for a logarithmic frequency spacing

Nblock = Inf;                 % one out of three consecutive harmonics is randomly eliminated
Spacing = 'lin';            % linear frequency spacing
MultiType = 'full';          % no even excited harmonics

total_periods = 15; % defination of total periods

[ExcitedHarm, N_odd_lin, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);

fs = NewDefFreq.fs;

% --------------------------------------------------------------------------------------------------
% generate orthogonal excitation

nu = 6;                                     % number of inputs
N = fs*2;                                  % number of time domain samples , each period lasts for 2 seconds and have 40000 data points in total
nh = length(ExcitedHarm);

% definition amplitude spectra inputs
AmplitudeExcitedHarm = 100000*ones(nu, nh);
% AmplitudeExcitedHarm(4:6, :) = 0.225; % relative amplitudes

% defition rms values inputs
RmsValues = [1; 1; 1; 0.1; 0.1; 0.1];

TheSignal = {};
for i = 1:1:total_periods
    TheSignal{i} = Calc_MIMO_Multisine(ExcitedHarm, N, AmplitudeExcitedHarm, RmsValues);
end

% these are the orthogonal exciation signals d used in 6 independent experiments
multisine_exp_1 = zeros(6,nh*total_periods);
% multisine_exp_2 = zeros(6,nh*total_periods);
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
Ts = 1/fs;

% r = randn(nu, N*total_periods);% white noise reference

%%
% Interact with simulink model, 6 independent experiments and obtain d,u
% and e for closed-loop identification
StopTime = length(multisine_exp_1)/fs - Ts;
for i = 1:1
   if i == 1
       EXCSIG = multisine_exp_1;
%        REFSIG = r;
%    elseif i == 2
%        EXCSIG = multisine_exp_2;
%    elseif i == 3
%        EXCSIG = multisine_exp_3;
%    elseif i == 4
%        EXCSIG = multisine_exp_4;
%    elseif i == 5
%        EXCSIG = multisine_exp_5;    
%    else
%        EXCSIG = multisine_exp_6;
   end
    
   data{i} = sim('system_identification_MIMO - test/simplified_model.slx');
end

%%
exclude = 0;
data{1,1}.u.signals.values = data{1,1}.u.signals.values(exclude*2*sampling_frequency+1:end,:); 
data{1,1}.e.signals.values = data{1,1}.e.signals.values(exclude*2*sampling_frequency+1:end,:); 
data{1,1}.d.signals.values = data{1,1}.d.signals.values(exclude*2*sampling_frequency+1:end,:); 

u_1 = data{1,1}.u.signals.values; u_1 = u_1';
e_1 = data{1,1}.e.signals.values; e_1 = e_1';
d_1 = data{1,1}.d.signals.values; d_1 = d_1';

%%
% construct the data
% for sensitivity
% L = length(d_1);
sensitivity.u = d_1;
sensitivity.y = u_1;
% sensitivity.r = r;
% sensitivity.N = [L L L L L L];
sensitivity.Ts = 1/fs;

% method
method.dof = 10;                        % degrees of freedom of the variance estimate
method.order = 2;                       % order local polynomial approximation
method.startfreq = DefFreq.fmin;              % defines the start frequency of the analysis 
method.stopfreq = DefFreq.fmax;                % defines the stop frequency of the analysis

[CZ_S, Z_S, freq, Sensitivity, CvecG_S, dof_S, CL_S] = ArbLocalPolyAnal(sensitivity, method);

%%
real_sensitivity = frd(loops_1.So,freq,FrequencyUnit='Hz');

mm = 0;
ny = 6;
nu = 6;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		semilogx(freq, db(squeeze(real_sensitivity.ResponseData(jj,ii,:))), 'r',freq, db(squeeze(Sensitivity(jj,ii,:))), 'b');
	end % ii
end % jj
subplot(ny,nu,3);
title('state-space sensitivity: red; identified sensitivity: blue;');
zoom on; shg

%%
% for process sensitivity
process_sensitivity.u = d_1;
process_sensitivity.y = e_1;
% sensitivity.r = zeros(6,600000);
% process_sensitivity.N = [L L L L L L];
process_sensitivity.Ts = 1/fs;

% method
method.dof = 10;                        % degrees of freedom of the variance estimate
method.order = 2;                       % order local polynomial approximation
method.startfreq = DefFreq.fmin;              % defines the start frequency of the analysis 
method.stopfreq = DefFreq.fmax;                % defines the stop frequency of the analysis

[CZ_PS, Z_PS, freq, Process_Sensitivity, CvecG_PS, dof_PS, CL_PS] = ArbLocalPolyAnal(process_sensitivity, method);
Process_Sensitivity(:,:,:) = -Process_Sensitivity(:,:,:);

%%
real_process_sensitivity = frd(loops_1.So*Plant_ct,freq,FrequencyUnit='Hz');
% bodemag(real_process_sensitivity,'r',real_process_sensitivity_2,'b',{1,2e3*2*pi})
% legend('Process Sensitivity')
% grid on

%%
mm = 0;
ny = 6;
nu = 6;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		semilogx(freq, db(squeeze(real_process_sensitivity.ResponseData(jj,ii,:))), 'r',freq, db(squeeze(Process_Sensitivity(jj,ii,:))), 'b');
	end % ii
end % jj
subplot(ny,nu,3);
title('state-space process sensitivity: red; identified sensitivity: blue;');
zoom on; shg

%%
% frd model for H
S = frd(Sensitivity,freq,FrequencyUnit='Hz');
PS = frd(Process_Sensitivity,freq,FrequencyUnit='Hz');
H = PS*inv(S);

H_real = freqresp(Plant_ct,freq,'Hz'); % frequency response of the real state-space plant

%%
mm = 0;
ny = 6;
nu = 6;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		semilogx(freq, db(squeeze(H_real(jj,ii,:))), 'r',freq, db(squeeze(H.ResponseData(jj,ii,:))), 'b');
	end % ii
end % jj
subplot(ny,nu,3);
title('state-space model: red; identified H: blue;');
zoom on; shg






