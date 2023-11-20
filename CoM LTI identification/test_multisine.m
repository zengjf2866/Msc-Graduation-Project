%%

clear;
clc;

%%
% excitation signal generation
% MIMO experiment with only odd frequency excited
% generate odd frequency first
% ExcitedHarm is the excited harmonics required for MIMO identification

DefFreq.fs      =   5.5; % sampling frequency of the generator
DefFreq.fres    =   0.01;% frequency spacing in Hz between the (odd) harmonics 
DefFreq.fmin    =   0.1;% lowest excited frequency in Hz (for lowpass design of odd multisines fmin = fres/2)
DefFreq.fmax    =   2;% largest excited frequency in Hz
% DefFreq.frat    =   % ratio between consecutive (odd) harmonics for a logarithmic frequency spacing

Nblock = Inf;                 % one out of three consecutive harmonics is randomly eliminated
Spacing = 'lin';            % linear frequency spacing
MultiType = 'full';          % no even excited harmonics

total_periods = 15; % defination of total periods

[ExcitedHarm, N_odd_lin, NewDefFreq] = HarmMultisine(DefFreq, Nblock, Spacing, MultiType);

fs = NewDefFreq.fs;

% --------------------------------------------------------------------------------------------------
% generate orthogonal excitation

nu = 2;                                     % number of inputs
ny = 3;
N = 5000;                                  % number of time domain samples , each period lasts for 2 seconds and have 40000 data points in total
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
multisine_exp_1 = zeros(2,nh*total_periods);
multisine_exp_2 = zeros(2,nh*total_periods);
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
f = 0.5:1:1999.5;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation known input - noisy output data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load MIMO_PlantNoiseModel;

% Contribution excitation to output %
u0 = multisine_exp_1;
y0 = zeros(ny, L);
for jj=1:ny
    for ii=1:nu
        y0(jj,:) = y0(jj,:) + filter(squeeze(B(jj,ii,:)), A, u0(ii,:));
    end % ii
end % jj

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling the ArbLocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data 
data.u = u0;                            % row index is the input number; the column index the time instant (in samples) 
data.y = y0;                             % row index is the output number; the column index the time instant (in samples) 
data.Ts = 1/fs;                         % sampling period

% method
method.dof = 10;                        % degrees of freedom of the variance estimate
method.order = 2;                       % order local polynomial approximation
method.startfreq = 0.1;              % defines the start frequency of the analysis 
method.stopfreq = 2;                % defines the stop frequency of the analysis
method.step = 1;

% local polynomial estimate FRF and its variance
[CZ, Z, freq, G, CvecG, dof, CL] = ArbLocalPolyAnal(data, method);

% estimated output noise power spectrum
S = CZ.n(1:ny,1:ny,:); 

% estimated variances FRM entries: keep the diagonal elements CvecG only
F = length(freq);
varG = zeros(ny, nu, F);
for kk=1:F
    varG(:, :, kk) = reshape(diag(CvecG(:, :, kk)), [ny, nu]);
end % kk
%%
q = exp(-sqrt(-1)*2*pi*freq/fs);        % z^(-1) as a function of the frequency 

% true FRM
G0 = zeros(ny, nu, F);
AA = polyval(fliplr(A), q.');           % denominator plant transfer function
for ii = 1:ny
	for jj = 1:nu
		G0(ii, jj, :) = polyval(fliplr(squeeze(B(ii, jj, :)).'), q.')./AA;
	end % jj
end % ii

mm = 0;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		plot(freq, db(squeeze(G0(jj,ii,:))), 'k');
	end % ii
end % jj
subplot(ny,nu,1);
title('identified G: red; true G: black');
zoom on; shg

%%

% true noise transfer function
H0 = zeros(ny, ny, F);
DD = polyval(fliplr(D), q.');           % denominator plant transfer function
for ii = 1:ny
	for jj = 1:ny
		H0(ii, jj, :) = polyval(fliplr(squeeze(C(ii, jj, :)).'), q.')./DD;
	end % jj
end % ii

% % true noise power spectrum
% S0 = zeros(ny, ny, F);
% CovE = Te*Te.';                         % covariance matrix of the driving white noise sources 
% for kk = 1:F
% 	S0(:,:,kk) = H0(:,:,kk) * CovE * H0(:,:,kk)';
% 	% remove small imaginary parts on the main diagonal
% 	S0(:,:,kk) = S0(:,:,kk) - diag(sqrt(-1)*imag(diag(S0(:,:,kk))));
% end % frequencies kk
%%
mm = 0;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		plot(freq, db(squeeze(G(jj,ii,:))), 'r', freq, db(squeeze(G0(jj,ii,:))), 'k');
	end % ii
end % jj
subplot(ny,nu,1);
title('identified G: red; true G: black');
zoom on; shg

%%
mm = 0;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		plot(freq, db(squeeze(G(jj,ii,:)-G0(jj,ii,:))), 'k--', freq, db(squeeze(varG(jj,ii,:)))/2, 'r--');
	end % ii
end % jj
subplot(ny,nu,1);
title('|G-G_0|: black --; var(G): red --');
zoom on; shg

%%
Tmatrix = ss([tf(squeeze(B(1,1,:))',A), tf(squeeze(B(1,2,:))',A);
           tf(squeeze(B(2,1,:))',A), tf(squeeze(B(2,2,:))',A);
           tf(squeeze(B(3,1,:))',A), tf(squeeze(B(3,2,:))',A)]);

bodemag(Tmatrix)

Ah = Tmatrix.A; Bh = Tmatrix.B; Ch = Tmatrix.C; Dh = Tmatrix.D;

pole(Tmatrix)
%%
% Interact with simulink model, 6 independent experiments and obtain d,u
% and e for closed-loop identification
clear data;
Ts = 1/fs;
StopTime = length(u0)/fs - Ts;
for i = 1:1
   if i == 1
       EXCSIG = u0;
   end    
   data{i} = sim('scripts_for_model\test_ss.slx');
end

%%
%%
% Extract the data
y_1 = data{1,1}.y.signals.values'; 
d_1 = data{1,1}.d.signals.values'; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calling the ArbLocalPolyAnal function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data 
data_1.u = u0;                            % row index is the input number; the column index the time instant (in samples) 
data_1.y = y_1;                             % row index is the output number; the column index the time instant (in samples) 
data_1.Ts = 1/fs;                         % sampling period

% method
method.dof = 10;                        % degrees of freedom of the variance estimate
method.order = 2;                       % order local polynomial approximation
method.startfreq = 0.1;              % defines the start frequency of the analysis 
method.stopfreq = 2;                % defines the stop frequency of the analysis

% local polynomial estimate FRF and its variance
[CZ, Z, freq, G_2, CvecG, dof, CL] = ArbLocalPolyAnal(data_1, method);

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

G_2 = frd(G_2,freq,FrequencyUnit = 'Hz');
mm = 0;
for jj = 1:ny
	for ii = 1:nu
		mm = mm+1;
		subplot(ny, nu, mm)
		plot(freq, db(squeeze(G_2.ResponseData(jj,ii,:))), 'r', freq, db(squeeze(G0(jj,ii,:))), 'k--');
	end % ii
end % jj
subplot(ny,nu,1);
title('identified G: red; G_true: black;');
zoom on; shg








