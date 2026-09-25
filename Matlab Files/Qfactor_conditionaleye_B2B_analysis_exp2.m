clc; clear; close all;

 
% SYSTEM PARAMETERS
 
Rb  = 25e9;                 % Bit rate [bit/s]
Ns  = 32;                   % Samples per bit
Fs  = Rb * Ns;              % Sampling frequency

 
% PRBS PARAMETERS
 
num_PRBS = 1;
deg      = 11;              % PRBS order
num_pol  = 1;

 
% MZM PARAMETERS (OOK)
 
Vpi   = 5;
Vbias = Vpi/2;              % Quadrature bias

 
% LASER PARAMETERS (EXPERIMENTAL)
 
LASER.linewidth = 1e6;      
LASER.RIN_dB    = -150;     % realistic commercial laser
LASER.phase0    = 0;
LASER.P0_dBm    = 10;

 
% PHOTODIODE PARAMETERS (NON-IDEAL)
 
PD.R  = 0.9;
PD.BW = 0.7 * Rb;           % IMPORTANT: non-ideal BW
PD.T  = 300;
PD.RL = 50;

 
% PRBS GENERATION
 
[data, err_flag] = PRBS_generator(num_PRBS, deg, num_pol);
if err_flag, error('PRBS generation failed'); end
bits = data(1,:);

 
% TRANSMITTER PULSE SHAPING (RRC)
 
symbols = bits - mean(bits);

PS.type     = 'RRC';
PS.rollOff  = 0.3;
PS.nTaps    = 256*Ns;

[nrz_shaped, ~] = pulseShaper(symbols, Ns, PS);
nrz_shaped = nrz_shaped / max(abs(nrz_shaped));
nrz = 0.5*(nrz_shaped + 1);

 
% MZM DRIVE
 
Vdrive = (Vpi/2) * (2*nrz - 1);

 
% LASER
 
[A_laser, ~] = laserCW(LASER, Fs, length(nrz));

 
% MZM MODULATION
 
Eout = A_laser .* cos((pi/(2*Vpi))*(Vdrive + Vbias));

 
% PHOTODIODE
 
[io, PD] = photodiode_PD(Eout, Fs, PD);
vout = io * PD.RL;

 
% MATCHED FILTER (RRC)
 
LPF.type           = 'RRC';
LPF.rollOff        = PS.rollOff;
LPF.nTaps          = 12*Ns;
LPF.implementation = 'conv';

[S_eye, ~] = LPF_apply(vout.', LPF, Fs, Rb);
s_eye = S_eye(:).';

 
% PREPARE SIGNAL
 
s_eye_raw = s_eye(20*Ns:end);
s_eye_plot = s_eye_raw / max(abs(s_eye_raw));

 
% EYE DIAGRAM
 
figure('Color','w');
eyediagram(s_eye_plot, 2*Ns, Ns);
title('OOK Eye Diagram (Experimental-Grade)');

 
% OFFLINE DSP
n_iter = 5;

DSP = offline_DSP_OOK( ...
        s_eye_raw, ...
        bits, ...
        Fs, ...
        Rb, ...
        n_iter );

% RESULTS
fprintf('\n--- Eye Metrics ---\n');
fprintf('Eye height        = %.4f V\n', DSP.eye_height);
fprintf('Eye width         = %.3f UI\n', DSP.eye_width_UI);
fprintf('Timing jitter RMS = %.4e UI\n', DSP.jitter_rms_UI);

fprintf('\n--- Performance Metrics ---\n');
fprintf('Q-factor          = %.3f\n', DSP.Q);
fprintf('BER (measured)    = %.3e\n', DSP.BER_measured);
fprintf('BER (Gaussian)    = %.3e\n', DSP.BER_gaussian);
