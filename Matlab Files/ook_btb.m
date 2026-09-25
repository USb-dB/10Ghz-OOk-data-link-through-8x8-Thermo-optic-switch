clc; clear; close all;

%% ================= SYSTEM PARAMETERS =================
Rb   = 25e9;                 % Bit rate [Hz]
Ns   = 32;                   % Samples per bit
Fs   = Rb * Ns;              % Sampling frequency
Nbit = 2^11 - 1;             % PRBS length = 2047

%% ================= PRBS ==============================
[data, err_flag] = PRBS_generator(1,11,1);
if err_flag, error('PRBS generation failed'); end
txBits = data(1,:);          % 0 / 1

%% ================= NRZ (INTENSITY SYMBOLS) ===========
txSym = txBits;

%% ================= INTENSITY NYQUIST SHAPING (RC) ====
rollOff = 0.3;
span    = 10;                              % in symbols
hRC = rcosdesign(rollOff, span, Ns, 'normal');

I_tx = conv(repelem(txSym, Ns), hRC, 'full');
I_tx = I_tx / max(I_tx);                   % normalize
I_tx(I_tx < 0) = 0;                        % physical intensity

%% ================= NON-IDEAL LASER ===================
LASER.linewidth = 1e6;                     % Hz
LASER.RIN_dB    = -150;                    % dB/Hz
LASER.phase0    = 0;
LASER.P0_dBm    = 10;

P0 = 10^((LASER.P0_dBm - 30)/10);           % W

[A_laser, ~] = laserCW(LASER, Fs, length(I_tx));

% Optical field (IM/DD abstraction)
Eout = sqrt(P0 * I_tx) .* A_laser;

%% ================= PHOTODIODE (YOUR FUNCTION) ========
PD.R  = 0.9;
PD.BW = 0.7 * Rb;            % used only for noise stats
PD.T  = 305;
PD.RL = 50;

[io, PD] = photodiode_PD(Eout, Fs, PD);
vout = io * PD.RL;           % electrical voltage

%% ================= ELECTRICAL LPF (CORRECT) ==========
BW_rx = 0.75 * Rb;              % MUST be >= (1+rollOff)*Rb/2
LPF_order = 4;

[b,a] = butter(LPF_order, BW_rx/(Fs/2));

% ZERO-PHASE filtering (no group delay!)
vout_lpf = filtfilt(b, a, vout);


%% ================= SYMBOL SAMPLING ===================
gd = (length(hRC)-1)/2;
startIdx = gd + Ns/2;
lastIdx  = startIdx + (Nbit-1)*Ns;

if lastIdx > length(vout_lpf)
    error('Not enough samples after LPF');
end

rxSym = vout_lpf(startIdx:Ns:lastIdx);

%% ================= THRESHOLD DECISION =================
thr = (mean(rxSym(txBits==1)) + mean(rxSym(txBits==0))) / 2;
rxBits = rxSym > thr;

%% ================= BER ===============================
BER = mean(rxBits ~= txBits);

fprintf('Bits transmitted = %d\n', length(txBits));
fprintf('Bits received    = %d\n', length(rxBits));
fprintf('BER              = %.3e\n', BER);

%% ================= EYE DIAGRAM =======================
% Prepare eye-diagram signal (oversampled, after LPF)
s_eye_raw = vout_lpf(startIdx-5*Ns : startIdx+200*Ns);
s_eye_plot = s_eye_raw / max(abs(s_eye_raw));

figure('Color','w');
eyediagram(s_eye_plot, 2*Ns, Ns);
title('OOK Eye Diagram (Experimental-Grade)');
