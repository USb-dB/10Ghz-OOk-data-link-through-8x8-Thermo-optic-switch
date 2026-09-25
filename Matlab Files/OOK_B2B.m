clc; clear; close all;

%% =============================================================
%  SYSTEM PARAMETERS
%% =============================================================
Rb  = 25e9;            % Bit rate [bit/s]
Tb  = 1/Rb;
Ns  = 32;              % Samples per bit
Fs  = Rb * Ns;         % Sampling frequency

%% =============================================================
%  PRBS PARAMETERS
%% =============================================================
num_PRBS = 1;
deg      = 13;         % Long PRBS for eye statistics
num_pol  = 1;

%% =============================================================
%  MZM PARAMETERS (OOK)
%% =============================================================
Vpi   = 5;             % Half-wave voltage [V]
Vbias = 0;             % NULL bias (critical for OOK)

%% =============================================================
%  LASER PARAMETERS
%% =============================================================
LASER.linewidth = 1e6;     % Laser linewidth [Hz]
LASER.RIN_dB    = -160;
LASER.phase0    = 0;
LASER.P0_dBm    = 10;      % CW power [dBm]

%% =============================================================
%  PHOTODIODE PARAMETERS
%% =============================================================
PD.R  = 0.9;            % Responsivity [A/W]
PD.BW = Rb;             % Electrical bandwidth
PD.T  = 300;            % Temperature [K]
PD.RL = 50;             % Load resistance [Ohm]

%% =============================================================
%  PRBS GENERATION
%% =============================================================
[data, err_flag] = PRBS_generator(num_PRBS, deg, num_pol);
if err_flag
    error('PRBS generation failed');
end
bits = data(1,:);                       % PRBS sequence
Nbit = length(bits);

%% =============================================================
%  ELECTRICAL PULSE SHAPING (TRANSMITTER)
%% =============================================================
symbols = bits;                         % OOK symbols (0/1)
symbols = symbols - mean(symbols);      % Zero-mean for shaping

PS.type     = 'RRC';
PS.rollOff = 0.7;                       % Small roll-off
PS.nTaps   = 256*Ns;                    % Long filter

[nrz_shaped, PS] = pulseShaper(symbols, Ns, PS);

% Normalize & map back to unipolar OOK
nrz_shaped = nrz_shaped / max(abs(nrz_shaped));
nrz = 0.5*(nrz_shaped + 1);             % Range: [0,1]

t = (0:length(nrz)-1)/Fs;

%% =============================================================
%  MZM DRIVE VOLTAGE
%% =============================================================
Vdrive = Vpi * nrz;                     % 0 → Vπ swing

%% =============================================================
%  CW LASER
%% =============================================================
nSamples = length(nrz);
[A_laser, LASER] = laserCW(LASER, Fs, nSamples);

%% =============================================================
%  MACH–ZEHNDER MODULATION
%% =============================================================
Eout = A_laser .* cos( (pi/(2*Vpi)) * (Vdrive + Vbias) );
Pout = abs(Eout).^2;

%% =============================================================
%  PHOTODIODE DETECTION
%% =============================================================
[io, PD] = photodiode_PD(Eout, Fs, PD);
vout = io * PD.RL;

%% =============================================================
%  RECEIVER ELECTRICAL LPF
%% =============================================================
LPF.type  = 'Bessel';
LPF.order = 5;
LPF.fc    = 0.75 * Rb;

[S_eye, LPF] = LPF_apply(vout.', LPF, Fs, Rb);
s_eye = S_eye(:).';

%% =============================================================
%  TIME-DOMAIN RECEIVER SIGNAL
%% =============================================================
figure('Color','w','Position',[200 100 900 350])
plot(t*1e9, s_eye, 'LineWidth',1.1)
xlabel('Time (ns)')
ylabel('Voltage (V)')
title('Receiver Electrical Signal After LPF')
grid on

%% =============================================================
%  EYE DIAGRAM (PROPER OOK EYE)
%% =============================================================
s_eye_dc   = s_eye - mean(s_eye);       % Remove DC
nDiscard   = 20*Ns;                     % Remove LPF transient
s_eye_trim = s_eye_dc(nDiscard:end);
s_eye_norm = s_eye_trim / max(s_eye_trim);

figure('Color','w')
eyediagram(s_eye_norm(:), 2*Ns, Ns);
title('OOK Eye Diagram (Pulse Shaped, Proper Eye Opening)')

% using eye diagram dense


s_eye_dc   = s_eye - mean(s_eye);
s_eye_trim = s_eye_dc(20*Ns:end);
s_eye_norm = s_eye_trim / max(abs(s_eye_trim));

rep = Ns;

[fh, scr] = eyediagram_dense( ...
    s_eye_norm, ...
    2*rep, ...                 % 2 UI
    "midbit", true, ...
    "histogram", true, ...     % <<< IMPORTANT
    "xresolution", 900, ...
    "yresolution", 600);

title('OOK Eye Diagram (Density / Persistence View)');




%% =============================================================
%  TRANSMITTER & OPTICAL SIGNAL PLOTS
%% =============================================================
figure('Color','w','Position',[200 100 900 700])

subplot(4,1,1)
stairs(bits,'LineWidth',1.2)
title('PRBS Bits'); ylim([-0.2 1.2]); grid on

subplot(4,1,2)
plot(t*1e9,nrz,'LineWidth',1.2)
title('Pulse-Shaped OOK Electrical Signal')
ylabel('Amplitude'); grid on

subplot(4,1,3)
plot(t*1e9,Vdrive,'LineWidth',1.2)
title('MZM Drive Voltage')
ylabel('Voltage (V)'); grid on

subplot(4,1,4)
plot(t*1e9,Pout,'LineWidth',1.2)
title('Optical OOK Power')
xlabel('Time (ns)'); ylabel('Power'); grid on

%% =============================================================
%  MZM TRANSFER FUNCTION
%% =============================================================
V = linspace(0,Vpi,1000);
T = cos(pi*(V+Vbias)/(2*Vpi)).^2;

figure('Color','w')
plot(V,T,'LineWidth',2)
xlabel('Voltage (V)')
ylabel('Normalized Optical Power')
title('MZM Transfer Function (OOK, Null Bias)')
grid on

%% =============================================================
%  NOISE ANALYSIS
%% =============================================================
fprintf('\n--- Photodiode Noise Statistics ---\n');
fprintf('Mean signal current    = %.3e A\n', PD.meanCurrent);
fprintf('Shot noise variance    = %.3e A^2\n', PD.shotNoiseVar);
fprintf('Thermal noise variance = %.3e A^2\n', PD.thermalNoiseVar);
fprintf('Total noise variance   = %.3e A^2\n', ...
        PD.shotNoiseVar + PD.thermalNoiseVar);

i_signal = PD.R * abs(Eout).^2;
i_noise  = io - i_signal;

figure('Color','w','Position',[200 100 900 450])

subplot(2,1,1)
plot(t*1e9, i_signal, 'LineWidth',1.1)
title('Photodiode Signal Current')
ylabel('Current (A)'); grid on

subplot(2,1,2)
plot(t*1e9, i_noise, 'LineWidth',1.1)
title('Photodiode Noise Current')
xlabel('Time (ns)'); ylabel('Current (A)'); grid on

figure('Color','w')
histogram(i_noise,100,'Normalization','pdf')
xlabel('Noise Current (A)')
ylabel('PDF')
title('Photodiode Noise Distribution')
grid on
