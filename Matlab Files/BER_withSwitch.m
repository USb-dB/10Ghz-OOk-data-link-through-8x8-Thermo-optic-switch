 clc; clear; close all;


%% ===================== DATA =====================

clc; clear; close all;

%% ======================================================
%                USER FILE
%% ======================================================
filename = 'TE0.dat';

%% ======================================================
%                LOAD DATA
%% ======================================================
T = readtable(filename, ...
    'Delimiter','\t', ...
    'FileType','text', ...
    'VariableNamingRule','preserve');

% Column names (as in your file)
wav_nm = T.("wavelength(1-1)");
TX     = T.("TX(1-7)");
phi    = T.("Angle(1-7)");

%% ======================================================
%                CLEAN DATA
%% ======================================================
valid = isfinite(wav_nm) & isfinite(TX) & isfinite(phi);
wav_nm = wav_nm(valid);
TX     = TX(valid);
phi    = phi(valid);

%% ======================================================
%                CONSTANTS
%% ======================================================
c = 3e8;
lambda0 = 1550e-9;
omega0  = 2*pi*c/lambda0;

%% ======================================================
%          WAVELENGTH → ANGULAR FREQUENCY
%% ======================================================
lambda = wav_nm * 1e-9;
omega  = 2*pi*c ./ lambda;

% Sort
[omega, idx] = sort(omega);
TX  = TX(idx);
phi = phi(idx);

% Enforce unique grid (MANDATORY)
[omega, uidx] = unique(omega,'stable');
TX  = TX(uidx);
phi = phi(uidx);
phi = phi/180;

%% ======================================================
%            BASE COMPLEX TRANSFER FUNCTION
%% ======================================================
H_raw  = TX .* exp(1j*phi);
Omega  = omega - omega0;

%% ======================================================
%          EXTRAPOLATION MODEL (PHYSICALLY SAFE)
%% ======================================================
% Magnitude and unwrapped phase
mag = abs(H_raw);
phi_u = unwrap(angle(H_raw));

% Group delay estimate
tau_g = mean(gradient(phi_u, Omega));

%% ======================================================
%         BUILD EXTRAPOLATED H(Ω) FUNCTION HANDLE
%% ======================================================
H_of_Omega = @(Om) ...
    interp1(Omega, mag, Om, 'linear', 'extrap') .* ...
    exp(1j * ...
        ( ...
        interp1(Omega, phi_u, Om, 'linear', 'extrap') + ...
        (Om < min(Omega)) .* tau_g .* (Om - min(Omega)) + ...
        (Om > max(Omega)) .* tau_g .* (Om - max(Omega)) ...
        ));

%% ======================================================
%           FINAL OUTPUT (DISCRETE GRID)
%% ======================================================
% You may choose any grid you want here
Omega_out = Omega;             % original grid
H = H_of_Omega(Omega_out);     % FINAL TRANSFER FUNCTION

%% ======================================================
%                EXPORT
%% ======================================================


fprintf('H generated successfully (%d points)\n', length(H));

%% ======================================================
%                SANITY CHECK PLOTS
%% ======================================================
figure('Color','w');

subplot(2,1,1)
plot(Omega_out/2/pi/1e9, abs(H),'LineWidth',1.5)
xlabel('\Omega / 2\pi (GHz)')
ylabel('|H|')
grid on
title('Magnitude Response')

subplot(2,1,2)
plot(Omega_out/2/pi/1e9, unwrap(angle(H)),'LineWidth',1.5)
xlabel('\Omega / 2\pi (GHz)')
ylabel('Phase (rad)')
grid on
title('Unwrapped Phase')



%% ================= USER TARGET ======================
BER_target_low  = 1e-5;
BER_target_high = 1e-4;

minErrors = 100;
maxBits   = 2e8;

SHOW_EYE = true;

%% ================= SYSTEM PARAMETERS =================
Rb   = 10e9;
Ns   = 32;
Fs   = Rb * Ns;

PRBS_order = 11;
Nbit_block = 2^PRBS_order - 1;

%% ================= TX: INTENSITY NYQUIST (RC) ========
rollOff = 0.3;
span    = 10;
hRC = rcosdesign(rollOff, span, Ns, 'normal');

%% ================= LASER PARAMETERS ==================
LASER.linewidth = 1e6;
LASER.RIN_dB    = -150;
LASER.phase0    = 0;

%% ================= PHOTODIODE PARAMETERS =============
PD.R  = 0.9;
PD.BW = 0.35 * Rb;
PD.T  = 300;
PD.RL = 5e4;

%% ================= ELECTRICAL LPF ====================
BW_rx = 0.75 * Rb;
[b,a] = butter(4, BW_rx/(Fs/2));

%% ================= POWER =============================
P_dBm_vec = -5;
P_store   = [];
BER_store = [];

%% ================= EYE STORAGE (NEW) =================
EYE_CAPTURED = false;                 %%%
EYE_BUF = [];                          %%%
Neye_sym = 200;                       %%%

%% ================= MAIN LOOP =========================
for P_dBm = P_dBm_vec

    fprintf('\nTesting P = %.1f dBm\n', P_dBm);

    P0 = 10^((P_dBm - 30)/10);
    LASER.P0_dBm = P_dBm;

    totalErr  = 0;
    totalBits = 0;
    blocks    = 0;

    while totalErr < minErrors && totalBits < maxBits

        %% ---- PRBS ----
        [data, err_flag] = PRBS_generator(1,PRBS_order,1);
        if err_flag, error('PRBS failed'); end
        txBits = data(1,:);

        %% ---- TX waveform ----
        I_tx = conv(repelem(txBits,Ns), hRC, 'full');
        I_tx = I_tx / max(I_tx);
        I_tx(I_tx < 0) = 0;

        %% ---- Laser ----
        [A_laser, ~] = laserCW(LASER, Fs, length(I_tx));
        Eout = sqrt(P0 * I_tx) .* A_laser;

        % --- SWITCH ---
        N = length(Eout);
        f = (-N/2:N/2-1)*(Fs/N);
        Omega_fft = 2*pi*f;
        Hsw = interp1(Omega,H,Omega_fft,'linear',0);
        Hsw = fftshift(Hsw);
    
        E = ifft(fft(Eout).*Hsw);


        %% ---- Photodiode ----
        [io, ~] = photodiode_PD(E, Fs, PD);
        vout = io * PD.RL;

        %% ---- RX LPF ----
        vout_lpf = filtfilt(b,a,vout);

        %% ---- Sampling phase ----
        gd = (length(hRC)-1)/2;
        startIdx = gd + Ns/2;

        sym_tmp = vout_lpf(startIdx:Ns:startIdx+(Nbit_block-1)*Ns);

        %% ---- AGC (signal only) ----
        A1 = mean(sym_tmp(txBits==1));
        A0 = mean(sym_tmp(txBits==0));
        vout_lpf = vout_lpf / (A1 - A0);

        %% ---- CAPTURE EYE ONCE (FIRST SNAPSHOT) ----
        if SHOW_EYE && ~EYE_CAPTURED
            eye_start = startIdx;
            eye_end   = min(startIdx + Neye_sym*Ns - 1, length(vout_lpf));
            EYE_BUF   = vout_lpf(eye_start:eye_end);
            EYE_CAPTURED = true;
        end

        %% ---- Decisions ----
        rxSym  = vout_lpf(startIdx:Ns:startIdx+(Nbit_block-1)*Ns);
        thr    = (mean(rxSym(txBits==1)) + mean(rxSym(txBits==0))) / 2;
        rxBits = rxSym > thr;

        %% ---- Errors ----
        errs = sum(rxBits ~= txBits);
        totalErr  = totalErr  + errs;
        totalBits = totalBits + Nbit_block;
        blocks    = blocks + 1;
    end

    BER_meas = totalErr / totalBits;

    fprintf('  Blocks=%d  Bits=%.2e  Errors=%d  BER=%.3e\n', ...
            blocks, totalBits, totalErr, BER_meas);

    P_store(end+1)   = P_dBm;
    BER_store(end+1) = BER_meas;
end

%% ================= BER PLOT ==========================
figure('Color','w');
semilogy(P_store, BER_store,'o-','LineWidth',1.8);
grid on;
xlabel('Received Optical Power (dBm)');
ylabel('Measured BER');
title('IM/DD OOK Receiver Sensitivity');
ylim([1e-6 1]);

yline(3.8e-3,'--r','HD-FEC');
yline(2e-4,'--b','KP4-FEC');

%% ================= EYE DIAGRAM =======================
if SHOW_EYE && ~isempty(EYE_BUF)
    figure('Color','w');
    eyediagram(EYE_BUF, 2*Ns, Ns);
    grid on;
    title(sprintf('OOK Eye Diagram @ %.1f dBm (BER ≈ %.2e)', ...
                  P_store(1), BER_store(1)));
end
