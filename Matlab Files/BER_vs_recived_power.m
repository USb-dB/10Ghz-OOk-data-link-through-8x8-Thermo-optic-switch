clc; clear; close all;

%% ================= USER TARGET ======================
BER_target_low  = 1e-5;
BER_target_high = 1e-4;

minErrors = 100;
maxBits   = 2e8;

SHOW_EYE = true;     % <<< force eye diagram

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
P_dBm_vec = -5;     % <<< single power point
P_store   = [];
BER_store = [];

found = false;

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

        %% ---- Photodiode ----
        [io, ~] = photodiode_PD(Eout, Fs, PD);
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

    if BER_meas >= BER_target_low && BER_meas <= BER_target_high
        found  = true;
        P_op   = P_dBm;
        BER_op = BER_meas;
    end
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
if SHOW_EYE

    % --- Use first simulated power point ---
    P_eye   = P_store(1);
    BER_eye = BER_store(1);

    %% ---- PRBS ----
    [data, ~] = PRBS_generator(1, PRBS_order, 1);
    txBits = data(1,:);

    %% ---- TX waveform ----
    I_tx = conv(repelem(txBits, Ns), hRC, 'full');
    I_tx = I_tx / max(I_tx);
    I_tx(I_tx < 0) = 0;

    %% ---- Laser ----
    [A_laser, ~] = laserCW(LASER, Fs, length(I_tx));
    Eout = sqrt(10^((P_eye-30)/10) * I_tx) .* A_laser;

    %% ---- Photodiode ----
    [io, ~] = photodiode_PD(Eout, Fs, PD);
    vout = io * PD.RL;

    %% ---- RX LPF ----
    vout = filtfilt(b, a, vout);

    %% ---- Sampling phase (group delay aligned) ----
    gd = (length(hRC)-1)/2;
    startIdx = gd + Ns/2;

    %% ---- AGC (based on symbol means) ----
    sym_tmp = vout(startIdx:Ns:startIdx+(Nbit_block-1)*Ns);
    A1 = mean(sym_tmp(txBits==1));
    A0 = mean(sym_tmp(txBits==0));
    vout = vout / (A1 - A0);

    %% ---- Eye extraction (SAFE indexing) ----
    Neye_sym = 200;                 % number of symbols in eye
    eye_start = startIdx;
    eye_end   = startIdx + Neye_sym*Ns - 1;
    eye_end   = min(eye_end, length(vout));

    s_eye = vout(eye_start:eye_end);

    %% ---- Plot eye diagram ----
    figure('Color','w');
    eyediagram(s_eye, 2*Ns, Ns);
    grid on;
    title(sprintf('OOK Eye Diagram @ %.1f dBm  (BER ≈ %.2e)', ...
                  P_eye, BER_eye));
end
