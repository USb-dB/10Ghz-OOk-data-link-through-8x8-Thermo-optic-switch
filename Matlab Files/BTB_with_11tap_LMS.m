clc; clear; close all;

%% ================= TARGET ===========================
minErrors = 100;
maxBits   = 5e7;

%% ================= SYSTEM ===========================
Rb = 10e9;
Ns = 32;
Fs = Rb*Ns;

PRBS_order = 11;
Nbit_block = 2^PRBS_order - 1;

%% ================= TX FILTER ========================
rollOff = 0.3;
span = 10;
hRC = rcosdesign(rollOff,span,Ns,'normal');

%% ================= LASER ============================
LASER.linewidth = 1e6;
LASER.RIN_dB = -150;
LASER.phase0 = 0;

%% ================= PHOTODIODE ======================
PD.R  = 0.9;
PD.BW = 0.30*Rb;
PD.RL = 5e4;
PD.T  = 300;

%% ================= RX LPF ==========================
BW_rx = 0.6*Rb;
[b,a] = butter(4,BW_rx/(Fs/2));

%% ================= OPTICAL NOISE ===================
OSNR_dB = 28;   % realistic, recoverable (adjust if needed)

%% ================= POWER SWEEP =====================
P_dBm_vec = -15:0.5:-5;

BER_fixed  = zeros(size(P_dBm_vec));
BER_opt    = zeros(size(P_dBm_vec));
SNR_eye_dB = zeros(size(P_dBm_vec));

%% ================= MAIN LOOP =======================
for p = 1:length(P_dBm_vec)

    P_dBm = P_dBm_vec(p);
    fprintf('\nTesting P = %.1f dBm\n',P_dBm);

    LASER.P0_dBm = P_dBm;
    P0 = db2pow(P_dBm-30);

    %% -------- FIXED SAMPLING ------------------------
    totalErr = 0; totalBits = 0;

    while totalErr < minErrors && totalBits < maxBits

        [data,~] = PRBS_generator(1,PRBS_order,1);
        txBits = data(1,:);

        I_tx = conv(repelem(txBits,Ns),hRC,'full');
        I_tx = max(I_tx,0);
        I_tx = I_tx/max(I_tx);

        [A_laser,~] = laserCW(LASER,Fs,length(I_tx));
        Eout = sqrt(P0*I_tx).*A_laser;

        % ===== OPTICAL AWGN (ASE-LIKE) =====
        Ps = mean(abs(Eout).^2);
        Pn = Ps / (10^(OSNR_dB/10));
        Eout = Eout + sqrt(Pn/2)*(randn(size(Eout))+1j*randn(size(Eout)));

        [io,~] = photodiode_PD(Eout,Fs,PD);
        v = filtfilt(b,a,io*PD.RL);

        gd = (length(hRC)-1)/2;
        startIdx = gd + Ns/2;

        rx = v(startIdx:Ns:startIdx+(Nbit_block-1)*Ns);

        thr = (mean(rx(txBits==1))+mean(rx(txBits==0)))/2;
        rxBits = rx > thr;

        totalErr  = totalErr  + sum(rxBits~=txBits);
        totalBits = totalBits + Nbit_block;
    end

    BER_fixed(p) = totalErr/totalBits;

    %% -------- OPTIMIZED SAMPLING --------------------
    bestBER   = 1;
    bestPhase = 1;

    for phase = 1:Ns

        totalErr = 0; totalBits = 0;

        while totalErr < minErrors && totalBits < maxBits/10

            [data,~] = PRBS_generator(1,PRBS_order,1);
            txBits = data(1,:);

            I_tx = conv(repelem(txBits,Ns),hRC,'full');
            I_tx = max(I_tx,0);
            I_tx = I_tx/max(I_tx);

            [A_laser,~] = laserCW(LASER,Fs,length(I_tx));
            Eout = sqrt(P0*I_tx).*A_laser;

            % ===== OPTICAL AWGN =====
            Ps = mean(abs(Eout).^2);
            Pn = Ps / (10^(OSNR_dB/10));
            Eout = Eout + sqrt(Pn/2)*(randn(size(Eout))+1j*randn(size(Eout)));

            [io,~] = photodiode_PD(Eout,Fs,PD);
            v = filtfilt(b,a,io*PD.RL);

            gd = (length(hRC)-1)/2;
            startIdx = gd + phase;

            rx = v(startIdx:Ns:startIdx+(Nbit_block-1)*Ns);

            thr = (mean(rx(txBits==1))+mean(rx(txBits==0)))/2;
            rxBits = rx > thr;

            totalErr  = totalErr  + sum(rxBits~=txBits);
            totalBits = totalBits + Nbit_block;
        end

        BERp = totalErr/totalBits;
        if BERp < bestBER
            bestBER   = BERp;
            bestPhase = phase;
        end
    end

    BER_opt(p) = bestBER;

    %% -------- EYE-SNR -------------------------------
    [data,~] = PRBS_generator(1,PRBS_order,1);
    txBits = data(1,:);

    I_tx = conv(repelem(txBits,Ns),hRC,'full');
    I_tx = max(I_tx,0);
    I_tx = I_tx/max(I_tx);

    [A_laser,~] = laserCW(LASER,Fs,length(I_tx));
    Eout = sqrt(P0*I_tx).*A_laser;

    Ps = mean(abs(Eout).^2);
    Pn = Ps / (10^(OSNR_dB/10));
    Eout = Eout + sqrt(Pn/2)*(randn(size(Eout))+1j*randn(size(Eout)));

    [io,~] = photodiode_PD(Eout,Fs,PD);
    v = filtfilt(b,a,io*PD.RL);

    gd = (length(hRC)-1)/2;
    startIdx = gd + bestPhase;
    rx = v(startIdx:Ns:startIdx+(Nbit_block-1)*Ns);

    sig1 = rx(txBits==1);
    sig0 = rx(txBits==0);
    eyeSNR = (mean(sig1)-mean(sig0))^2/(var(sig1)+var(sig0));
    SNR_eye_dB(p) = 10*log10(eyeSNR);

    fprintf('  BER fixed = %.3e | BER opt = %.3e | Eye-SNR = %.2f dB\n', ...
            BER_fixed(p),BER_opt(p),SNR_eye_dB(p));
end

%% ================= PLOTS ===========================
figure('Color','w');
semilogy(P_dBm_vec,BER_fixed,'--o','LineWidth',1.6); hold on;
semilogy(P_dBm_vec,BER_opt,'-s','LineWidth',2);
grid on;
xlabel('Received Optical Power (dBm)');
ylabel('Measured BER');
title('BER with Optical AWGN (Same Testbench)');
legend('Fixed sampling','Optimized sampling','Location','SouthWest');
ylim([1e-6 1]);

figure('Color','w');
plot(P_dBm_vec,SNR_eye_dB,'-o','LineWidth',2);
grid on;
xlabel('Received Optical Power (dBm)');
ylabel('Eye SNR (dB)');
title('Eye-SNR vs Power with Optical Noise');
