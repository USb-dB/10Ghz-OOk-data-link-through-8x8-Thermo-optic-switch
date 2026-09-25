clc; clear; close all;

%% ================= TARGET ===========================
minErrors = 100;
maxBits   = 5e7;

%% ================= SYSTEM ===========================
Rb = 10e9;
Ns = 32;
Fs = Rb*Ns;

PRBS_order = 11;
Nbit_block = 2^PRBS_order - 1;     % 2047 bits
Nblocks = 600;                     % <<< REQUIRED
Nbits_total = Nblocks * Nbit_block;

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

%% ================= OPERATING POINT =================
P_dBm = -5;
P0 = db2pow(P_dBm-30);
fprintf('\nOperating point = %.1f dBm\n',P_dBm);

%% ================= OPTICAL NOISE ===================
OSNR_dB = 28;

%% ================= NN PARAMETERS ===================
Ntaps = 9;
NN_hidden = [32 16];
delay = floor(Ntaps/2);

%% ================= GENERATE 600 PRBS BLOCKS ========
fprintf('Generating %d PRBS-%d blocks...\n',Nblocks,PRBS_order);

txBits_all = zeros(1,Nbits_total);
idx = 1;
for k = 1:Nblocks
    [data,~] = PRBS_generator(1,PRBS_order,1);
    txBits_all(idx:idx+Nbit_block-1) = data(1,:);
    idx = idx + Nbit_block;
end

%% ================= TRANSMITTER ======================
I_tx = conv(repelem(txBits_all,Ns),hRC,'full');
I_tx = max(I_tx,0);
I_tx = I_tx / max(I_tx);

[A_laser,~] = laserCW(LASER,Fs,length(I_tx));
Eout = sqrt(P0*I_tx) .* A_laser;

Ps = mean(abs(Eout).^2);
Pn = Ps / (10^(OSNR_dB/10));
Eout = Eout + sqrt(Pn/2)*(randn(size(Eout))+1j*randn(size(Eout)));

%% ================= PHOTODETECTION ==================
[io,~] = photodiode_PD(Eout,Fs,PD);
v = filtfilt(b,a,io*PD.RL);

%% ================= SAMPLING ========================
gd = (length(hRC)-1)/2;
startIdx = gd + round(Ns/2);
rx = v(startIdx:Ns:startIdx+(Nbits_total-1)*Ns);

%% ================= FIXED SAMPLING BER ===============
thr = (mean(rx(txBits_all==1)) + mean(rx(txBits_all==0)))/2;
rxBits_fixed = rx > thr;
BER_fixed = mean(rxBits_fixed ~= txBits_all);

%% ================= NN TRAINING =====================
fprintf('\nTraining NN using %d bits...\n',Nbits_total);

Xtr = zeros(length(rx),Ntaps);
for n = 1:length(rx)
    ii = max(1,n-delay):min(length(rx),n+delay);
    Xtr(n,delay+1-(n-ii(1)):delay+1+(ii(end)-n)) = rx(ii);
end

mu = mean(Xtr(:));
sg = std(Xtr(:));
Xtr = (Xtr - mu)/sg;

net = fitnet(NN_hidden,'trainscg');
net.performFcn = 'mse';
net.trainParam.epochs = 300;
net.trainParam.showCommandLine = true;
net.trainParam.showWindow = false;

net.divideParam.trainRatio = 0.8;
net.divideParam.valRatio   = 0.2;
net.divideParam.testRatio  = 0;

[net,tr] = train(net,Xtr',txBits_all);

%% ================= NN BER ==========================
rxBits_nn = net(Xtr') > 0.5;
BER_NN = mean(rxBits_nn(:) ~= txBits_all(:));


%% ================= EYE SNR =========================
sig1 = rx(txBits_all==1);
sig0 = rx(txBits_all==0);
EyeSNR_dB = 10*log10((mean(sig1)-mean(sig0))^2/(var(sig1)+var(sig0)));

%% ================= RESULTS =========================
fprintf('\nRESULTS @ -5 dBm (PRBS-%d × %d blocks)\n',PRBS_order,Nblocks);
fprintf('BER fixed = %.3e\n',BER_fixed);
fprintf('BER NN    = %.3e\n',BER_NN);
fprintf('Eye SNR   = %.2f dB\n',EyeSNR_dB);

%% ================= EYE DIAGRAM =====================
figure('Color','w');
eyediagram(v(1:10000)/max(abs(v)),2*Ns,Ns);
title('OOK Eye Diagram @ -5 dBm (600 PRBS-11 blocks)');
