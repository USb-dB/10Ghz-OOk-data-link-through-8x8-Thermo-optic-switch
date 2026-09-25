clc; clear; close all;

%% ================= PARAMETERS =================
Rb   = 25e9;
Ns   = 32;
Fs   = Rb*Ns;
Nbit = 2^11 - 1;

%% ================= PRBS ======================
[data, err_flag] = PRBS_generator(1,11,1);
if err_flag, error('PRBS failed'); end
txBits = data(1,:);   % 0/1

%% ================= NRZ =======================
txSym = txBits;

%% ================= INTENSITY NYQUIST (RC) ===
rollOff = 0.3;
span    = 10;                     % in symbols
hRC = rcosdesign(rollOff, span, Ns, 'normal');

I_tx = conv(repelem(txSym, Ns), hRC, 'full');
I_tx = I_tx / max(I_tx);          % normalize
I_tx(I_tx < 0) = 0;               % physical intensity

%% ================= OPTICAL FIELD ============
Etx = sqrt(I_tx);

%% ================= PHOTODIODE ===============
I_rx = abs(Etx).^2;

%% ================= SYMBOL SAMPLING ==========
gd = (length(hRC)-1)/2;
startIdx = gd + Ns/2;
rxSym = I_rx(startIdx:Ns:startIdx+(Nbit-1)*Ns);

%% ================= DECISION =================
thr = (max(rxSym) + min(rxSym))/2;
rxBits = rxSym > thr;

%% ================= BER ======================

%% ================= BER ===============================
BER = BER_eval(txBits, rxBits);
fprintf('Bits transmitted = %d\n', length(txBits));
fprintf('Bits received    = %d\n', length(rxBits));
fprintf('BER              = %.3e\n', BER);
