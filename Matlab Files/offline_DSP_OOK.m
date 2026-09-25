function DSP = offline_DSP_OOK(rx_osc, tx_bits, Fs_osc, Rs, n_iter)

%% ================== 1. RESAMPLING ==================
Fs_dsp = 2 * Rs;                              % 2 samples/symbol
[rx_rs, ~] = applyResample(rx_osc, Fs_osc, Fs_dsp);

rx_rs = rx_rs(:).';
Ns    = round(Fs_dsp / Rs);                   % Ns = 2

%% ================== 1b. SAMPLING JITTER (REALISTIC) ==================
sigma_jitter_UI = 0.01;                       % ~1% UI RMS
rx_rs = addSamplingJitter(rx_rs, sigma_jitter_UI, Ns);

%% ================== 2. TIMING RECOVERY (Gardner TED) ==================
mu_tr = 0.01;
tau   = 0;

Nsym  = floor(length(rx_rs) / Ns);
rx_tr = builtin('zeros', 1, Nsym);

for k = 2 : Nsym-1
    idx = (k-1)*Ns + round(tau) + 1;

    if idx-1 < 1 || idx+1 > length(rx_rs)
        continue;
    end

    xE = rx_rs(idx-1);
    xM = rx_rs(idx);
    xL = rx_rs(idx+1);

    e   = (xL - xE) * xM;
    tau = tau + mu_tr * e;

    rx_tr(k) = xM;
end

rx_tr   = rx_tr(3:end);
tx_bits = tx_bits(1:length(rx_tr));

%% ================== 3. n-th ORDER THRESHOLD DETECTION ==================
thr = mean(rx_tr);

for ii = 1:n_iter
    s1 = rx_tr(rx_tr > thr);
    s0 = rx_tr(rx_tr <= thr);

    mu1 = mean(s1);
    mu0 = mean(s0);

    thr = (mu1 + mu0)/2;
end

sigma1 = std(s1);
sigma0 = std(s0);

%% ================== 3b. THRESHOLD DRIFT (EXPERIMENTAL) ==================
thr = thr + 0.03 * thr * randn;   % ~3% slow drift

%% ================== 4. HARD DECISION ==================
rx_bits = rx_tr > thr;

%% ================== 5. BIT ALIGNMENT (MANDATORY) ==================
[c,lags] = xcorr(double(rx_bits), double(tx_bits));
[~,idx]  = max(c);
lag      = lags(idx);

if lag > 0
    rx_bits = rx_bits(1+lag:end);
    tx_bits = tx_bits(1:length(rx_bits));
elseif lag < 0
    tx_bits = tx_bits(1-lag:end);
    rx_bits = rx_bits(1:length(tx_bits));
end

%% ================== 6. BER (Measured) ==================
BER_measured = mean(rx_bits ~= tx_bits);

%% ================== 7. Q FACTOR ==================
Q = (mu1 - mu0) / (sigma1 + sigma0);

%% ================== 8. BER (Gaussian Estimate) ==================
BER_gaussian = 0.5 * erfc(Q / sqrt(2));

%% ================== 9. EYE METRICS (EXPERIMENTAL) ==================
rx_rs_eye = rx_rs(1 : floor(length(rx_rs)/Ns)*Ns);
rx_eye    = reshape(rx_rs_eye, Ns, []);

eye_height = mu1 - mu0;

% Percentile-based eye width (robust, experimental)
upper = prctile(rx_eye, 90, 2);
lower = prctile(rx_eye, 10, 2);
eye_open = upper - lower;

idx = find(eye_open > 0.5*max(eye_open));
eye_width_UI = (idx(end) - idx(1)) / Ns;

%% ================== 10. JITTER ESTIMATION ==================
jitter_samples = [];

for n = 1:size(rx_eye,2)
    w = rx_eye(:,n);
    c = find(diff(w > thr) ~= 0);
    if ~isempty(c)
        jitter_samples(end+1,1) = c - round(Ns/2); %#ok<SAGROW>
    end
end

jitter_rms_UI = std(jitter_samples) / Ns;

%% ================== 11. OUTPUT ==================
DSP.mu1 = mu1;
DSP.mu0 = mu0;
DSP.threshold = thr;

DSP.eye_height   = eye_height;
DSP.eye_width_UI = eye_width_UI;

DSP.Q             = Q;
DSP.BER_measured  = BER_measured;
DSP.BER_gaussian  = BER_gaussian;

DSP.jitter_rms_UI = jitter_rms_UI;
DSP.rx_bits       = rx_bits;

end
