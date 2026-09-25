function [A,LASER] = laserCW(LASER,Fs,nSamples)

%% Defaults
if ~isfield(LASER,'linewidth'), LASER.linewidth = 0; end
if ~isfield(LASER,'RIN_dB'),    LASER.RIN_dB = -inf; end
if ~isfield(LASER,'phase0'),    LASER.phase0 = 0; end
if ~isfield(LASER,'P0_dBm'),    LASER.P0_dBm = 30; end

lw      = LASER.linewidth;
RIN_dB  = LASER.RIN_dB;
ph0     = LASER.phase0;
P0_dBm  = LASER.P0_dBm;

%% Phase noise (Wiener)
phVar = 2*pi*lw/Fs;
phNoise = cumsum(sqrt(phVar)*randn(1,nSamples));

%% Optical power
P0 = db2pow(P0_dBm-30);

%% RIN (frequency-shaped, experimental)
intVar = 10^(RIN_dB/10) * Fs * P0^2;

white = randn(1,nSamples);
rinNoise = filter(1, [1 -0.995], white);
rinNoise = rinNoise / std(rinNoise);

intNoise = sqrt(intVar) * rinNoise;

%% Instantaneous power (clamped)
Pinst = max(P0 + intNoise, 0.1*P0);

%% Optical field
A = sqrt(Pinst) .* exp(1j*(ph0 + phNoise));

%% Outputs
LASER.phaseVar   = phVar;
LASER.phaseNoise = phNoise;
LASER.intNoise   = intNoise;
LASER.RIN_dB     = RIN_dB;

end
