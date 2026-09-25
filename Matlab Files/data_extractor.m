clc; clear; close all;

%% ===================== LOAD RAW DATA =====================
data_te0 = readtable( ...
    'TE0.dat', ...
    'Delimiter','\t', ...
    'FileType','text', ...
    'VariableNamingRule','preserve' ...
);

% Extract columns
wav_nm = data_te0.("wavelength(1-1)");
TX     = data_te0.("tx(1-1)");
phi    = data_te0.("Angle(1-1)");

%% ===================== CLEAN DATA ========================
% Remove NaN / Inf / empty rows
valid = isfinite(wav_nm) & isfinite(TX) & isfinite(phi);

wav_nm = wav_nm(valid);
TX     = TX(valid);
phi    = phi(valid);

%% ===================== PHYSICAL CONSTANTS ================
c = 3e8;
lambda0 = 1550e-9;
omega0  = 2*pi*c/lambda0;

%% ===================== CONVERT TO FREQUENCY ==============
lambda = wav_nm * 1e-9;
omega  = 2*pi*c ./ lambda;

% Sort by increasing omega
[omega, idx] = sort(omega);
TX  = TX(idx);
phi = phi(idx);

% Remove duplicate frequency points (MANDATORY)
[omega, uniqIdx] = unique(omega, 'stable');
TX  = TX(uniqIdx);
phi = phi(uniqIdx);

%% ===================== COMPLEX TRANSFER FUNCTION =========
Omega = omega - omega0;          % Detuning (rad/s)
H     = TX .* exp(1j * phi);     % Complex response

%% ===================== FINAL CONSISTENCY CHECK ===========
assert(all(isfinite(Omega)), 'Omega contains NaN/Inf');
assert(all(isfinite(H)),     'H contains NaN/Inf');
assert(issorted(Omega),      'Omega must be sorted');

fprintf('TE0 data cleaned: %d points retained\n', length(Omega));

%% ===================== EXPORT DATA =======================

% ---- 1️⃣ MATLAB-native (recommended) ----
save('TE0_clean.mat','Omega','H','TX','phi','wav_nm');

% ---- 2️⃣ Human-readable ASCII ----
TE0_export = table( ...
    wav_nm, ...
    Omega, ...
    TX, ...
    phi, ...
    real(H), ...
    imag(H), ...
    'VariableNames', ...
    {'lambda_nm','Omega_rad_per_s','TX','phi_rad','H_real','H_imag'} ...
);

writetable(TE0_export,'TE0_clean.txt','Delimiter','\t');
writetable(TE0_export,'TE0_clean.csv');

%% ===================== QUICK SANITY PLOTS ================
figure('Color','w');

subplot(3,1,1)
plot(wav_nm, TX,'LineWidth',1.2)
xlabel('Wavelength (nm)')
ylabel('|H|')
grid on

subplot(3,1,2)
plot(wav_nm, phi,'LineWidth',1.2)
xlabel('Wavelength (nm)')
ylabel('Phase (rad)')
grid on

subplot(3,1,3)
plot(Omega/2/pi/1e9, unwrap(angle(H)),'LineWidth',1.2)
xlabel('Detuning (GHz)')
ylabel('Unwrapped Phase (rad)')
grid on
