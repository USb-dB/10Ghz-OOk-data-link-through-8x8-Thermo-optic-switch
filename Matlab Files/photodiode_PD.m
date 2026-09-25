function [io, PD] = photodiode_PD(s, Fs, PD)
% ---------------------------------------------------------
% Photodiode model (IM/DD)
%
% io(t) = R|s(t)|^2 + is(t) + it(t)
%
% INPUTS:
%   s   : received optical field (complex) [W^0.5]
%   Fs  : sampling frequency [Hz]
%   PD  : structure with photodiode parameters
%
% OUTPUTS:
%   io  : photodiode output current [A]
%   PD  : updated structure with noise statistics
% ---------------------------------------------------------

%% Default photodiode parameters
if ~isfield(PD,'R')
    PD.R = 0.8;               % Responsivity [A/W]
end
if ~isfield(PD,'BW')
    PD.BW = Fs/2;             % Electrical bandwidth [Hz]
end
if ~isfield(PD,'T')
    PD.T = 300;               % Temperature [K]
end
if ~isfield(PD,'RL')
    PD.RL = 50;               % Load resistance [Ohm]
end
if ~isfield(PD,'q')
    PD.q = 1.602e-19;         % Electron charge [C]
end
if ~isfield(PD,'kB')
    PD.kB = 1.38e-23;         % Boltzmann constant [J/K]
end

%% Optical power
P = abs(s).^2;                % Optical power [W]

%% Signal current
i_sig = PD.R * P;             % Signal current [A]

%% Shot noise
shot_var = 2 * PD.q * mean(i_sig) * PD.BW;
i_shot = sqrt(shot_var) * randn(size(i_sig));

%% Thermal noise
thermal_var = 4 * PD.kB * PD.T * PD.BW / PD.RL;
i_thermal = sqrt(thermal_var) * randn(size(i_sig));

%% Total photodiode output current
io = i_sig + i_shot + i_thermal;

%% Store noise statistics
PD.shotNoiseVar    = shot_var;
PD.thermalNoiseVar = thermal_var;
PD.meanCurrent     = mean(i_sig);

end

