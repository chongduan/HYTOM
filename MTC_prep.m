function [ zout ] = MTC_prep( z0, T1x, T2a, f, ka, G)
%   [ zout ] = MTC_prep( z0, T1x, T2a, f, ka, G)
%   Simulate the effect of long MTC-prep pulse on Mz

%   Using the EPG-X framework to simulate the effect of MTC-prep pulses on
%   Mz of free water pool and bound-water pool

%   z0......    Initial z0[z0f, z0b]
%   T1x.....    T1x[T1_f, T1_b]
%   T2a.....    T2a, i.e., T2 of the free water only
%   f.......    fraction of bound water
%   ka......    exchange rate constant
%   G.......    MT bound water SL linewidth at the RF frequence

%   Chong Duan, 2018/02/13


% Check number of inputs.
if nargin > 6
    error('TooManyInputs');
end

TR = 1;     % Divide the MT-prep pulse into 1 ms unit (could be smaller) 
% angle = 800;
TauRF = 20;
Repeat = 8;

%%% Sequences
npulse = round(TauRF/TR)*Repeat; % 20 msec pulse, repeat 8 times 
phi = RF_phase_cycle(npulse,'balanced');
FA = zeros(npulse,1);   

% Calculate the MTC-prep power
% angle = d2r(angle);
% gam = 267.5221 * 1e-3;       % rad /ms /uT
% b1 = angle/(gam * TauRF);   % uT

% RF`mtc:[ 0 ]:B1 read from scanner
b1 = 15.2;

% Divide the 20ms MTC pulse into TR modules, hard pulse approximation
b1sqrdtau = b1^2*TR;
B1SqrdTau = b1sqrdtau*ones(npulse,1);

[~, ~, Zn] = EPGX_GRE_MT_HYTOM(FA, phi, B1SqrdTau, TR, T1x, T2a, f, ka, G,...
    'kmax',inf, 'zinit', z0);

% Take Z0 at end of MTC-prep, assuming Mxy completely dephased with crusher
% gradient
zout = squeeze(Zn(1,end,:))';

end