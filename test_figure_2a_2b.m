%% This test script produces figure 2a and figure 2b in the paper

%%% addpath
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%% bSSFP signal evolution

%%% Case#1: single pool bSSFP
%%% Sequences (from scanner)
TR = 2.76;
alpha = 35;
npulse = 1000; 
phi = RF_phase_cycle(npulse,'balanced');
FA = d2r(alpha)*ones(npulse,1);

%%% Relaxation parameters: single pool (from Robson MRM paper)
T1 = 1175;
T2 = 54.4;

%%% Linear ramp excitation by D.G. Nishimura, Analysis and Reduction of the
%%% Transient Response in SSFP Imaging, Proc. Intl. Soc. Mag. Reson. Med (2000)
FA(1:10) = linspace(0, d2r(alpha), 10);

%%% Initialise magnetization for the first bSSFP
z0 = 1;
[s, Fn, Zn] = EPG_GRE_Chong(FA, phi, TR, T1, T2, 'kmax',...
    inf, 'zinit', z0);

%%% Save signal
%  Need iFFT to get off-resonance profile
mxys = ifftshift(size(Fn,1)*(ifft(ifftshift(Fn,1),[],1)),1);
phiTR= linspace(-pi, pi,size(mxys,1));

%% bSSFP signal evolution
%%% case #2: two pool bSSFP

%%% Relaxation parameters: MT
T1_MT = [1175 1175];
f = 0.07;    % Robson et al. MRM 2013
ka = 4.1e-3;  % msec-1
T2_MT = T2;    

%%% 
M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;

% RF saturation factor for MT, assuming T2r = 8.5us as described in Robson
% et al., MRM 2013, take the mean from -1kHz to 1kHz. 
G = 15.1;         % us   

heartRate = 1000; % ms per beat
kcenter = 50;   % 10 startup + 40 echeos

b1 = 1.42;      % uT, B1 rms read from scanner
gam = 267.5221*1e-3; %< rad /ms /uT
% trf = d2r(alpha)/(gam*b1);% ms
trf = 0.4; % ms
b1sqrdtau = b1^2*trf;
B1SqrdTau = b1sqrdtau*ones(npulse,1);

% setup linear startup echoes
B1SqrdTau(1:10) = B1SqrdTau(1:10).*linspace(0, 1, 10)';

z0 = [(1-f) f];
[s_mt, Fn_mt, Zn_mt] = EPGX_GRE_MT_HYTOM(FA,phi,B1SqrdTau,...
                TR,T1_MT,T2_MT,f,ka,G,'kmax',inf, 'zinit', z0);
            
mxys_mt = ifftshift(size(Fn_mt,1)*(ifft(ifftshift(Fn_mt,1),[],1)),1);

%% bSSFP signal evolution
%%% case #3:two pool + MT-prep bSSFP
pseudo_prep = struct('flip', 0, 't_delay', 0, 'B1SqrdTau', 0, 'MTC', true);
z0 = [(1-f) f];
[s_mtprep, Fn_mtprep, Zn_mtprep] = EPGX_GRE_MT_HYTOM(FA, phi, B1SqrdTau, TR, T1_MT,...
    T2_MT,f,ka,G,'kmax',inf, 'zinit', z0, 'prep', pseudo_prep);

mxys_mtprep = ifftshift(size(Fn_mtprep,1)*(ifft(ifftshift(Fn_mtprep,1),[],1)),1);

%%  plot bSSFP to steady-state
figure
% [~,idx0] = min(abs(phiTR - pi));
[~, idx_min] = min(std(abs(mxys(:, 11:50)),0,2));
plot(11:npulse, abs(mxys(idx_min,11:end)),'r-', 'linewidth',1.5); hold on
plot(11:npulse, abs(mxys_mt(idx_min,11:end)),'g-','linewidth',1.5);
plot(11:npulse, abs(mxys_mtprep(idx_min,11:end)),'b-','linewidth',1.5);
legend({'Single pool', 'Two pool', 'Two pool with MT-prep'}, 'Location','best')
legend boxoff
ylabel('Signal')
xlabel('TR number')
set(gca, 'Fontsize',8, 'xlim', [0 1000], 'ylim', [0.05 0.35], ...
    'Ytick', 0.05:0.1:0.35, 'Xtick', 0:250:1000)

% add insert
gg = get(gcf, 'children');
axes(gg(2))
ax = axes('Position',[0.5 0.4 0.2 0.15]);
box on
plot(900:npulse, abs(mxys(idx_min,900:end)), 'r-', 'linewidth', 1.5); hold on
plot(900:npulse, abs(mxys_mt(idx_min,900:end)), 'g-', 'linewidth', 1.5);
plot(900:npulse, abs(mxys_mtprep(idx_min,900:end)), 'b-', 'linewidth',1.5);
% ylim([0.07 0.11])
set(ax, 'xticklabel',[], 'yticklabel', [])
annotation('arrow',[0.7 0.85],[0.4 0.3]);

% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3.46 3];
% print('bSSFP_to_ss','-dtiff','-r900')

%% Plot bSSFP off-resonance profile
figure
y = abs(mxys(:,89)); % take kspace center
y_mt = abs(mxys_mt(:,89));
y_mtprep = abs(mxys_mtprep(:,89));
plot(phiTR, y, 'r-', 'linewidth',1.5); hold on
plot(phiTR, y_mt, 'g-', 'linewidth',1.5); 
plot(phiTR, y_mtprep, 'b-', 'linewidth',1.5)
xlim([-pi pi])
ylabel('Signal')
xlabel('\phi, rad')
legend({'Single pool', 'Two pool', 'Two pool with MT-prep'}, 'Location','best')
legend boxoff
set(gca, 'Fontsize',8, 'xlim', [-pi pi], 'ylim', [0 0.25], ...
    'Ytick', 0:0.1:0.25, 'Xtick', -3:1:3)

% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3.46 3];
% print('bSSFP_profile','-dtiff','-r900')