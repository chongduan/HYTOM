%% This test script generats figure 2c in the paper

%%% addpath
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%% MOLLI T1 measurement

%%% Case #1: single pool MOLLI
%%% Sequences (from scanner)
TR = 2.76;
alpha = 35;

%%% Relaxation parameters: single pool (from Robson MRM paper)
T1 = 1175;
T2 = 54.4;

heartRate = 1000; % ms per beat
kcenter = 40; % 10 startup + 30 echeos

npulse = 79; % 
phi = RF_phase_cycle(npulse,'balanced');

FA = d2r(alpha)*ones(npulse,1);
FA(1:10) = linspace(0, d2r(alpha), 10);

% Initiate signal cell (two MOLLI blocks)
Sig = cell(2,1);

% Delay time for each block
Inv_delay = [160, 212];

for MOLLI_block = 1:2
    switch MOLLI_block
        case 1
            n_bSSFP = 5;
            % set inversion pulse
            prep = struct('flip', pi, 't_delay', Inv_delay(1));
        case 2
            n_bSSFP = 3;
            % set inversion pulse
            prep = struct('flip', pi, 't_delay', Inv_delay(2));
    end
    
    for i = 1:n_bSSFP
        % single pool
        if i == 1
            %%% Initialise magnetization for the first bSSFP
            z0 = 1;
            [s, Fn, Zn] = EPG_GRE_Chong(FA, phi, TR, T1, T2, 'kmax',...
                inf, 'zinit', z0, 'prep', prep);
        else
            [s, Fn, Zn] = EPG_GRE_Chong(FA, phi, TR, T1, T2, 'kmax',...
                inf, 'zinit', z0);
        end
        
        %%% Save signal
        mxys = size(Fn,1)*ifftshift(ifft(ifftshift(Fn,1),[],1),1);
        
        %%% Select dephasing paramter
        phiTR= linspace(-pi,pi,size(mxys,1));
        [~,idx0] = min(abs(phiTR));
        
%         [~,idx50] = min(abs(phiTR-pi/2));
        
        %     % plot
        %     figure; plot(abs(mxys(idx0,:)))
        
        Sig{MOLLI_block}(i, :) = mxys(idx0, kcenter);

        mxys_signal{MOLLI_block}(i,:,:) = mxys;
        %%% Apply free precession to z0, assuming all transverse components
        %%% completely dephased
        
        % First take Z0 at end of TSE shot
        z0 = squeeze(Zn(1,end));
        
        % Now work out delay to add after each bSSFP
        Tdelay = heartRate - npulse*TR;
        Xi = exp(-Tdelay/T1);
        Zoff = (1-exp(-Tdelay/T1));
        
        % Now evolve it by amount due to recovery period
        z0 = Xi*z0 + Zoff;
    end
end

Signal = [Sig{1}; Sig{2}];
            
% Phsse signal based on the last one
Signal_phased_s = Signal * conj(Signal(end))/abs(Signal(end));

% plot signal 
TI_s = [Inv_delay(1) + heartRate*(0:4) Inv_delay(2) + heartRate*(0:2)] + TR*kcenter;


%% Case #2: two pool MOLLI

%%% Relaxation parameters: MT
T1_MT = [1175 1175];
f = 0.07;    % Robson et al. MRM 2013
ka = 4.1e-3;  % msec-1
T2_MT = 54.4;    

%%% 
M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;

R1a = 1/T1_MT(1);
R1b = 1/T1_MT(2);
R2a = 1/T2_MT;

% RF saturation factor for MT, assuming T2r = 8.5us as described in Robson
% et al., MRM 2013, take the mean from -1kHz to 1kHz. 
G = 15.1;         % us   

heartRate = 1000; % ms per beat
% kcenter = 50;   % 10 startup + 31 echeos

npulse = 79;    % 79 echoes + 10 startup
phi = RF_phase_cycle(npulse,'balanced');

b1 = 1.42;      % uT, B1 rms read from scanner
gam = 267.5221*1e-3; %< rad /ms /uT
% trf = d2r(alpha)/(gam*b1);% ms
trf = 0.44; % ms
b1sqrdtau = b1^2*trf;
B1SqrdTau = b1sqrdtau*ones(npulse,1);

% setup linear startup echoes
FA = d2r(alpha)*ones(npulse,1);
FA(1:10) = linspace(0, d2r(alpha), 10);
B1SqrdTau(1:10) = B1SqrdTau(1:10).*linspace(0, 1, 10)';

% Initiate signal cell (two MOLLI blocks)
Sig = cell(2,1);

% Delay time for each block
Inv_delay = [160, 212];

for MOLLI_block = 1:2
    switch MOLLI_block
        case 1
            n_bSSFP = 5;
            % set inversion pulse
            prep = struct('flip', pi, 't_delay', Inv_delay(1), 'B1SqrdTau', b1sqrdtau);
        case 2
            n_bSSFP = 3;
            % set inversion pulse
            prep = struct('flip', pi, 't_delay', Inv_delay(2), 'B1SqrdTau', b1sqrdtau);
    end
    
    for i = 1:n_bSSFP
        if i == 1
            %%% Initialise magnetization for the first bSSFP
            z0 = [(1-f) f];
            [s, Fn, Zn] = EPGX_GRE_MT_HYTOM(FA,phi,B1SqrdTau,...
                TR,T1_MT,T2_MT,f,ka,G,'kmax',inf, 'zinit', z0, 'prep', prep);
        else
            [s, Fn, Zn] = EPGX_GRE_MT_HYTOM(FA,phi,B1SqrdTau,...
                TR,T1_MT,T2_MT,f,ka,G,'kmax',inf, 'zinit', z0);
        end
        
        %%% Save signal
        mxys = size(Fn,1)*ifftshift(ifft(ifftshift(Fn,1),[],1),1);
        
        %%% Select dephasing paramter
        phiTR= linspace(-pi,pi,size(mxys,1));
        [~,idx0] = min(abs(phiTR));
        
        %     % plot
        %     figure; plot(abs(mxys(idx0,:)))
        
        Sig{MOLLI_block}(i, :) = mxys(idx0, kcenter);
        
        %%% Apply free precession to z0, assuming all transverse components
        %%% completely dephased
        
        % First take Z0 at end of bSSFP
        z0 = squeeze(Zn(1,end,:));
        
        % Now work out delay to add after each bSSFP
        Tdelay = heartRate - npulse*TR;
        
        % Evolve Z0 between bSSFP readout
        L = [[-R1a-ka kb];[ka -R1b-kb]];
        C = [R1a*(1-f) R1b*f]';
        Xi = expm(L*Tdelay);
        I=eye(2);
        Zoff = (Xi - I)*inv(L)*C;
        
        %%% Now evolve it by amount due to recovery period
        z0 = Xi*z0 + Zoff;
    end
end

Signal = [Sig{1}; Sig{2}];
            
% Phsse signal based on the last one
Signal_phased = Signal * conj(Signal(end))/abs(Signal(end));

% plot signal 
TI_t = [Inv_delay(1) + heartRate*(0:4) Inv_delay(2) + heartRate*(0:2)] + TR*kcenter;

%% Case #3: two pool + MT-prep (i.e., HYTOM)

%%% Relaxation parameters: MT
T1_MT = [1175 1175];
f = 0.07;    % Robson et al. MRM 2013
ka = 4.1e-3;  % msec-1
T2_MT = 54.4;  

%%% 
M0b = f;
M0a = (1-f);
kb = ka * M0a/M0b;

R1a = 1/T1_MT(1);
R1b = 1/T1_MT(2);
R2a = 1/T2_MT;

%%% RF saturation factor for MT
G = 12;% us   

heartRate = 1000; % ms per beat
% kcenter = 50; % 10 startup + 40 echeos

npulse = 79; % 79 echoes + 10 startup
phi = RF_phase_cycle(npulse,'balanced');

b1 = 5.1; % uT
% gam = 267.5221*1e-3; %< rad /ms /uT
% trf = d2r(alpha)/(gam*b1);% ms
trf = 0.44;
b1sqrdtau = b1^2*trf;
B1SqrdTau = b1sqrdtau*ones(npulse,1);

% setup linear startup echoes
FA = d2r(alpha)*ones(npulse,1);
FA(1:10) = linspace(0, d2r(alpha), 10);
B1SqrdTau(1:10) = B1SqrdTau(1:10).*linspace(0, 1, 10)';

% Initiate signal cell (two MOLLI blocks)
Sig = cell(2,1);

% Delay time for each block
Inv_delay = [0, 52];

for MOLLI_block = 1:2
    switch MOLLI_block
        case 1
            n_bSSFP = 5;
            % set inversion pulse, ignore the MT effect of the inversion
            prep = struct('flip', pi, 't_delay', Inv_delay(1), 'B1SqrdTau', b1sqrdtau, 'MTC', true);
        case 2
            n_bSSFP = 3;
            % set inversion pulse
            prep = struct('flip', pi, 't_delay', Inv_delay(2), 'B1SqrdTau', b1sqrdtau, 'MTC', true);
    end
    
    for i = 1:n_bSSFP
        if i == 1
            %%% Initialise magnetization for the first bSSFP
            z0 = [(1-f) f];
            [s, Fn, Zn] = EPGX_GRE_MT_HYTOM(FA,phi,B1SqrdTau,...
                TR,T1_MT,T2_MT,f,ka,G,'kmax',inf, 'zinit', z0, 'prep', prep);
        else
            % no inversion pulse, only MT-prep 
            pseudo_prep = struct('flip', 0, 't_delay', 0, 'B1SqrdTau', ...
                0, 'MTC', true);
            [s, Fn, Zn] = EPGX_GRE_MT_HYTOM(FA,phi,B1SqrdTau,TR,T1_MT,...
                T2_MT,f,ka,G,'kmax',inf, 'zinit', z0, 'prep', pseudo_prep);
        end
        %%% Save signal
        mxys = size(Fn,1)*ifftshift(ifft(ifftshift(Fn,1),[],1),1);
        
        %%% Select dephasing paramter
        phiTR= linspace(-pi,pi,size(mxys,1));
        [~,idx0] = min(abs(phiTR));
        
        %     % plot
        %     figure; plot(abs(mxys(idx0,:)))
        
        Sig{MOLLI_block}(i, :) = mxys(idx0, kcenter);
        
        %%% Apply free precession to z0, assuming all transverse components
        %%% completely dephased
        
        % First take Z0 at end of bSSFP
        z0 = squeeze(Zn(1,end,:));
        
        % Now work out delay to add after each bSSFP
        Tdelay = heartRate - npulse*TR - 160; % 160ms MTC-prep time
        
        % Evolve Z0 between bSSFP readout
        L = [[-R1a-ka kb];[ka -R1b-kb]];
        C = [R1a*(1-f) R1b*f]';
        Xi = expm(L*Tdelay);
        I=eye(2);
        Zoff = (Xi - I)*inv(L)*C;
        
        %%% Now evolve it by amount due to recovery period
        z0 = Xi*z0 + Zoff;
    end
end

Signal = [Sig{1}; Sig{2}];
            
% Phsse signal based on the last one
Signal_phased_MTprep = Signal * conj(Signal(end))/abs(Signal(end));

% % plot signal 
% TI = [0 + 160 + heartRate*(0:4) 100 + 160 + heartRate*(0:2)] + TR*kcenter;
TI = [Inv_delay(1) + 160 + TR*kcenter + heartRate*(0:4) ...
    Inv_delay(2) + 160 + TR*kcenter + heartRate*(0:2)];

%% Plot the relaxation curves

figure
[x1_s, ~] = fitFunctionT1IR_3p(real(Signal_phased_s),TI_s);
xd = 0:10:5000;
y_fit = x1_s(1) - x1_s(2)*exp(-xd(:)/x1_s(3));
h2 = plot(xd, y_fit, 'r-', 'linewidth',1.5); hold on
h1 = plot(TI_s, real(Signal_phased_s),'ro', 'linewidth',1.5, 'MarkerSize', 4);
% x1_s(3)*(x1_s(2)/x1_s(1) - 1)

[x1, ~] = fitFunctionT1IR_3p(real(Signal_phased),TI);
xd = 0:10:5000;
y_fit = x1(1) - x1(2)*exp(-xd(:)/x1(3));
h4 = plot(xd, y_fit, 'g-', 'linewidth',1.5);
h3 = plot(TI, real(Signal_phased),'go','linewidth',1.5, 'MarkerSize', 4); 
% x1(3)*(x1(2)/x1(1) - 1)

[x1, ~] = fitFunctionT1IR_3p(real(Signal_phased_MTprep),TI);
xd = 0:10:5000;
y_fit = x1(1) - x1(2)*exp(-xd(:)/x1(3));
h6 = plot(xd, y_fit, 'b-', 'linewidth',1.5);
h5 = plot(TI, real(Signal_phased_MTprep),'bo','linewidth',1.5, 'MarkerSize', 4); 
% x1(3)*(x1(2)/x1(1) - 1)

legend([h2, h4, h6], {'Single pool', 'Two pool', 'Two with MT-prep'}, 'Location','east')
legend boxoff

xlabel('Inversion time, ms')
ylabel('Signal')

set(gca, 'Fontsize',8, 'xlim', [0 5000], 'ylim', [-0.3 0.2], ...
    'Ytick', -0.3:0.1:0.2, 'Xtick', 0:1000:5000)

% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 3.46 3];
% print('Longitudinal Relax','-dtiff','-r900')
