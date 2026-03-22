% mprage_hcp_260321.m
% Generate 3D MPRAGE .seq files matching the HCP T1w protocol (0.7 mm iso).
%
% Outputs:
%   mprage_hcp.seq        - fully sampled
%   mprage_hcp_grappa.seq - GRAPPA R=2, 32 ACS lines
%
% Usage (after sourcing env.sh):
%   matlab -batch "addpath(pwd,'extern/pulseq/matlab'); mprage_hcp_260321"

%% System limits (conservative; works on Prisma-class hardware)
sys = mr.opts('MaxGrad', 30, 'GradUnit', 'mT/m', ...
    'MaxSlew', 130, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% HCP T1w MPRAGE parameters
alpha         = 8;           % flip angle [deg]
rfLen         = 100e-6;      % hard (block) pulse duration [s]
TI            = 1.0;         % inversion time [s]
TRout         = 2.4;         % outer TR (inversion interval) [s]
ro_os         = 1;           % readout oversampling
ro_spoil      = 3;           % RO spoiling [multiples of kmax]
rfSpoilingInc = 84;          % RF spoiling phase increment [deg]

% Readout duration from BW ~ 210 Hz/pixel.
% Dwell must be a multiple of 500 ns so that dwell*N_adc lands on
% the 10 us block duration raster (with N_ro=320).
BW_per_pixel = 210;          % [Hz/pixel]
N_ro         = 320;          % readout matrix (base resolution)
dwellRaster  = 500e-9;       % lcm(adcRaster=100ns, blockRaster/N_ro)
dwell  = round(1/(BW_per_pixel * N_ro) / dwellRaster) * dwellRaster;  % 15.0 us
ro_dur = dwell * N_ro * ro_os;  % 4.800 ms  (BW ~ 208 Hz/pixel)

% Sagittal 0.7 mm isotropic
% Encoding axes: x = L/R (partition), y = A/P (phase), z = S/I (readout)
phaseOS       = 0.10;           % 10% phase oversampling (avoids aliasing at FOV edges)
N_phase_base  = 320;            % base phase-encoding matrix
N_phase       = round(N_phase_base * (1 + phaseOS));  % 352 with oversampling
fov_phase     = 224e-3 * (1 + phaseOS);               % FOV extended to match

fov = [179.2e-3 fov_phase 224e-3];  % [x y z] FOV [m]
N   = [256  N_phase 320];           % [x y z] matrix

ax.d1 = 'z';  % readout        (S/I)
ax.d2 = 'x';  % partition enc  (L/R, inner loop)
ax.d3 = 'y';  % phase enc      (A/P, outer loop — GRAPPA applied here)
ax.n1 = strfind('xyz', ax.d1);  % 3
ax.n2 = strfind('xyz', ax.d2);  % 1
ax.n3 = strfind('xyz', ax.d3);  % 2

% GRAPPA
accelFactorPE = 2;
ACSnum        = 32;

%% Pulse and gradient construction
% Wrapped in a local function-like block so objects are rebuilt fresh
% for each sequence variant (avoids stale .id fields from registerXxxEvent).

deltak = 1 ./ fov;

% These values are invariant across variants, compute once for reporting
[rf0, rf180_0, adc0, gro10, groSp0, groPre0, gpe10, gpe20, gslSp0, TRinner0, TIdelay0, TRoutDelay0] = ...
    build_mprage_events(sys, alpha, rfLen, ro_dur, ro_os, ro_spoil, TI, TRout, ax, N, fov, deltak);

pe1Steps = ((0:N(ax.n2)-1) - N(ax.n2)/2) / N(ax.n2) * 2;
pe2Steps = ((0:N(ax.n3)-1) - N(ax.n3)/2) / N(ax.n3) * 2;

TE  = mr.calcDuration(rf0) - rf0.shape_dur/2 - rf0.delay ...
    + mr.calcDuration(gro10) - mr.calcDuration(adc0)/2;
ESP = mr.calcDuration(rf0) + mr.calcDuration(gro10);

fprintf('TRinner = %.3f ms, TE = %.3f ms, ESP = %.3f ms\n', TRinner0*1e3, TE*1e3, ESP*1e3);
fprintf('TIdelay = %.1f ms, TRoutDelay = %.1f ms\n', TIdelay0*1e3, TRoutDelay0*1e3);
fprintf('BW = %.1f Hz/pixel (dwell = %.1f ns)\n', 1/(adc0.dwell*N(ax.n1)), adc0.dwell*1e9);

% Phase resolution ratio (for GRAPPA metadata)
phaseResolution = fov(ax.n3)/N(ax.n3) / (fov(ax.n2)/N(ax.n2));

%% GRAPPA sampling pattern
centerLineIdx = floor(N(ax.n3)/2) + 1;
PEsamp_u = [];
cnt = 1;
for ii = 1:N(ax.n3)
    if mod(ii - centerLineIdx, accelFactorPE) == 0
        PEsamp_u(cnt) = ii;
        cnt = cnt + 1;
    end
end
PEsamp_ACS = (centerLineIdx - ACSnum/2) : (centerLineIdx + floor(ACSnum-1)/2);
PEsamp_grappa     = union(PEsamp_u, PEsamp_ACS);
PEsamp_grappa_INC = diff([PEsamp_grappa, PEsamp_grappa(end)]);

fprintf('GRAPPA: %d / %d PE lines sampled (R=%d, %d ACS)\n', ...
    length(PEsamp_grappa), N(ax.n3), accelFactorPE, ACSnum);

%% Generate both sequence variants
for use_grappa = [false, true]
    if use_grappa
        seqName  = 'mprage_hcp_grappa';
        seqLabel = 'mp_hcp_gt';
        PEsamp     = PEsamp_grappa;
        PEsamp_INC = PEsamp_grappa_INC;
    else
        seqName  = 'mprage_hcp';
        seqLabel = 'mprage_hcp';
        PEsamp     = 1:N(ax.n3);
        PEsamp_INC = ones(1, N(ax.n3));
    end
    nPEsamp = length(PEsamp);
    fprintf('\n=== %s (%d outer PE steps) ===\n', seqName, nPEsamp);

    % Rebuild all pulse/gradient objects fresh (avoids stale .id fields)
    [rf, rf180, adc, gro1, groSp, groPre, gpe1, gpe2, gslSp, TRinner, TIdelay, TRoutDelay] = ...
        build_mprage_events(sys, alpha, rfLen, ro_dur, ro_os, ro_spoil, TI, TRout, ax, N, fov, deltak);

    seq = mr.Sequence(sys);

    % Labels
    lblIncPar   = mr.makeLabel('INC', 'PAR', 1);
    lblResetPar = mr.makeLabel('SET', 'PAR', 0);

    % Pre-register invariant events
    gslSp.id = seq.registerGradEvent(gslSp);
    groSp.id = seq.registerGradEvent(groSp);
    gro1.id  = seq.registerGradEvent(gro1);
    [~, rf.shapeIDs]         = seq.registerRfEvent(rf);
    [rf180.id, rf180.shapeIDs] = seq.registerRfEvent(rf180);
    lblIncPar.id = seq.registerLabelEvent(lblIncPar);

    if use_grappa
        lblSetRef    = mr.makeLabel('SET', 'REF', true);
        lblSetIma    = mr.makeLabel('SET', 'IMA', true);
        lblResetRef  = mr.makeLabel('SET', 'REF', false);
        lblResetIma  = mr.makeLabel('SET', 'IMA', false);
        lblSetRef.id   = seq.registerLabelEvent(lblSetRef);
        lblSetIma.id   = seq.registerLabelEvent(lblSetIma);
        lblResetRef.id = seq.registerLabelEvent(lblResetRef);
        lblResetIma.id = seq.registerLabelEvent(lblResetIma);

        % Noise prescan
        seq.addBlock(adc, mr.makeLabel('SET','LIN',0), mr.makeLabel('SET','NOISE',true), ...
            lblResetIma, lblResetRef);
        seq.addBlock(mr.makeLabel('SET','NOISE',false));
        % Initial line counter
        seq.addBlock(mr.makeLabel('SET','LIN', PEsamp(1)-1));
    end

    tic;
    for pe_idx = 1:nPEsamp
        % --- GRAPPA REF/IMA flags ---
        if use_grappa
            is_acs = ismember(PEsamp(pe_idx), PEsamp_ACS);
            is_usamp = ismember(PEsamp(pe_idx), PEsamp_u);
            if is_acs && is_usamp
                seq.addBlock(lblSetIma, lblSetRef);
            elseif is_acs
                seq.addBlock(lblResetIma, lblSetRef);
            else
                seq.addBlock(lblResetIma, lblResetRef);
            end
        end

        % --- Inversion pulse + TI delay ---
        seq.addBlock(rf180);
        seq.addBlock(mr.makeDelay(TIdelay), gslSp);

        % Pre-register PE2 gradients for this outer step
        gpe2e = mr.scaleGrad(gpe2, pe2Steps(PEsamp(pe_idx)));
        gpe2e.id = seq.registerGradEvent(gpe2e);
        gpe2r = mr.scaleGrad(gpe2, -pe2Steps(PEsamp(pe_idx)));
        gpe2r.id = seq.registerGradEvent(gpe2r);

        rf_phase = 0;
        rf_inc   = 0;

        % --- Inner loop: partition encoding + GRE readout ---
        for ii = 1:N(ax.n2)
            rf.phaseOffset  = rf_phase/180*pi;
            adc.phaseOffset = rf_phase/180*pi;
            rf_inc   = mod(rf_inc   + rfSpoilingInc, 360.0);
            rf_phase = mod(rf_phase + rf_inc,         360.0);

            if ii == 1
                seq.addBlock(rf);
            else
                seq.addBlock(rf, groSp, mr.scaleGrad(gpe1, -pe1Steps(ii-1)), gpe2r, lblIncPar);
            end
            seq.addBlock(adc, gro1, mr.scaleGrad(gpe1, pe1Steps(ii)), gpe2e);
        end

        % --- End of inversion block: rewind last PE + TR dead-time ---
        seq.addBlock(groSp, mr.scaleGrad(gpe1, -pe1Steps(N(ax.n2))), gpe2r, ...
            mr.makeDelay(TRoutDelay), ...
            mr.makeLabel('INC', 'LIN', PEsamp_INC(pe_idx)), lblResetPar);
    end
    fprintf('Block generation: %.1f s\n', toc);

    % --- Timing check ---
    [ok, error_report] = seq.checkTiming();
    if ok
        fprintf('Timing check passed\n');
    else
        fprintf('Timing check FAILED:\n');
        fprintf([error_report{:}]);
        fprintf('\n');
    end

    % --- Sequence metadata ---
    seq.setDefinition('FOV', fov);
    seq.setDefinition('Name', seqLabel);
    seq.setDefinition('OrientationMapping', 'SAG');
    seq.setDefinition('ReceiverGainHigh', 1);
    if use_grappa
        seq.setDefinition('ReadoutOversamplingFactor', ro_os);
        seq.setDefinition('kSpaceCenterLine', centerLineIdx - 1);
        seq.setDefinition('PhaseResolution', phaseResolution);
    end

    % --- Write .seq file ---
    fname = [seqName '.seq'];
    seq.write(fname);
    fprintf('Written: %s\n', fname);
end

fprintf('\nDone.\n');


function [rf, rf180, adc, gro1, groSp, groPre, gpe1, gpe2, gslSp, TRinner, TIdelay, TRoutDelay] = ...
        build_mprage_events(sys, alpha, rfLen, ro_dur, ro_os, ro_spoil, TI, TRout, ax, N, fov, deltak)
    rf = mr.makeBlockPulse(alpha*pi/180, sys, 'Duration', rfLen, 'use', 'excitation');
    rf180 = mr.makeAdiabaticPulse('hypsec', sys, 'Duration', 10.24e-3, 'dwell', 1e-5, ...
        'use', 'inversion');

    gro = mr.makeTrapezoid(ax.d1, ...
        'Amplitude', N(ax.n1)*deltak(ax.n1)/ro_dur, ...
        'FlatTime', ceil((ro_dur + sys.adcDeadTime)/sys.gradRasterTime)*sys.gradRasterTime, ...
        'system', sys);
    adc = mr.makeAdc(N(ax.n1)*ro_os, 'Duration', ro_dur, 'Delay', gro.riseTime, 'system', sys);

    groPre = mr.makeTrapezoid(ax.d1, ...
        'Area', -gro.amplitude*(adc.dwell*(adc.numSamples/2 + 0.5) + 0.5*gro.riseTime), ...
        'system', sys);

    gpe1 = mr.makeTrapezoid(ax.d2, 'Area', -deltak(ax.n2)*(N(ax.n2)/2), 'system', sys);
    gpe2 = mr.makeTrapezoid(ax.d3, 'Area', -deltak(ax.n3)*(N(ax.n3)/2), 'system', sys);
    gslSp = mr.makeTrapezoid(ax.d3, 'Area', max(deltak.*N)*4, 'Duration', 10e-3, 'system', sys);

    [gro1, groSp] = mr.splitGradientAt(gro, gro.riseTime + gro.flatTime);
    if ro_spoil > 0
        groSp = mr.makeExtendedTrapezoidArea(gro.channel, gro.amplitude, 0, ...
            deltak(ax.n1)/2*N(ax.n1)*ro_spoil, sys);
    end

    rf.delay = mr.calcDuration(groSp, gpe1, gpe2);
    [groPre, ~, ~] = mr.align('right', groPre, gpe1, gpe2);
    gro1.delay = mr.calcDuration(groPre);
    adc.delay  = gro1.delay + gro.riseTime;
    gro1 = mr.addGradients({gro1, groPre}, 'system', sys);
    TRinner = mr.calcDuration(rf) + mr.calcDuration(gro1);

    pe1Steps = ((0:N(ax.n2)-1) - N(ax.n2)/2) / N(ax.n2) * 2;
    TIdelay = round((TI ...
        - (find(pe1Steps==0) - 1)*TRinner ...
        - (mr.calcDuration(rf180) - mr.calcRfCenter(rf180) - rf180.delay) ...
        - rf.delay - mr.calcRfCenter(rf)) / sys.blockDurationRaster) * sys.blockDurationRaster;
    TRoutDelay = TRout - TRinner*N(ax.n2) - TIdelay - mr.calcDuration(rf180);
end
