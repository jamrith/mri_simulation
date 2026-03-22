% 260321_write_mprage_grappa.m
% Patched copy of extern/pulseq/matlab/demoSeq/writeMPRAGE_grappa.m
% Removes seq.plot()/figure calls for headless batch execution on HPC.
%
% Generates: mprage_grappa.seq (GRAPPA R=2, 32 ACS, sagittal, 192x240x256)
%
% Usage (after sourcing env.sh):
%   matlab -batch "addpath('extern/pulseq/matlab'); run('260321_write_mprage_grappa.m')"

addpath(fullfile(fileparts(mfilename('fullpath')), 'extern', 'pulseq', 'matlab'));

% set system limits
sys = mr.opts('MaxGrad', 24, 'GradUnit', 'mT/m', ...
    'MaxSlew', 100, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq = mr.Sequence(sys);
alpha = 7;
ro_dur = 5120e-6;
ro_os = 2;                        % readout oversampling
ro_spoil = 3;
TI = 1.1;
TRout = 2.5;
rfSpoilingInc = 84;
rfLen = 100e-6;
ax = struct;
% sagittal
fov = [192 240 256]*1e-3;
N = [192 240 256];
phaseResoluion = fov(3)/N(3) / (fov(2)/N(2));
ax.d1 = 'z';
ax.d2 = 'x';

ax.d3 = setdiff('xyz', [ax.d1 ax.d2]);
ax.n1 = strfind('xyz', ax.d1);
ax.n2 = strfind('xyz', ax.d2);
ax.n3 = strfind('xyz', ax.d3);

%% Pulses and gradients

rf = mr.makeBlockPulse(alpha*pi/180, sys, 'Duration', rfLen, 'use', 'excitation');
rf180 = mr.makeAdiabaticPulse('hypsec', sys, 'Duration', 10.24e-3, 'dwell', 1e-5,...
    'use', 'inversion');

deltak = 1./fov;
gro = mr.makeTrapezoid(ax.d1,'Amplitude',N(ax.n1)*deltak(ax.n1)/ro_dur,'FlatTime',ceil((ro_dur+sys.adcDeadTime)/sys.gradRasterTime)*sys.gradRasterTime,'system',sys);
adc = mr.makeAdc(N(ax.n1)*ro_os,'Duration',ro_dur,'Delay',gro.riseTime,'system',sys);
groPre = mr.makeTrapezoid(ax.d1,'Area',-gro.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gro.riseTime),'system',sys);
gpe1 = mr.makeTrapezoid(ax.d2,'Area',-deltak(ax.n2)*(N(ax.n2)/2),'system',sys);
gpe2 = mr.makeTrapezoid(ax.d3,'Area',-deltak(ax.n3)*(N(ax.n3)/2),'system',sys);
gslSp = mr.makeTrapezoid(ax.d3,'Area',max(deltak.*N)*4,'Duration',10e-3,'system',sys);
[gro1, groSp] = mr.splitGradientAt(gro,gro.riseTime+gro.flatTime);
if ro_spoil>0
    groSp=mr.makeExtendedTrapezoidArea(gro.channel,gro.amplitude,0,deltak(ax.n1)/2*N(ax.n1)*ro_spoil,sys);
end

rf.delay=mr.calcDuration(groSp, gpe1, gpe2);
[groPre,~,~]=mr.align('right',groPre,gpe1,gpe2);
gro1.delay=mr.calcDuration(groPre);
adc.delay=gro1.delay+gro.riseTime;
gro1=mr.addGradients({gro1,groPre},'system',sys);
TRinner=mr.calcDuration(rf)+mr.calcDuration(gro1);

pe1Steps=((0:N(ax.n2)-1)-N(ax.n2)/2)/N(ax.n2)*2;
pe2Steps=((0:N(ax.n3)-1)-N(ax.n3)/2)/N(ax.n3)*2;

TIdelay=round((TI-(find(pe1Steps==0)-1)*TRinner-(mr.calcDuration(rf180)-mr.calcRfCenter(rf180)-rf180.delay)-rf.delay-mr.calcRfCenter(rf))/sys.blockDurationRaster)*sys.blockDurationRaster;
TRoutDelay=TRout-TRinner*N(ax.n2)-TIdelay-mr.calcDuration(rf180);

TE = mr.calcDuration(rf) - rf.shape_dur/2 -rf.delay + mr.calcDuration(gro1) - mr.calcDuration(adc)/2;
ESP = mr.calcDuration(rf) + mr.calcDuration(gro1);

%% Labels

lblIncPar=mr.makeLabel('INC','PAR', 1);
lblResetPar=mr.makeLabel('SET','PAR', 0);

lblSetRefScan = mr.makeLabel('SET','REF', true);
lblSetRefAndImaScan = mr.makeLabel('SET','IMA', true);
lblResetRefScan = mr.makeLabel('SET','REF', false);
lblResetRefAndImaScan = mr.makeLabel('SET','IMA', false);

%% Pre-register

gslSp.id=seq.registerGradEvent(gslSp);
groSp.id=seq.registerGradEvent(groSp);
gro1.id=seq.registerGradEvent(gro1);
[~, rf.shapeIDs]=seq.registerRfEvent(rf);
[rf180.id, rf180.shapeIDs]=seq.registerRfEvent(rf180);

lblSetRefScan.id=seq.registerLabelEvent(lblSetRefScan);
lblSetRefAndImaScan.id=seq.registerLabelEvent(lblSetRefAndImaScan);
lblResetRefScan.id=seq.registerLabelEvent(lblResetRefScan);
lblResetRefAndImaScan.id=seq.registerLabelEvent(lblResetRefAndImaScan);

%% GRAPPA sampling pattern

nY = N(ax.n3);
accelFactorPE = 2;
ACSnum = 32;
centerLineIdx = floor(nY/2) + 1;
count = 1;
clear PEsamp_u;
for i = 1:nY
    if ( mod(i-centerLineIdx, accelFactorPE)==0 )
        PEsamp_u(count) = i;
        count = count + 1;
    end
end
minPATRefLineIdx = centerLineIdx - ACSnum/2;
maxPATRefLineIdx = centerLineIdx + floor(ACSnum-1)/2;
PEsamp_ACS = minPATRefLineIdx : maxPATRefLineIdx;
PEsamp = union(PEsamp_u, PEsamp_ACS);
nPEsamp = length(PEsamp);
PEsamp_INC = diff([PEsamp, PEsamp(end)]);

fprintf('GRAPPA: %d / %d PE lines (R=%d, %d ACS)\n', nPEsamp, nY, accelFactorPE, ACSnum);

% reverse gradient polarity to match Siemens product sequence
groSp = mr.scaleGrad(groSp, -1);
gro1 = mr.scaleGrad(gro1, -1);
gpe1.amplitude = -gpe1.amplitude;
gslSp.amplitude = -gslSp.amplitude;

%% Build sequence

tic;

% Noise prescan
seq.addBlock(adc, mr.makeLabel('SET', 'LIN', 0), mr.makeLabel('SET', 'NOISE', true), lblResetRefAndImaScan, lblResetRefScan);
seq.addBlock(mr.makeLabel('SET', 'NOISE', false));

% Initial line counter
seq.addBlock(mr.makeLabel('SET', 'LIN', PEsamp(1)-1));

for count=1:nPEsamp
    % PAT labels
    if ismember(PEsamp(count), PEsamp_ACS)
        if ismember(PEsamp(count), PEsamp_u)
            seq.addBlock(lblSetRefAndImaScan, lblSetRefScan);
        else
            seq.addBlock(lblResetRefAndImaScan, lblSetRefScan);
        end
    else
        seq.addBlock(lblResetRefAndImaScan, lblResetRefScan);
    end

    seq.addBlock(rf180);
    seq.addBlock(mr.makeDelay(TIdelay), gslSp);
    rf_phase = 0;
    rf_inc = 0;
    gpe2je = mr.scaleGrad(gpe2, pe2Steps(PEsamp(count)));
    gpe2je.id = seq.registerGradEvent(gpe2je);
    gpe2jr = mr.scaleGrad(gpe2, -pe2Steps(PEsamp(count)));
    gpe2jr.id = seq.registerGradEvent(gpe2jr);

    for i=1:N(ax.n2)
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        if (i==1)
            seq.addBlock(rf);
        else
           seq.addBlock(rf, groSp, mr.scaleGrad(gpe1, -pe1Steps(i-1)), gpe2jr, lblIncPar);
        end
        seq.addBlock(adc, gro1, mr.scaleGrad(gpe1, pe1Steps(i)), gpe2je);
    end
    seq.addBlock(groSp, mr.scaleGrad(gpe1, -pe1Steps(i)), gpe2jr, ...
        mr.makeDelay(TRoutDelay), mr.makeLabel('INC','LIN', PEsamp_INC(count)), lblResetPar);
end
fprintf('Sequence ready (blocks generation took %g seconds)\n', toc);

%% Timing check
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Write
seq.setDefinition('FOV', fov);
seq.setDefinition('Name', 'mp_gt');
seq.setDefinition('ReadoutOversamplingFactor', ro_os);
seq.setDefinition('OrientationMapping', 'SAG');
seq.setDefinition('kSpaceCenterLine', centerLineIdx-1);
seq.setDefinition('PhaseResolution', phaseResoluion);

seq.write('mprage_grappa.seq');
fprintf('Written: mprage_grappa.seq\n');
