% 260321_write_mprage.m
% Patched copy of extern/pulseq/matlab/demoSeq/writeMPRAGE.m
% Removes seq.plot() calls for headless batch execution on HPC.
%
% Generates: mprage.seq (fully sampled 3D MPRAGE, sagittal, 192x240x256)
%
% Usage (after sourcing env.sh):
%   matlab -batch "addpath('extern/pulseq/matlab'); run('260321_write_mprage.m')"

addpath(fullfile(fileparts(mfilename('fullpath')), 'extern', 'pulseq', 'matlab'));

% set system limits (slew rate 130 and max_grad 30 work on Prisma; Trio requires 20us RF-ringdown)
sys = mr.opts('MaxGrad', 24, 'GradUnit', 'mT/m', ...
    'MaxSlew', 100, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 20e-6, ...
    'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

seq=mr.Sequence(sys);           % Create a new sequence object
alpha=7;                        % flip angle
ro_dur=5017.6e-6; % BW=200Hz/pix
ro_os=1;                        % readout oversampling
ro_spoil=3;                     % additional k-max excursion for RO spoiling
TI=1.1;
TRout=2.5;
rfSpoilingInc=84;              % RF spoiling increment
rfLen=100e-6;
ax=struct; % encoding axes
% sagittal
fov=[192 240 256]*1e-3;         % Define FOV and resolution
N = [192 240 256];              % matrix sizes
ax.d1='z'; % the fastest dimension (readout)
ax.d2='x'; % the second-fast dimension (the inner pe loop)

ax.d3=setdiff('xyz',[ax.d1 ax.d2]); % automatically set the slowest dimension
ax.n1=strfind('xyz',ax.d1);
ax.n2=strfind('xyz',ax.d2);
ax.n3=strfind('xyz',ax.d3);

%% Pulses and gradients

rf = mr.makeBlockPulse(alpha*pi/180,sys,'Duration',rfLen, 'use', 'excitation');
rf180 = mr.makeAdiabaticPulse('hypsec',sys,'Duration',10.24e-3,'dwell',1e-5,...
    'use', 'inversion');

deltak=1./fov;
gro = mr.makeTrapezoid(ax.d1,'Amplitude',N(ax.n1)*deltak(ax.n1)/ro_dur,'FlatTime',ceil((ro_dur+sys.adcDeadTime)/sys.gradRasterTime)*sys.gradRasterTime,'system',sys);
adc = mr.makeAdc(N(ax.n1)*ro_os,'Duration',ro_dur,'Delay',gro.riseTime,'system',sys);
groPre = mr.makeTrapezoid(ax.d1,'Area',-gro.amplitude*(adc.dwell*(adc.numSamples/2+0.5)+0.5*gro.riseTime),'system',sys);
gpe1 = mr.makeTrapezoid(ax.d2,'Area',-deltak(ax.n2)*(N(ax.n2)/2),'system',sys);
gpe2 = mr.makeTrapezoid(ax.d3,'Area',-deltak(ax.n3)*(N(ax.n3)/2),'system',sys);
gslSp = mr.makeTrapezoid(ax.d3,'Area',max(deltak.*N)*4,'Duration',10e-3,'system',sys);
[gro1,groSp]=mr.splitGradientAt(gro,gro.riseTime+gro.flatTime);
if ro_spoil>0
    groSp=mr.makeExtendedTrapezoidArea(gro.channel,gro.amplitude,0,deltak(ax.n1)/2*N(ax.n1)*ro_spoil,sys);
end

rf.delay=mr.calcDuration(groSp,gpe1,gpe2);
[groPre,~,~]=mr.align('right',groPre,gpe1,gpe2);
gro1.delay=mr.calcDuration(groPre);
adc.delay=gro1.delay+gro.riseTime;
gro1=mr.addGradients({gro1,groPre},'system',sys);
TRinner=mr.calcDuration(rf)+mr.calcDuration(gro1);

pe1Steps=((0:N(ax.n2)-1)-N(ax.n2)/2)/N(ax.n2)*2;
pe2Steps=((0:N(ax.n3)-1)-N(ax.n3)/2)/N(ax.n3)*2;

TIdelay=round((TI-(find(pe1Steps==0)-1)*TRinner-(mr.calcDuration(rf180)-mr.calcRfCenter(rf180)-rf180.delay)-rf.delay-mr.calcRfCenter(rf))/sys.blockDurationRaster)*sys.blockDurationRaster;
TRoutDelay=TRout-TRinner*N(ax.n2)-TIdelay-mr.calcDuration(rf180);

lblIncLin=mr.makeLabel('INC','LIN', 1);
lblIncPar=mr.makeLabel('INC','PAR', 1);
lblResetPar=mr.makeLabel('SET','PAR', 0);

%% Pre-register and build sequence

gslSp.id=seq.registerGradEvent(gslSp);
groSp.id=seq.registerGradEvent(groSp);
gro1.id=seq.registerGradEvent(gro1);
[~, rf.shapeIDs]=seq.registerRfEvent(rf);
[rf180.id, rf180.shapeIDs]=seq.registerRfEvent(rf180);
lblIncPar.id=seq.registerLabelEvent(lblIncPar);

tic;
for j=1:N(ax.n3)
    seq.addBlock(rf180);
    seq.addBlock(TIdelay,gslSp);
    rf_phase=0;
    rf_inc=0;
    gpe2je=mr.scaleGrad(gpe2,pe2Steps(j));
    gpe2je.id=seq.registerGradEvent(gpe2je);
    gpe2jr=mr.scaleGrad(gpe2,-pe2Steps(j));
    gpe2jr.id=seq.registerGradEvent(gpe2jr);
    for i=1:N(ax.n2)
        rf.phaseOffset=rf_phase/180*pi;
        adc.phaseOffset=rf_phase/180*pi;
        rf_inc=mod(rf_inc+rfSpoilingInc, 360.0);
        rf_phase=mod(rf_phase+rf_inc, 360.0);
        if (i==1)
            seq.addBlock(rf);
        else
            seq.addBlock(rf,groSp,mr.scaleGrad(gpe1,-pe1Steps(i-1)),gpe2jr,lblIncPar);
        end
        seq.addBlock(adc,gro1,mr.scaleGrad(gpe1,pe1Steps(i)),gpe2je);
    end
    seq.addBlock(groSp,mr.makeDelay(TRoutDelay),lblResetPar,lblIncLin);
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
seq.setDefinition('Name', 'mprage');
seq.setDefinition('OrientationMapping', 'SAG');
seq.setDefinition('ReceiverGainHigh',1);

seq.write('mprage.seq');
fprintf('Written: mprage.seq\n');
