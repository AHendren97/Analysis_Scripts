data    = nan(1,13);
data(1) = output.peak;
if (output.peak < 0)
    warning('IMPOSSIBLE NEGATIVE PEAK VOLTAGE')
end
data(2) = output.afterhyp;
if (output.afterhyp > 0)
    warning('IMPOSSIBLE POSITIVE AFTERHYPERPOLARISATION')
end
data(3) = output.amp;
if (output.amp < 0)
    warning('IMPOSSIBLE NEGATIVE AP AMPLITUDE')
end
data(4) = output.thresh;
if (output.thresh > 0)
    warning('IMPOSSIBLE POSITIVE THRESHOLD POTENTIAL')
end
data(5) = output.half;
if (output.half < 0)
    warning('IMPOSSIBLE NEGATIVE HALF VALUE')
end
data(6) = output.rise;
if (output.rise < 0)
    warning('IMPOSSIBLE NEGATIVE RISE')
end
data(7) = output.fall;
if (output.fall > 0)
    warning('IMPOSSIBLE POSITIVE FALL VALUE')
end
data(8) = output.IR;
if (output.IR < 0)
    warning('IMPOSSIBLE NEGATIVE INPUT RESISTANCE')
end
data(9) = output.Rh;
if (output.Rh < 0)
    warning('IMPOSSIBLE NEGATIVE RHEOBASE')
end
data(10) = output.sag_mV(1);
if ~strcmpi(output.Vm_stability,'stable')
    warning ('UNSTABLE RESTING MEMBRANE POTENTIAL');
end
data(11) = output.Vm;
if (output.Vm > 0)
    warning('IMPOSSIBLE POSITIVE RMP')
end
data(12) = max(output.numSpikes);
if (max(output.numSpikes) < 0)
    warning('IMPOSSIBLE NEGATIVE MAX FREQUENCY')
end
data(13) = output.adaptIndex;
if (output.adaptIndex < 0)
    warning('IMPOSSIBLE NEGATIVE ADAPTATION INDEX')
end