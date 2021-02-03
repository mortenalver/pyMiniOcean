function plotVerticalMixingFromFile(filename, sample);


[os, nSamplesInFile, depth, layerDepths] = loadState(filename, sample);
nSamplesInFile
sz = size(os.T);
imax = sz(1);
jmax = sz(2);

kmax = sz(3);
% Prepare OceanState and setup variables to complete the state:
dz = diff([0;layerDepths]);
init_simparam;
adaptDepthField;
sp.imax = imax;
sp.jmax = jmax;
sp.kmax = kmax;
os = os.updateCellHeights(kmm, dzz, sp);
os = os.updateRho(kmm, sp);
os = os.calcVerticalMixingCoeffRichardson(kmm, sp);
% done

plotVerticalMixing(os, sp, depth, layerDepths);
