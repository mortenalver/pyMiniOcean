
function [os, nSamples, time, depth, layerDepths] = loadState(filename, sample);

os = OceanState();

ncid = netcdf.open(filename,'NOWRITE');
timeD = netcdf.inqDimID(ncid,'time');
[~, nSamples] = netcdf.inqDim(ncid, timeD);
if sample < 0
   sample = nSamples-1;
end
xcD = netcdf.inqDimID(ncid,'xc');
[~, imax] = netcdf.inqDim(ncid,xcD);
ycD = netcdf.inqDimID(ncid,'yc');
[~, jmax] = netcdf.inqDim(ncid,ycD);
zcD = netcdf.inqDimID(ncid,'zc');
[~, kmax] = netcdf.inqDim(ncid,zcD);

time = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'time'), [sample], [1]);
os.U = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'U'), [0 0 0 sample], [imax jmax kmax 1]);
os.U = os.U(1:end-1,:,:);
os.V = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'V'), [0 0 0 sample], [imax jmax kmax 1]);
os.V = os.V(:,1:end-1,:);
os.W = zeros(imax,jmax,kmax+1);
try
    os.W(:,:,1:kmax) = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'W'), [0 0 0 sample], [imax jmax kmax 1]);
catch
    warning('No vertical speeds found in file.')
end
os.E = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'E'), [0 0 sample], [imax jmax 1]);
os.T = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'T'), [0 0 0 sample], [imax jmax kmax 1]);
os.S = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'S'), [0 0 0 sample], [imax jmax kmax 1]);
try
    os.X = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'X'), [0 0 0 sample], [imax jmax kmax 1]);
catch
    warning('No passive tracer values found in file.')
end
try
    os.K_v = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'K_v'), [0 0 0 sample], [imax jmax kmax-1 1]);
catch
try
    os.windU = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'windU'), [0 0 sample], [imax jmax 1]);
    os.windU = os.windU(1:end-1,:);
    os.windV = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'windV'), [0 0 sample], [imax jmax 1]);
    os.windV = os.windV(:,1:end-1);
catch
    warning('No wind values found in file.')
    os.windU = zeros(imax,jmax);
    os.windV = zeros(imax,jmax);

end
    
end
layerDepths = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'zc'));
depth = netcdf.getVar(ncid, netcdf.inqVarID(ncid, 'depth'));