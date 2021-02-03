function plotVerticalMixing(os, sp, depth, layerDepths);

xval = 1:sp.imax;
yval = 1:sp.jmax;
jmid = floor(sp.jmax/2);


density = permute(os.rho(:,jmid,:),[3 1 2]);
dryCells = density<100;
density(dryCells) = NaN;
figure,
subplot(2,2,1); 
mix = permute(os.K_v(:,jmid,:),[3 1 2]);
mix(dryCells(2:end,:)) = NaN;
pcolor(xval, layerDepths(2:end), mix), shading flat, colorbar;
set(gca,'YDir','reverse'), title('Vertical mixing coeff');
subplot(2,2,2); 

pcolor(xval, layerDepths, density), shading flat, colorbar;
set(gca,'YDir','reverse'), title('Density');
subplot(2,2,3); 
value = permute(os.T(:,jmid,:),[3 1 2]); value(dryCells) = NaN;
pcolor(xval, layerDepths, value), shading flat, colorbar;
set(gca,'YDir','reverse'), title('Temperature');
subplot(2,2,4); 
value = permute(os.S(:,jmid,:),[3 1 2]); value(dryCells) = NaN;
pcolor(xval, layerDepths, value), shading flat, colorbar;
set(gca,'YDir','reverse'), title('Salinity');

