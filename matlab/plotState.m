function plotState(os, imax, jmax, kmax, depth, layerDepths);

global h1 h2 h3;

if isempty(h1)
    h1 = figure;
    h2 = figure;
end
try
    figure(h1);
catch
    h1 = figure;
    h2 = figure;
    h3 = figure;
    figure(h1);
end
    
inds = depth==0;
imid = floor(imax/2);
jmid = floor(jmax/2);
xval = 1:imax;
yval = 1:jmax;
%snitt = permute(os.T(imid,:,:),[3 2 1]);
snitt = permute(os.T(:,jmid,:),[3 1 2]);
subplot(3,3,1), pcolor(xval, layerDepths, snitt), colorbar, shading flat
set(gca,'YDir','reverse');
title('Temperature');
%snitt = permute(os.S(imid,:,:),[3 2 1]);
snitt = permute(os.S(:,jmid,:),[3 1 2]);
subplot(3,3,2), pcolor(xval, layerDepths, snitt), colorbar, shading flat
set(gca,'YDir','reverse');
title('Salinity');
pad = 1;
inds2 = depth(1+pad:end-pad,1+pad:end-pad)==0;
snitt = permute(os.T(1+pad:end-pad,1+pad:end-pad,1),[2 1 3]);
snitt(inds2') = NaN;
subplot(3,3,4), pcolor(snitt), colorbar, shading flat
title('Temperature k=1');
snitt = permute(os.S(1+pad:end-pad,1+pad:end-pad,1),[2 1 3]);
snitt(inds2') = NaN;
subplot(3,3,5), pcolor(snitt), colorbar, shading flat
title('Salinity k=1');

inds2 = depth(1+pad:end-pad,1+pad:end-pad)==0;
snitt = permute(os.T(1+pad:end-pad,1+pad:end-pad,4),[2 1 3]);
snitt(inds2') = NaN;
subplot(3,3,7), pcolor(snitt), colorbar, shading flat
title('Temperature k=4');
snitt = permute(os.S(1+pad:end-pad,1+pad:end-pad,4),[2 1 3]);
snitt(inds2') = NaN;
subplot(3,3,8), pcolor(snitt), colorbar, shading flat
title('Salinity k=4');

if ~isempty(os.X)
    snitt = permute(os.X(1+pad:end-pad,1+pad:end-pad,1),[2 1 3]);
    snitt(inds2') = NaN;
    subplot(3,3,6), pcolor(snitt), colorbar, shading flat
    %set(gca,'YDir','reverse');
    title('Tracer k=1');
    %caxis([0 0.01]);
    
    snitt = permute(os.X(1+pad:end-pad,1+pad:end-pad,4),[2 1 3]);
    snitt(inds2') = NaN;
    subplot(3,3,9), pcolor(snitt), colorbar, shading flat
    %set(gca,'YDir','reverse');
    title('Tracer k=4');
    %caxis([0 0.01]);    
end
subplot(3,3,3),
% snitt = permute(os.U(:,jmid,:),[3 1 2]);
% pcolor(xval(1:end-1), layerDepths, snitt), colorbar, shading flat
% set(gca,'YDir','reverse');
% title('U current');
% ca = caxis;
% caxis([-max(abs(ca)) max(abs(ca))]);
snitt = os.E';
snitt(inds') = NaN;
pcolor(snitt), colorbar, shading flat;
title('Elevation');
subs = 1;
border=1;
quivMult=2;
%kval = [1 5 8 10];
kval = [1 2 3 5];
contourLevs = [50 75];%[200 500 1000 1500];
for i=1:length(kval)
    k = kval(i);
    if k > kmax
        continue;
    end
    figure(h2);
    subplot(2,2,i);
    [uu,vv, t1, t2] = interpolateUV(os.U(:,:,k), os.V(:,:,k), 1,border);
    [uus,vvs, t1s, t2s] = interpolateUV(os.U(:,:,k), os.V(:,:,k), subs,border);
    dpt = depth(border+1:subs:end-border,border+1:subs:end-border);
   
    %inds = uu==0 & vv==0;
    if k>1
        inds = dpt<layerDepths(k-1);
    else
        inds = dpt<=0;
    end
    %uu(inds) = NaN;
    %vv(inds) = NaN;
    
    hold off;
    pcolor(t1', t2', sqrt(uu.^2 + vv.^2)'), shading flat, colorbar
    cax = caxis;
    %hold on, quiver(t1'+subs*0.5,t2'+subs*0.5, uu',vv',4,'k'); %mesh(U(:,:,1)');
    hold on, quiver(t1s'+subs*0.5,t2s'+subs*0.5, uus',vvs',quivMult,'k'); %mesh(U(:,:,1)');
    %[contr, hh] = contour(dpt',contourLevs,'k');
    caxis(cax);
    title(['Current k=' num2str(k)]);
    axis image
%     figure(h3);
%     subplot(2,2,i);
%     ww = os.W(border+1:subs:end-border,border+1:subs:end-border, k);
%     pcolor(t1',t2',ww'), shading flat, colorbar
%     hold on
%     %[contr, hh] = contour(dpt',contourLevs,'k');
%     %caxis([-0.05 0.05]);
%     title(['Vertical speed k=' num2str(k)]);
    
end

% figure(h2);
% subplot(2,2,4);
% hold off;
% pcolor(t1', t2', sqrt(os.windU.^2 + os.windV.^2)'), shading flat, colorbar
% hold on, quiver(t1'+subs*0.5,t2'+subs*0.5, os.windU',os.windV',2,'k');
% title('Wind');
% border = 1;
% [uu,vv, t1, t2] = interpolateUV(os.U(:,:,1), os.V(:,:,1), subs,border);
% pcolor(t1', t2', atan2(uu', vv')), shading flat, colorbar


%figure,quiver(os.windU(:,1:end-1)', os.windV(1:end-1,:)')




