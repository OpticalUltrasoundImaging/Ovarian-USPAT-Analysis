function [weighted_mean, fwhm] = FindFWHM(x)

x(isinf(x)) = 0;
xmax = max(x(:));  xmin = min(x(:));
xnorm = (x-xmin)/(xmax-xmin);   % NORMALIZE TO (0,1)
N_bin = 50;
[binidx,edges] = discretize(xnorm,N_bin);
dedges = edges(2) - edges(1);
hist = zeros(N_bin,1);
for i = 1:N_bin
hist(i) = sum(sum(binidx==i));
end
hist = hist/max(hist(:));

[~,histmaxidx] = max(hist(:));
fwhmidx1 = find(flipud(hist(1:histmaxidx))>=0.5,1,'first');
if isempty(fwhmidx1)
    fwhmidx1 = 1;
end
fwhmidx2 = find((hist(histmaxidx:end))<=0.5,1,'first');
if isempty(fwhmidx2)
    fwhmidx2 = N_bin;
end
fwhm = (fwhmidx2 - fwhmidx1 + 1)*dedges*(xmax-xmin) + xmin;
xseg = (fwhmidx1:fwhmidx2)*dedges*(xmax-xmin) + xmin;
weighted_mean = xseg*hist(fwhmidx1:fwhmidx2)/sum(hist(fwhmidx1:fwhmidx2));

end
