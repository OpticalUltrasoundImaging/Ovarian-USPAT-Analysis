function Zmoment = computeShapeMoments(d)
bw_mask = poly2mask(d(:,1),d(:,2),1024,2048);
x_range = max(d(:,1)) - min(d(:,1));
z_range = max(d(:,2)) - min(d(:,2));
L_roi = max(x_range,z_range) + 50;

if x_range >= z_range
    x_roi = [max(1,min(d(:,1))-25) , min(2048,max(d(:,1))+25)] ;
    z_os = round((L_roi - z_range)/2);
    z_roi = [max(1,min(d(:,2))-z_os) , min(1024,max(d(:,2))+z_os)] ;
else
    z_roi = [max(1,min(d(:,2))-25) , min(1024,max(d(:,2))+25)] ;
    x_os = round((L_roi - x_range)/2);
    x_roi = [max(1,min(d(:,1))-x_os) , min(2048,max(d(:,1))+x_os)] ; 
end
I = bw_mask(z_roi(1):z_roi(2) , x_roi(1):x_roi(2));
%% COMPUTE ZERNIKE MOMENTS (3)
N = 128;
I_zernike = imresize(I,[N,N]);
I_zernike = I_zernike/max(I_zernike(:));
x = 1:N; y = x;
[X,Y] = meshgrid(x,y);
R = sqrt((2.*X-N-1).^2+(2.*Y-N-1).^2)/N;
Theta = atan2((N-1-2.*Y+2),(2.*X-N+1-2));
R = (R<=1).*R;
cnt = nnz(R)+1;
ZfeatureList = [4,2; 6,0; 9,7]; % (4,2) OVAL, (6,0) IRREGULARITY, (9,7) BOUNDARY DISTINCTNESS
Zmoment = zeros(1,size(ZfeatureList,1));
for i_feature = 1:size(ZfeatureList,1)
    n = ZfeatureList(i_feature,1); m = ZfeatureList(i_feature,2);
    rad = zeros(size(R));
    for s = 0:(n-abs(m))/2
        c = (-1)^s*factorial(n-s)/(factorial(s)*factorial((n+abs(m))/2-s)*factorial((n-abs(m))/2-s));
        rad = rad + c*R.^(n-2*s);
    end
    Product = I_zernike(x,y).*rad.*exp(-1i*m*Theta);
    Z = sum(Product(:));
    Z = (n+1)*Z/cnt;
    Zmoment(i_feature) = abs(Z);
end
end