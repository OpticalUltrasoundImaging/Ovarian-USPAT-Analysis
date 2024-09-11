function [Psi , G] = computeShapeIndices(d,L_US)
% FUNCTION FIRST COMPUTES GEOMETRICAL PROPERTIES SAVED IN STRUCT G
% THEN EXTRACTS 10 GEOMETRICAL FEATURES SAVED IN STRUCT PSI
%% EXTRACT GEOMETRICAL PARAMETERS
N_pts = size(d,1);
z_coord = d(:,2)/1024*(1+L_US);      % [cm]
x_coord = (d(:,1)-1024)/1024*(1+L_US); % [cm]
z_coord_tmp = [z_coord(2:end);z_coord(1)];
x_coord_tmp = [x_coord(2:end);x_coord(1)];
G.P = sum(sqrt((z_coord-z_coord_tmp).^2+(x_coord-x_coord_tmp).^2)); % P = perimeter
clearvars z_coord_tmp x_coord_tmp
G.A = polyarea(x_coord , z_coord); % A = area

CH = convhull(x_coord , z_coord);
x_coord_convhull = x_coord(CH);
z_coord_convhull = z_coord(CH);
z_coord_tmp = z_coord_convhull(2:end);
x_coord_tmp = x_coord_convhull(2:end);
G.P_convhull = sum(sqrt((z_coord_convhull(1:end-1)-z_coord_tmp).^2+(x_coord_convhull(1:end-1)-x_coord_tmp).^2)); % P_convhull = perimeter of convex hull
G.A_convhull = polyarea(x_coord_convhull , z_coord_convhull); % A = area of the convex hull
clearvars z_coord_tmp x_coord_tmp

figure
% scatter(x_coord,z_coord,'MarkerFaceColor','r','MarkerEdgeColor','r')
plot([x_coord; x_coord(1)],[z_coord; z_coord(1)],'MarkerFaceColor','r','MarkerEdgeColor','r')
set(gca,'xtick','','ytick','')
hold on
%scatter(x_coord_convhull,z_coord_convhull,'MarkerFaceColor','b','MarkerEdgeColor','b')
plot(x_coord_convhull,z_coord_convhull,'color','b','LineWidth',2)
hold off
set(gca,'xtick','','ytick','')

[xc,zc] = centroid(polyshape(x_coord, z_coord)); % [xc,zc] is the centroid
R_hist = sqrt((xc-x_coord).^2 + (zc-z_coord).^2);
G.Rmax = max(R_hist); G.Rmin = min(R_hist);

bw_mask = poly2mask(d(:,1),d(:,2),1024,2048);
Dferet = bwferet(bw_mask,"MaxFeretProperties");
EDferet = bwferet(bw_mask,"MinFeretProperties");
%bw_mask = imresize(bw_mask,[1024,round(2048/2/L_US*(1+L_US))]);
s = regionprops(bw_mask,{'MajorAxisLength','MinorAxisLength','Orientation'});
d_bw = (1+L_US)/1024; % [cm]
G.Lma = s.MajorAxisLength*d_bw; % Lma = major axis length
G.Lsa = s.MinorAxisLength*d_bw; % Lsa = minor axis length
G.D = Dferet.MaxDiameter * d_bw;
G.ED = EDferet.MinDiameter * d_bw;

tf_dist = bwdist(~bw_mask);
rho_inscribe = max(tf_dist(:));
[zcc , xcc] = find(tf_dist == rho_inscribe);
G.rho_i = rho_inscribe*d_bw*1.0; % rho_inscribe
zcc = zcc*d_bw;
xcc = (xcc - ceil(size(bw_mask,2)/2))*d_bw;
[rho_circumscribe,CC_circumscribe,~] = FindMinEncloseCircle([x_coord_convhull , z_coord_convhull]); % rho_circumscribe
G.rho_e = rho_circumscribe;
% figure
% scatter(x_coord,z_coord,'MarkerFaceColor','r','MarkerEdgeColor','r')
% hold on
% scatter(xc,zc,'MarkerEdgeColor','b')
% viscircles(CC_circumscribe, G.rho_e,'Color','b')
% viscircles([xcc, zcc], G.rho_i,'Color','b')
% hold off
% set(gca,'xtick','','ytick','')
%% EXTRACT SHAPE INDICES
Psi.R = pi * G.Rmax * G.Rmin / G.A;
Psi.MA = pi/4 * G.Lma * G.Lsa / G.A;
Psi.D = pi/2 * G.ED * G.D / G.A;
z_pix_convhull = round(z_coord_convhull*1024/(1+L_US));
x_pix_convhull = round(x_coord_convhull*1024/(1+L_US)+1024);
bw_mask_ch = poly2mask(x_pix_convhull, z_pix_convhull,1024,2048);
bw_mask_not = (bw_mask_ch - bw_mask)>0;
bw_mask_not = imopen(bw_mask_not,strel('disk',1));
%figure; imagesc(bw_mask_not)
CC = bwconncomp(bw_mask_not);
numPixels = cellfun(@numel,CC.PixelIdxList);
Ncc = sum(numPixels>100);
Psi.Ncc = 1/(1+Ncc);

Psi.Extension1 = G.ED / G.D;
Psi.Extension2 = G.rho_i/G.rho_e;
Psi.Deficit = 1 - pi*(G.rho_e - G.rho_i)^2/G.P/G.P;
Psi.Circularity = G.Rmax / G.Rmin;
Psi.Convexity1 = G.A / G.A_convhull;
Psi.Convexity2 = G.P_convhull / G.P;

end
