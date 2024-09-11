function [NGLDM , Q] = ComputeNGLDM(I, d, a, Nlevel)
% FUNCTION EXTRACTS 14 FEATURES STORED IN NGLDM STRUCT 
% I is a double image, d = radius of neighborhood, a = threshold, usually set to 0
% QUANTIZE I
I = double(I);
%I = I/max(I(:));
quant_levels = linspace(0,1,Nlevel); quant_levels = quant_levels(2:end);
I_discrete = imquantize(I,quant_levels);
clearvars quant_levels
% COMPUTE NGLDM MATRIX
L_neigh = 2*d+1;
N_pix_neigh = L_neigh*L_neigh;
middleIdx = ceil(N_pix_neigh);
B = im2col(I_discrete,[L_neigh,L_neigh]); % convert I into column-wise arrangement, where each col contains all elements in a (L_neigh x L_neigh) bounding box
Cdiff = abs(bsxfun(@minus, B, B(middleIdx,:)));
C = (Cdiff<=a);
D = sum(C,1) - 1;
Ng = B(middleIdx,:).' + 1; % Ng lists all the gray-tone levels      
Nr = D.' + 1; % Nr lists the co-occurence frequency at corresponding gray-tones
Q = accumarray([Ng,Nr], 1, [Nlevel, N_pix_neigh]);
% COMPUTE TEXTURAL FEATURES
[rr,ss] = meshgrid(1:N_pix_neigh,1:Nlevel);
ss_squared = ss.*ss;
rr_squared = rr.*rr;

tmp = Q./ss_squared;
Qsum = sum(Q(:));                                                          % normalization factor
NGLDM.LDE = sum(tmp(:))/Qsum;                                              % (1) low dependence emphasis
tmp = Q.*ss_squared;
NGLDM.HDE = sum(tmp(:))/Qsum;                                              % (2) high dependence emphasis

tmp = Q./rr_squared;
NGLDM.LGLCE = sum(tmp(:))/Qsum;                                            % (3) low gray level count emphasis
tmp = Q.*rr_squared;
NGLDM.HGLCE = sum(tmp(:))/Qsum;                                            % (4) high gray level count emphasis

tmp = Q./rr_squared./ss_squared;
NGLDM.LDLGLE = sum(tmp(:))/Qsum;                                           % (5) low dependence low gray level emphasis
tmp = Q.*rr_squared./ss_squared;
NGLDM.LDHGLE = sum(tmp(:))/Qsum;                                           % (6) low dependence high gray level emphasis
tmp = Q./rr_squared.*ss_squared;
NGLDM.HDLGLE = sum(tmp(:))/Qsum;                                           % (7) high dependence low gray level emphasis
tmp = Q.*rr_squared.*ss_squared;
NGLDM.HDHGLE = sum(tmp(:))/Qsum;                                           % (8) high dependence high gray level emphasis

NGLDM.GLNUN = sum(sum(Q,2).^2)/Qsum/Qsum;                                  % (9) normalized gray level count non-uniformity
NGLDM.DCNUN = sum(sum(Q,1).^2)/Qsum/Qsum;                                  % (10) normalized dependence count non-uniformity

Qnorm = Q/Qsum;
mu_gl = sum(sum(rr.*Qnorm));
mu_ct = sum(sum(ss.*Qnorm));
NGLDM.GLV = sum(sum((rr-mu_gl).^2.*Qnorm));                                 % (11) gray level variance
NGLDM.DCV = sum(sum((ss-mu_ct).^2.*Qnorm));                                 % (12) dependence count variance
NGLDM.DCE = entropy(Qnorm);                                                % (13) dependence count entropy
NGLDM.DCENE = sum(sum(Qnorm.*Qnorm));                                      % (14) dependence count energy
end