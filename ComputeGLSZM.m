function [textures, GLSZM] = ComputeGLSZM(I, Nlevel)
% FUNCTION EXTRACTS 13 FEATURES STORED IN textures STRUCT
I = double(I);
%I = I/max(I(:));
quant_levels = linspace(0,1,Nlevel); quant_levels = quant_levels(2:end);
I_discrete = imquantize(I,quant_levels); % I_discrete has integer values from 1 to Nlevel
clearvars quant_levels
% INITIALIZATION
nInit = numel(I_discrete);
GLSZM = zeros(Nlevel,nInit);
[H,W] = size(I_discrete);
%% COMPUTE GLSZM
%fwb = waitbar(0,'COMPUTING GLSZM ...');
for i = 1:Nlevel
    temp = zeros(H,W);
    temp(I_discrete == i) = 1;
    connObjects = bwconncomp(temp,8);
    nZone = length(connObjects.PixelIdxList);
    for j = 1:nZone
        col = length(connObjects.PixelIdxList{j});
        GLSZM(i,col) = GLSZM(i,col) + 1;
    end
    %perc = i/Nlevel;
    %msg = sprintf('COMPUTING GLSZM LEVEL %d / %d.',i,Nlevel);
    %waitbar(perc,fwb,msg);
end
%close(fwb)
stop = find(sum(GLSZM),1,'last');
GLSZM(:,(stop+1):end) = []; % REMOVE EMPTY COLUMNS
GLSZM = GLSZM./sum(GLSZM(:)); % NORMALIZE GLSZM
%% COMPUTE GLSZM FEATURES
% USEFUL MATRICES, VECTORS AND QUANTITIES
sz = size(GLSZM); % Size of GLSZM
cVect = 1:sz(2); rVect = 1:sz(1);% Row and column vectors
[cMat,rMat] = meshgrid(cVect,rVect); % Column and row indicators for each entry of the GLSZM
pg = sum(GLSZM,2)'; % Gray-Level Run-Number Vector
pr = sum(GLSZM); % Run-Length Run-Number Vector
ug = (pg*rVect')/(sz(1)*sz(2));
ur = (pr*cVect')/(sz(1)*sz(2));
textures.SZE = pr*(cVect.^(-2))';                                           % (1) small zone emphasis
textures.LZE = pr*(cVect.^2)';                                              % (2) large zone emphasis
textures.GLN = sum(pg.^2);                                                  % (3) gray-level non-uniformity
textures.ZSN = sum(pr.^2);                                                  % (4) zone-size non-uniformity
textures.ZP = sum(pg)/(pr*cVect');                                          % (5) zone percentage
textures.LGZE = pg*(rVect.^(-2))';                                          % (6) low gray-level zone emphasis
textures.HGZE = pg*(rVect.^2)';                                             % (7) high gray-level zone emphasis
textures.SZLGE = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^(-2))));               % (8) SMALL-zone LOW-gray-level emphasis
textures.SZHGE = sum(sum(GLSZM.*(rMat.^2).*(cMat.^(-2))));                  % (9) SMALL-zone HIGH-gray-level emphasis
textures.LZLGE = sum(sum(GLSZM.*(rMat.^(-2)).*(cMat.^2)));                  % (10) LARGE-zone LOW-gray-level emphasis
textures.LZHGE = sum(sum(GLSZM.*(rMat.^2).*(cMat.^2)));                     % (11) LARGE-zone HIGH-gray-level emphasis
GLV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        GLV = GLV + (GLSZM(g,r)*g-ug)^2;
    end
end
textures.GLV = GLV/(sz(1)*sz(2));                                           % (12) gray level variance
ZSV = 0;
for g = 1:sz(1)
    for r = 1:sz(2)
        ZSV = ZSV + (GLSZM(g,r)*r-ur)^2;
    end
end
textures.ZSV = ZSV/(sz(1)*sz(2));                                           % (13) zone size variance

end