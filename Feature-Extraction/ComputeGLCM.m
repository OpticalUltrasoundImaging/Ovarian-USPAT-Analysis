function [GLSDM, Q] = ComputeGLCM(I,Nlevel)
% FUNCTION EXTRACTS 13 FEATURES STORED IN GLSDM STRUCT
% COMPUTE GLCM
I = double(I);
%I = I/max(I(:));
offsets = [0 1; -1 1; -1 0; -1 -1];   % offsets are in [0 45 90 135]
Q = graycomatrix(I,'NumLevels',Nlevel,'Offset',offsets);
Q = mean(Q,3);
Q = Q/sum(Q(:)); % normalize Q
% EXTRACT GLCM FEATURES
sizeQ = size(Q);
[col,row] = meshgrid(1:sizeQ(1), 1:sizeQ(2));
mu_row = sum(row(:).*Q(:));         std_row = sqrt(sum( (row(:)-mu_row).^2.*Q(:) ));
mu_col = sum(col(:).*Q(:));         std_col = sqrt(sum( (col(:)-mu_col).^2.*Q(:) ));
rowQ = sum(Q,2); colQ = sum(Q,1);
% p_x+y
Q90 = rot90(Q);
p_xpy = zeros(2*sizeQ(1)-1 , 1);
k = 1;
for idx = (1-sizeQ(1)):(sizeQ(1)-1)
    p_xpy(k) = sum(diag(Q90,idx));
    k = k+1;
end
% p_x-y
p_xmy = zeros(sizeQ(1),1);    p_xmy(1) = sum(diag(Q,0));
k = 2;
for idx = 1:(sizeQ(1)-1)
    p_xmy(k) = sum([diag(Q,idx);diag(Q,-idx)]);
    k = k+1;
end

GLSDM.ASM = sum(Q(:).^2);                                                   % (1) angular second moment or energy
tmp = (abs(row-col).^2).*Q;
GLSDM.CONTRAST = sum(tmp(:));                                               % (2) contrast
tmp = sum( (row(:)-mu_row).*(col(:)-mu_col).*Q(:) );
GLSDM.CORR = tmp/std_row/std_col;                                           % (3) correlation
GLSDM.VAR = sum( (row(:)-mean(Q(:))).^2 .*Q(:) );                           % (4) variance
GLSDM.IDM = sum( Q(:)./(1+(row(:) - col(:)).^2) );                          % (5) inverse difference moment
GLSDM.SAV = sum( (2:(2*sizeQ(1)))'.*p_xpy );                                % (6) sum average
GLSDM.SENT = - sum( p_xpy(p_xpy~=0).*log(p_xpy(p_xpy~=0)));                 % (7) sum entropy
GLSDM.SVAR = sum(((2:(2*sizeQ(1)))' - GLSDM.SENT).^2 .*p_xpy);              % (8) sum variance
GLSDM.ENT = - sum( Q(Q~=0).*log2(Q(Q~=0)));                                 % (9) entropy
GLSDM.DVAR = sum(((0:sizeQ(1)-1)' - mean(p_xmy)).^2.* p_xmy);               % (10) difference variance
GLSDM.DENT = - sum( p_xmy(p_xmy~=0).*log(p_xmy(p_xmy~=0)));                 % (11) difference entropy

logrc = log2(rowQ*colQ);
logrc(isinf(logrc)) = 0;
HXY1 = -sum(Q(:).*logrc(:));
numerator = GLSDM.ENT - HXY1;
logc = log2(colQ); logc(isinf(logc)) = 0;    HX = -sum(colQ.*logc);
logr = log2(rowQ); logr(isinf(logr)) = 0;    HY = -sum(rowQ.*logr);
GLSDM.IMC1 = numerator/max([HX,HY]);                                        % (11) information measure of correlation 1
HXY2 = -sum(sum((rowQ*colQ).*logrc));
GLSDM.IMC2 = ((1-exp(-2*(HXY2-GLSDM.ENT)))).^(1/2);                         % (12) information measure of correlation 2
qtmp = zeros(sizeQ(1),sizeQ(2));
for idx = 1:sizeQ(2)
    qtmp(idx,:) = sum((repmat(Q(idx,:),sizeQ(1),1).*Q)./...
        repmat(rowQ(idx).*colQ, sizeQ(1),1),2,"omitnan");
end
qtmp(isnan(qtmp))=0;
eigenvec = eig(qtmp);
eigenvec(eigenvec==max(eigenvec))=[];
if ~any(eigenvec) || isempty(eigenvec)
    GLSDM.MAXCORR = 0;
else
    GLSDM.MAXCORR = abs(sqrt(max(eigenvec)));                               % (13) maximum correlation coefficient
end
end



