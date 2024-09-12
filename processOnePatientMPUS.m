function [] = processOnePatientMPUS(folderName, patientID, all_area_dirs, idx_area, info, parentFolder)
patientNum = sprintf('%03d',str2double(patientID(2:end)));
areaNum = all_area_dirs(idx_area);
saveFolder = strcat(patientNum,'_',areaNum{1,1});
savePath = strcat(parentFolder,'\',saveFolder);
if ~exist(savePath)
    mkdir(parentFolder,saveFolder);
end

%% LOAD RF DATA AND COMPUTE BEAMFORMING
idx_list = [1];
N_repeat = length(idx_list)*4;
R_raw_hist = zeros(info.Nfocus,info.Nsc,N_repeat);
R_env_raw_hist = zeros(info.Nfocus,info.Nsc,N_repeat);
f = [0 0.1 0.1 1]; m = [ 0 0 1 1];
DC_cancle = fir2(64,f,m);
for i = 1:N_repeat
    id_add_tmp = floor(i/length(idx_list)); if mod(i,length(idx_list))==0; id_add_tmp = id_add_tmp - 1; end
    idx_area_tmp = idx_area + id_add_tmp;
    area_id = all_area_dirs(idx_area_tmp);
    pname =string(fullfile(folderName,'\',area_id,'\'));
    idx_tmp = mod(i,length(idx_list)); if idx_tmp == 0; idx_tmp = length(idx_list); end
    ChannelFile = [num2str(idx_list(idx_tmp)),'_layer0_idx',num2str(idx_list(idx_tmp)),'_BDATA_RF.mat'];
    ChannelFile =  fullfile(pname,ChannelFile);
    ChannelData = load (ChannelFile);
    clearvars pname ChannelFile

    info.Nfocus = double(ChannelData.AlignedSampleNum);
    info.d_sample = 0:info.pixel_d:(info.Nfocus-1)*info.pixel_d;

    % REORDERING DATA
    Reconstruction = zeros(ChannelData.AlignedSampleNum, info.N_ch, info.Nsc);
    for sc = 1:info.Nsc
        ti = sprintf('AdcData_scanline%03d_thi0_sa0',sc-1);
        imgtmp = double(ChannelData.(ti));
        %img = [imgtmp(:,:,1);imgtmp(:,:,2)];
        img = reorder_chdata_ch2dly(sc-1, info.N_ele, info.Nsc, info.N_ch, imgtmp);
        Reconstruction(:,:,sc) = img';
    end
    clearvars imgtmp img fwd

    data_total = floor(ChannelData.AlignedSampleNum);
    data_total1 = 1500;
    rx_HalfCh = info.N_ch*0.5;
    rx_ch_mtx = [-rx_HalfCh:rx_HalfCh-1];

    RxMux = zeros(info.Nsc, info.N_ch);
    SCvsEle = info.d_theta/info.one_angle; % ratio of incremental scan angle / tx pitch

    for sc = 1:info.Nsc
        idx = floor((sc-1)*SCvsEle) + 1;
        rx_idx_tmp = idx + rx_ch_mtx;
        rx_idx = rx_idx_tmp((rx_idx_tmp > 0) & (rx_idx_tmp <= info.N_ele));
        RxMux(sc,:) = rx_idx_tmp;
    end
    % DAS BEAMFORMING (RAW DATA -> BEAM DATA
    half_rx_ch = info.N_ch*info.pitch*0.5;
    n_sample = double([0:data_total-1]');
    d_sample = n_sample * info.pixel_d;
    RF_Sum = zeros(data_total,info.Nsc);
    RF_beam = zeros(data_total,info.Nsc);
    fwd = waitbar(0,'Beamforming scanline ... ');
    for sc = 1:info.Nsc
        perc = sc/info.Nsc;
        msg = sprintf('Beamforming scanline # %0.2f',perc);
        waitbar(perc,fwd,msg);
        ChData = Reconstruction(:,:,sc);
        MuxTbl = RxMux(sc,:);
        sc_ang = info.ScanAngle(sc)*pi/180;
        sc_pos = sc_ang*info.ROC;
        RF_scanline_tofadjusted = zeros(data_total,info.N_ch);
        active_ch = [];
        for ch = 1:info.N_ch
            ElePos = MuxTbl(ch);
            if(ElePos > 0 && ElePos <= info.N_ele)
                active_ch = [active_ch,ch];
                tmp = ChData(:,ch);
                fil_tmp = conv(tmp, DC_cancle, 'same');
                % CALCULATE DELAY
                ele_ang = info.EleAngle(ElePos)*pi/180;     ch_ang = abs(sc_ang - ele_ang);
                x = info.ROC*sin(ch_ang);                   z = info.ROC*(1 - cos(ch_ang));
                tx_t= d_sample/info.c;                      rx_t = sqrt(x.^2 + (d_sample + z).^2)/info.c;
                TOF = tx_t + rx_t;
                RxReadPointer = round(TOF*info.fs);
                idx1 = find(RxReadPointer > data_total);    RxReadPointer(idx1) = data_total;
                idx2 = find(RxReadPointer < 1);             RxReadPointer(idx2) = 1;
                rf_tmp = fil_tmp(RxReadPointer);
                rf_tmp(idx1) = 0;                           rf_tmp(idx2) = 0;
                % APODIZATION
                Rx_apod_idx = zeros(data_total, 1);         Rx_apo_r = zeros(data_total, 1);
                Xr = abs(sc_pos - ele_ang*info.ROC);
                h_aper_size =  d_sample/info.RxFnum*0.5;
                idx1 = h_aper_size >= half_rx_ch;           idx2 = find(h_aper_size < half_rx_ch);
                Rx_apod_idx(idx1) = Xr./half_rx_ch;         Rx_apod_idx(idx2) = Xr./h_aper_size(idx2);
                idx3 = find(Rx_apod_idx >= 1);              idx4 = find(Rx_apod_idx < 1);
                Rx_apo_r(idx3) = 0;                         Rx_apo_r(idx4) = 1;
                RF_scanline_tofadjusted(:,ch) = (rf_tmp.*Rx_apo_r)';
            end
        end

        r0 = 2;
        start_left = round(info.N_ch/2-r0);                 start_right = round(info.N_ch/2+r0+1);
        idx_left = round(max(active_ch(1),start_left - (1:data_total)*64/data_total));
        idx_right = round(min(active_ch(end),start_right + (1:data_total)*64/data_total));
        row_sum_abs = zeros(data_total,1);                  row_abs_sum = zeros(data_total,1);
        for t = 1:data_total
            rf_temp = RF_scanline_tofadjusted(t,idx_left(t):idx_right(t));
            row_sum_abs(t) = sum(rf_temp)^2;
            abs_rf = abs(rf_temp);
            row_abs_sum(t) = sum(abs_rf.^2)*length(rf_temp);
        end
        CF = row_sum_abs./row_abs_sum;                      CF(isnan(CF)) = 1;

        sum_tmp = sum(RF_scanline_tofadjusted,2);
        RF_beam(:,sc) = sum_tmp;                            RF_Sum(:,sc) = hilbert(sum_tmp.*CF);
    end
    close(fwd)
    clearvars fwd msg perc
    clearvars idx1 idx2 idx3 idx4 idx_left idx_right h_aper_size fil_tmp ChData CF abs_rf active_ch
    clearvars rf_tmp row_abs_sum row_sum_abs Rx_apo_r Rx_apod_idx rx_ch_mtx rx_idx rx_idx_tmp rx_t RxReadPointer
    clearvars x z Xr t ti start_right start_left
    clearvars tmp TOF tx_t sum_tmp
    RF_env = abs(RF_Sum);
    R_env_raw_hist(:,:,i) = RF_env(1:info.Nfocus,:);
    R_raw_hist(:,:,i) = RF_beam(1:info.Nfocus,:);
    Reconstruction_hist(:,:,:,i) = Reconstruction;
    fprintf('%d / %d.\n',i,N_repeat)
end
clearvars Reconstruction rf_temp idx_list all_area_dirs

%% DYNAMIC RANGING
dB_US = 67;                                     min_dB = 10^(-dB_US/20);
RF_env = abs(RF_Sum);                           RF_env_raw = RF_env;
RF_env_norm = RF_env./max(max(RF_env));
idx = find(RF_env_norm < min_dB);
RF_log = (20/dB_US)*log10(RF_env_norm)+1;
RF_log(idx) = 0;
clearvars idx
MPUS.RF_LOG = RF_log;
% CONVERT FROM POLAR TO CARTESIAN COORDINATES
RF_cartesian_B = polar2cartImgPAT(RF_log,[info.ROC*100,info.ROC*100 + max(d_sample)*100],[-1*abs(sc_ang),abs(sc_ang)])';
f = figure;
imshow(RF_cartesian_B,[],'border','tight'); colormap('gray')
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_B_FULL.tiff');
saveas(gcf,imname)
close(f)

%% COMPUTE FEATURE IN A RECTANGULAR BOUNDING BOX
% THIS BOUNDING BOX IS NOT THE SAME AS THE FINAL ROIS. ONLY USED TO SAVE
% COMPUTATION TIME. THIS BOUNDING BOX IS ALWAYS STRICTLY BIGGER THAN THE
% FINAL ROI
roi_x = [-4.25,4.25];
roi_z = [0.25,4];
if info.Nfocus == 3840;  roi_z(2) = 5.75;
elseif info.Nfocus == 3328; roi_z(2) = 5.25;
elseif info.Nfocus == 5376; roi_z(2) = 8.75;
elseif info.Nfocus == 4352; roi_z(2) = 7.25;
elseif info.Nfocus == 4864; roi_z(2) = 7.75;
elseif info.Nfocus < 2305; roi_z(2) = 3.25; end
RF_cartesian_B2 = polar2cartImgPAT_crop(RF_log,[info.ROC*100,(info.ROC*100 + max(d_sample)*100)],[-1*abs(sc_ang),abs(sc_ang)],...
    roi_x,roi_z+1)';
RF_cartesian_B2 = imresize(RF_cartesian_B2,[512,768]);
f = figure;
imshow(RF_cartesian_B2,[],'border','tight'); colormap('gray')
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_B_ROI.tiff');
saveas(gcf,imname)
close(f)
disp('>>>>>>>> B MODE BEAMFORMING COMPLETE.')
MPUS.roi_x = roi_x; MPUS.roi_z = roi_z;
clearvars RF_cartesian_B2 RF_cartesian_B

%% TAI (TISSUE ATTENUATION)
[fc_map, TAI, fc_map_var] = TAI_CURVE(R_raw_hist,info);
mask = imgaussfilt(GENERATE_US_MASK(RF_log,info),2);
MPUS.mask = mask;
RF_cartesian_tai = polar2cartImgPAT_crop(TAI,[info.ROC*100,info.ROC*100 + max(d_sample)*100],[-1*abs(sc_ang),abs(sc_ang)],...
    roi_x,roi_z+1)';
RF_cartesian_tai = imbilatfilt(imresize(RF_cartesian_tai,[512,768]),0.5,3);
MPUS.FC = fc_map;
MPUS.TAI = TAI;
MPUS.FC_VAR = fc_map_var;
TAIMap = turbo(128);
TAIMap = [[0,0,0];TAIMap];
f = figure;
imshow(RF_cartesian_tai,[0.0,1],'border','tight'); colormap(TAIMap)
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_TAI.tiff');
saveas(gcf,imname)
close(f)
disp('>>>>>>>> TAI CALCULATION COMPLETE.')
clearvars fc_map TAI

%% TSI (TISSUE SCATTERING, NAKAGAMI IMAGING)
naka_map = TSI_CURVE(R_env_raw_hist,info);
MPUS.TSI = naka_map;
RF_cartesian_tsi = polar2cartImgPAT_crop(naka_map,[info.ROC*100,info.ROC*100 + max(d_sample)*100],[-1*abs(sc_ang),abs(sc_ang)],...
    roi_x,roi_z+1)';
RF_cartesian_tsi = imbilatfilt(imresize(RF_cartesian_tsi,[512,768]),0.5,1);
JetMap = turbo(64); JetMap = [[0,0,0];JetMap];
f = figure;
imshow(RF_cartesian_tsi,[0.0,0.5],'border','tight'); colormap(JetMap)
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_TSI.tiff');
saveas(gcf,imname)
close(f)
disp('>>>>>>>> TSI CALCULATION COMPLETE.')
f = figure;
imshow((2*RF_cartesian_tsi + RF_cartesian_tai),[0.1,1.25],'border','tight'); colormap(JetMap)
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_TTI.tiff');
saveas(gcf,imname)
close(f)
clearvars naka_map RF_cartesian_tai RF_cartesian_tsi

%% H-SCAN (PARKER 2016 PHYS IN MED BIOL)
[H_ratio,H_mode] = HMODE(R_raw_hist,info);
MPUS.Hratio = H_ratio; MPUS.Hmode = H_mode;
factor = [1.25,0.75,1.25];
for i = 1:3
    RF_cartesian_H_tmp = polar2cartImgPAT_crop(H_mode(:,:,i),[info.ROC*100,info.ROC*100 + max(d_sample)*100],[-1*abs(sc_ang),abs(sc_ang)],roi_x,roi_z+1)';
    RF_cartesian_H(:,:,i) = (imbilatfilt(imresize(RF_cartesian_H_tmp,[512,768]),0.1,1)).*factor(i);
end
RF_cartesian_H = im2uint8(mat2gray(RF_cartesian_H));
HSV = rgb2hsv(RF_cartesian_H);
HSV(:, :, 2) = HSV(:, :, 2)*1.25;
HSV(HSV > 1) = 1;
RF_cartesian_H = hsv2rgb(HSV);
[~,rr] = meshgrid(1:info.Nsc,1:size(H_ratio,1)); tgc = 1./(1+0.5e-3*rr);
H_ratioc = (H_ratio + 1).*tgc;
H_mode_cart = polar2cartImgPAT_crop(H_ratioc,[info.ROC*100,info.ROC*100 + max(d_sample)*100],[-1*abs(sc_ang),abs(sc_ang)],...
    roi_x,roi_z+1)';
H_mode_cart = imbilatfilt(imresize(H_mode_cart,[512,768]),0.25,2);

f = figure;
imshow(RF_cartesian_H,'border','tight'); %colormap(HMap)
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_H.tiff');
saveas(gcf,imname)
close(f)
f = figure;
imshow(H_mode_cart,[0.4,1.25],'border','tight'); colormap('turbo')
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_H_ratio.tiff');
saveas(gcf,imname)
close(f)
disp('>>>>>>>> H MODE CALCULATION COMPLETE.')
clearvars RF_cartesian_H H_ratio H_mode factor

%% SOS (SPEED OF SOUND ESTIMATION)
lambdaG = 1e-4;
lambdaR = 2e-5;
N_iter = 4e4;
[CAVG,CLOC] = SOS_CURVE(Reconstruction_hist, RxMux, info, lambdaG, lambdaR, N_iter);

CAVGc = CAVG;
Ly = 100; sy = 45; y = exp(-((1:Ly)-Ly/2).^2/2/sy/sy); y = y/sum(y);
Ly3 = 30; sy3 = 15; y3 = exp(-((1:Ly3)-Ly3/2).^2/2/sy3/sy3); y3 = y3/sum(y3);
Ly2 = 10; sy2 = 3; y2 = exp(-((1:Ly2)-Ly2/2).^2/2/sy2/sy2); y2 = y2/sum(y2);
CLOC_tmp = padarray(CLOC,[sy,sy2],"replicate","both");
CLOCc = deconvlucy(deconvlucy(deconvlucy(medfilt2(CLOC_tmp),y',5),y3',5),y2',5);
CLOCc = CLOCc(sy+1:end-sy,sy2+1:end-sy2);
Lz = size(CLOCc,1);
CLOCc = imresize(CLOCc,[round(1.4*Lz),info.Nsc],'cubic');
CLOCc = CLOCc(1:Lz,:);

MPUS.CAVG = CAVGc;
MPUS.CLOC = CLOCc;
if info.Nfocus == 2816; zsos = 5.25;
elseif info.Nfocus == 2304; zsos = 4.25;
elseif info.Nfocus < 2000; zsos = 3.25;
elseif info.Nfocus == 4352; zsos = 8.25;
elseif info.Nfocus == 4864; zsos = 8.75;
elseif info.Nfocus == 3328; zsos = 6.25;
else; zsos = 7.25; end
RF_cartesian_Cavg = polar2cartImgPAT_crop(CAVGc,[info.ROC*100 + 0.25,info.ROC*100 + zsos],[-1*abs(sc_ang),abs(sc_ang)],...
    roi_x,roi_z+1)';
RF_cartesian_Cavg = imbilatfilt(imresize(RF_cartesian_Cavg,[512,768]),1,2);
RF_cartesian_Cloc = polar2cartImgPAT_crop(CLOCc,[info.ROC*100 + 0.25,info.ROC*100 + zsos],[-1*abs(sc_ang),abs(sc_ang)],...
    roi_x,roi_z+1)';
RF_cartesian_Cloc = imbilatfilt(imresize(RF_cartesian_Cloc,[512,768]),5,3);

f = figure;
imshow(RF_cartesian_Cavg,[0.0,1.7],'border','tight'); colormap('turbo')
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_CAVG.tiff');
saveas(gcf,imname)
close(f)
f = figure;
imshow(RF_cartesian_Cloc,[0.0,1.7],'border','tight'); colormap('turbo')
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_CLOC.tiff');
saveas(gcf,imname)
close(f)
disp('>>>>>>>> SOS CALCULATION COMPLETE.')

%% SLSC BF
[RF_env_slsc, RF_slsc_log] = SLSC(Reconstruction_hist(:,:,:,1), RxMux, info);
MPUS.SLSC_LOG = RF_slsc_log;

RF_cartesian_SLSC = polar2cartImgPAT_crop(RF_slsc_log,[info.ROC*100,(info.ROC*100 + max(d_sample)*100)],[-1*abs(sc_ang),abs(sc_ang)],...
    roi_x,roi_z+1)';
RF_cartesian_SLSC = imresize(RF_cartesian_SLSC,[512,768]);
f = figure;
imshow(RF_cartesian_SLSC,[0,1.05],'border','tight'); colormap('gray')
set(gca,'xtick','','ytick','')
imname = strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_SLSC_ROI.tiff');
saveas(gcf,imname)
close(f)
disp('>>>>>>>> SLSC CALCULATION COMPLETE.')

%% SAVE RESULTS
save(strcat(savePath,'\',patientNum,'_',areaNum{1,1},'_MPUS_ROI_results.mat'),"MPUS","-v7.3","-nocompression")

end