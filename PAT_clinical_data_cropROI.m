clear
clc
warning('off','all')
patientID = 49;
parentFolder = 'C:\Users\Yixiao Lin\Box\MPUS results v4';
patientResultsFolder = strcat(parentFolder,'\',sprintf('%03d',patientID));
ROIResultsFolder = 'C:\Users\Yixiao Lin\Desktop\Eghbal Yun QPAT simulation\MPUS_ROI';
areaFiles = dir(patientResultsFolder);  all_area_dirs = {areaFiles([areaFiles.isdir]).name};    all_area_dirs = all_area_dirs(3:end);
N_scan = length(all_area_dirs);
clearvars areaFiles
fprintf('>>>>>>>> Patient %d has %d scans.\n', patientID,N_scan)
sos_zstart = 0.225;
sos_lookup = [1792, 3.25; 2304, 4.25; 2816, 5.25; 3328, 6.25; 3840, 7.25; 4864, 8.75; 5376, 9.5];
US_pixel_d = 1540/2/40e4; % [cm]
PA_pixel_d = 1540/40e4;   % [cm]
dB_PAT = 35;        min_dB = 10^(-dB_PAT/20);
dB_PAT_roi = 25;    min_dB_roi = 10^(-dB_PAT_roi/20);
sc_ang = 150.7/2*pi/180;
%%
ROI_stat_all = zeros(N_scan,116);
for idx_scan = 1:N_scan
    area_id_tmp = all_area_dirs{idx_scan};
    savePath = fullfile(ROIResultsFolder,area_id_tmp);
    if exist(savePath,'file')
        disp('>>>>>>>> An ROI from this image already exists')
    end
    timestamp_tmp = str2double(area_id_tmp(end-5:end));
    US_results_tmp = load(strcat(patientResultsFolder,'\',area_id_tmp,'\',area_id_tmp,'_MPUS_ROI_results'));
    US_results_tmp = US_results_tmp.MPUS;
    PA_results_tmp = load(strcat(patientResultsFolder,'\',area_id_tmp,'\',area_id_tmp,'_MPPA_ROI_results'));
    PA_results_tmp = PA_results_tmp.MPPA;
    %% CROP ROI FROM ULTRASOUND SCAN
    L_US = US_pixel_d*size(US_results_tmp.RF_LOG,1);
    L_PA = PA_pixel_d*size(PA_results_tmp.PA,1);
    convert_factor_US = L_PA / max(L_US,L_PA);
    convert_factor_PA = L_US / max(L_US,L_PA);
    us_polar = US_results_tmp.RF_LOG; us_polar = us_polar(1:round(size(us_polar,1)*convert_factor_US),:);
    pa_polar = PA_results_tmp.PA;     pa_polar = pa_polar(1:round(size(pa_polar,1)*convert_factor_PA),:);  pa_polar(1:265,:)=0;     pa_raw_polar = pa_polar;
    [zz,xx] = meshgrid(1:size(pa_polar,1),1:size(pa_polar,2));
    tgc = 1+0.165e-2*zz;
    tgc2 = 1+0.25*(64-abs(xx-64));
    pa_polar = tgc2'.*tgc'.*pa_polar;
    pa_polar(1501:end,:)=0;
    pa_polar = pa_polar/max(pa_polar(:));
    imagesc(pa_polar)
    dB_PAT = 25;        min_dB = 10^(-dB_PAT/20);
    idx = find(pa_polar < min_dB);
    pa_polar = (20/dB_PAT)*log10(pa_polar)+1;
    pa_polar(idx) = 0;
    clearvars idx
    figure;imagesc(pa_polar)
    us_cartesian = polar2cartImgPAT(us_polar,[1,1 + L_US],[-1*sc_ang,sc_ang])'; us_cartesian = imresize(us_cartesian,[1024,2048]);
    pa_cartesian = polar2cartImgPAT(pa_polar,[0.9,0.9 + L_US],[-1*sc_ang,sc_ang])'; pa_cartesian = imresize(pa_cartesian,[1024,2048]);
    paRAW_cartesian = polar2cartImgPAT(pa_raw_polar,[0.9,0.9 + L_US],[-1*sc_ang,sc_ang])'; paRAW_cartesian = imresize(paRAW_cartesian,[1024,2048]);
    pa_cartesian_hot = double2hot((imbilatfilt(pa_cartesian,0.5,3)),0.85,0.025,'red');
    figure;imagesc(pa_cartesian_hot)
    %% FREEHAND CONTOUR SEGMENTATION
    f = figure('Position',[100,400,1000,400]);
    imshow(imsharpen(imbilatfilt(us_cartesian,0.75,1)),[0,0.9], 'border', 'tight')
    colormap('gray')
    set(gca,'xtick','','ytick','','fontsize',14)
    ovary_contour = drawfreehand('Color','r','FaceAlpha',0);
    d = round(ovary_contour.Position);
    saveas(gcf,'C:\Users\Yixiao\Desktop\Eghbal Yun QPAT simulation\MPUSPAT_results\eg_image.tif')
    close(f)
    [Psi,~] = computeShapeIndices(d,L_US);  psi_cell = (struct2cell(Psi));
    Zmoment = computeShapeMoments(d);
    shape_array = zeros(1,13);
    for i = 1:10
        shape_array(i) = psi_cell{i};
    end
    shape_array(11:13) = Zmoment;
    ROI.US_SHAPE = shape_array;
    clearvars psi_cell Zmoment

    %% FOR US PROCESSING
    f = figure('Position',[300,400,1400,400]);
    ax(1) = subplottight(1,2,1); imshow(us_cartesian,[0,0.9], 'border', 'tight');  colormap(ax(1),'gray'); set(gca,'xtick','','ytick','','fontsize',14); title('US')
    ax(2) = subplottight(1,2,2);
    us = imshow(us_cartesian,[0,0.95], 'border', 'tight');    colormap(ax(2),'gray');   hold on;     
    pa = imshow(pa_cartesian_hot, 'border', 'tight');   hold off;        alpha(pa,min(1,pa_cartesian*2.25)); set(gca,'xtick','','ytick','','fontsize',14); title('US + PA')
    roi_boundary = ginput(2); % THE FIRST POINT IS TOP LEFT, THE SECOND POINT IS BOTTOM RIGHT!
    waitfor(roi_boundary)
    roi_boundary = round(roi_boundary); roi_boundary(1,:) = roi_boundary(1,:) - 5; roi_boundary(2,:) = roi_boundary(2,:) + 5;
    % CONVERT COORDINATES TO PHYSICAL DIMENSIONS
    z_range = roi_boundary(:,2)/1024*(1+L_US); % [cm]
    x_range = (roi_boundary(:,1)-1024)/1024*(1+L_US); % [cm]
    fprintf('US ROI #%d coordinates are z: [%.2f,%.2f], x: [%.2f,%.2f]\n',idx_scan, z_range(1),z_range(2),x_range(1),x_range(2))
    saveas(gcf,'C:\Users\Yixiao\Desktop\Eghbal Yun QPAT simulation\MPUSPAT_results\eg_NORMAL_image.tif')
    close(f)

    ROI.UScoordinates = roi_boundary;
    ROI.USdimensionsCM  = [x_range, z_range];
    ROI.US = us_cartesian(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    Nlevel = 32;
    [GLCM , ~] = ComputeNGLDM(ROI.US, 2, 0, Nlevel); glcm_array = cell2mat(struct2cell(GLCM));
    [GLSDM, ~] = ComputeGLCM(ROI.US,Nlevel);         glsdm_array = cell2mat(struct2cell(GLSDM));
    [GLSZM, ~] = ComputeGLSZM(ROI.US,Nlevel);        glszm_array = cell2mat(struct2cell(GLSZM));
    ROI.US_RADIOMICS = [glcm_array' , glsdm_array' , glszm_array'];
    clearvars GLCM GLSDM GLSZM glcm_array glsdm_array glszm_array
    
    %% EXTRACT MPUS ROI (FC, TAI, TSI, HMODE, CAVG, SLSC)
    ROI_MPUS = zeros(1,12);
    us_fc = US_results_tmp.FC; us_fc = us_fc(1:round(size(us_fc,1)*convert_factor_US),:);
    us_fc = polar2cartImgPAT(us_fc,[0.9,0.9 + L_US],[-1*sc_ang,sc_ang])'; us_fc = imresize(us_fc,[1024,2048]);
    ROI_US_Fc = us_fc(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    ROI_MPUS(1) = mean(ROI_US_Fc(:)); ROI_MPUS(2) = std(ROI_US_Fc(:));
    clearvars us_fc ROI_US_fc
    us_tai = US_results_tmp.TAI; us_tai = us_tai(1:round(size(us_tai,1)*convert_factor_US),:);
    us_tai = polar2cartImgPAT(us_tai,[0.9,0.9 + L_US],[-1*sc_ang,sc_ang])'; us_tai = imresize(us_tai,[1024,2048]);
    ROI_US_TAI = us_tai(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPUS(3), ROI_MPUS(4)] = FindFWHM(ROI_US_TAI);
    clearvars us_tai ROI_US_TAI
    us_tsi = US_results_tmp.TSI; us_tsi = us_tsi(1:round(size(us_tsi,1)*convert_factor_US),:);
    us_tsi = polar2cartImgPAT(us_tsi,[0.9,0.9 + L_US],[-1*sc_ang,sc_ang])'; us_tsi = imresize(us_tsi,[1024,2048]);
    us_tsi(isinf(us_tsi)) = 0; us_tsi(isnan(us_tsi)) = 0;
    ROI_US_TSI = us_tsi(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPUS(5), ROI_MPUS(6)] = FindFWHM(ROI_US_TSI);
    clearvars us_tsi ROI_US_TSI
    us_hmode = US_results_tmp.Hmode; us_hmode = us_hmode(1:round(size(us_hmode,1)*convert_factor_US),:,:);
    us_hmode = (us_hmode(:,:,1)-us_hmode(:,:,3))./(us_hmode(:,:,1)+us_hmode(:,:,3));
    us_hmode = polar2cartImgPAT(us_hmode,[0.9,0.9 + L_US],[-1*sc_ang,sc_ang])'; us_hmode = imresize(us_hmode,[1024,2048]);
    us_hmode(isinf(us_hmode)) = 0; us_hmode(isnan(us_hmode)) = 0;
    ROI_US_HMODE = us_hmode(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPUS(7), ROI_MPUS(8)] = FindFWHM(ROI_US_HMODE);
    clearvars us_hmode ROI_US_HMODE
    us_cavg = US_results_tmp.CAVG;
    L_US_pix = size(US_results_tmp.RF_LOG,1);
    L_sos = 5.25;
    for i = 1:size(sos_lookup,1)
        if L_US_pix == sos_lookup(i,1)
            L_sos = sos_lookup(i,2);
        end
    end
    US_dr = 1540/2/7.1e4;
    prepad = round(sos_zstart/US_dr);
    postpad = round((L_US - L_sos)/US_dr);
    us_cavg_pad = [zeros(prepad,128); us_cavg; zeros(postpad,128)];
    us_cavg = polar2cartImgPAT(us_cavg_pad,[0.9 , 0.9 + L_US],[-1*sc_ang,sc_ang])'; us_cavg = imresize(us_cavg,[1024,2048]);
    ROI_US_CAVG = us_cavg(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPUS(9), ROI_MPUS(10)] = FindFWHM(ROI_US_CAVG);
    clearvars us_cavg ROI_US_CAVG
    us_slsc = US_results_tmp.SLSC_LOG; us_slsc = us_slsc(1:round(size(us_slsc,1)*convert_factor_US),:);
    us_slsc = polar2cartImgPAT(us_slsc,[0.9,0.9 + L_US],[-1*sc_ang,sc_ang])'; us_slsc = imresize(us_slsc,[1024,2048]);
    ROI_US_SLSC = us_slsc(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPUS(11), ROI_MPUS(12)] = FindFWHM(ROI_US_SLSC);
    clearvars us_slsc ROI_US_SLSC
    ROI.USMP = ROI_MPUS;
    clearvars ROI_MPUS
    
    img_plot = ROI_US_CAVG;
    img_plot = (medfilt2(imbilatfilt((img_plot),0.75,3),[5,5]));
    [zz,~] = meshgrid(1:size(img_plot,1),1:size(img_plot,2));
    tgc = 1+0.165e-2*zz;
    img_plot = tgc'.*img_plot;
    img_plot_hot = double2hot((img_plot),1.125,0.1,'turbo'); % SOS:[1.15,0.1], TAI:[0.65,0.03], TSI:[0.5,0.03]
    figure;imagesc(img_plot_hot)
    figure;
    us = imshow((ROI.US),[0,2.5], 'border', 'tight');    colormap('gray');   hold on;
    pa = imshow(img_plot_hot, 'border', 'tight');   hold off;        alpha(pa,min(1,0.0+img_plot*1.25));
    set(gca,'xtick','','ytick','','fontweight','bold','fontsize',12)
    %colorbar()
    
    %% FOR PA PROCESSING
    f = figure('Position',[300,400,1400,400]);
    ax(1) = subplottight(1,2,1); imshow(us_cartesian,[0,0.9], 'border', 'tight');  colormap(ax(1),'gray'); set(gca,'xtick','','ytick','','fontsize',14); title('  ')
    ax(2) = subplottight(1,2,2);
    us = imshow(us_cartesian,[0,0.95], 'border', 'tight');    colormap(ax(2),'gray');   hold on;     
    pa = imshow(pa_cartesian_hot, 'border', 'tight');   hold off;        alpha(pa,min(1,pa_cartesian*2.25)); set(gca,'xtick','','ytick','','fontsize',14); title('COREG')
    roi_boundary = ginput(2); % THE FIRST POINT IS TOP LEFT, THE SECOND POINT IS BOTTOM RIGHT!
    waitfor(roi_boundary)
    roi_boundary = round(roi_boundary); roi_boundary(1,:) = roi_boundary(1,:) - 5; roi_boundary(2,:) = roi_boundary(2,:) + 5;
    % CONVERT COORDINATES TO PHYSICAL DIMENSIONS
    z_range = roi_boundary(:,2)/1024*(1+L_US); % [cm]
    x_range = (roi_boundary(:,1)-1024)/1024*(1+L_US); % [cm]
    fprintf('PA ROI #%d coordinates are z: [%.2f,%.2f], x: [%.2f,%.2f]\n',idx_scan, z_range(1),z_range(2),x_range(1),x_range(2))
    close(f)
    ROI.PAcoordinates = roi_boundary;
    ROI.PAdimensionsCM  = [x_range, z_range];

    pa_roi_tmp = paRAW_cartesian(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [zz,xx] = meshgrid(1:size(pa_roi_tmp,1),1:size(pa_roi_tmp,2));
    tgc = 1+0.165e-2*zz;
    pa_roi_tmp = tgc'.*pa_roi_tmp;
    pa_roi_tmp = pa_roi_tmp/max(pa_roi_tmp(:));
    dB_PAT_roi = 25;    min_dB_roi = 10^(-dB_PAT_roi/20);
    idx = find(pa_roi_tmp < min_dB_roi);
    pa_roi_tmp = (20/dB_PAT_roi)*log10(pa_roi_tmp)+1;
    pa_roi_tmp(idx) = 0;
    ROI.PA = pa_roi_tmp;

    figure;imagesc(pa_roi_tmp)
    clearvars idx pa_roi_tmp
    [GLCM , ~] = ComputeNGLDM(ROI.PA, 2, 0, Nlevel); glcm_array = cell2mat(struct2cell(GLCM));
    [GLSDM, ~] = ComputeGLCM(ROI.PA,Nlevel);         glsdm_array = cell2mat(struct2cell(GLSDM));
    [GLSZM, ~] = ComputeGLSZM(ROI.PA,Nlevel);        glszm_array = cell2mat(struct2cell(GLSZM));
    ROI.PA_RADIOMICS = [glcm_array' , glsdm_array' , glszm_array'];
    clearvars GLCM GLSDM GLSZM glcm_array glsdm_array glszm_array

    %% EXTRACT MPPA ROI (FSLOPE, FINTERCEPT, THB, SO2)
    ROI_MPPA = zeros(1,8);
    pa_fslope = PA_results_tmp.Fslope; pa_fslope = pa_fslope(1:round(size(pa_fslope,1)*convert_factor_PA),:);
    pa_fslope = polar2cartImgPAT(pa_fslope,[1,1 + L_PA],[-1*sc_ang,sc_ang])'; pa_fslope = imresize(pa_fslope,[1024,2048]);
    ROI_PA_Fslope = pa_fslope(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPPA(1), ROI_MPPA(2)] = FindFWHM(ROI_PA_Fslope);
    clearvars pa_fslope ROI_PA_Fslope
    pa_fintercept = PA_results_tmp.Fintercept; pa_fintercept = pa_fintercept(1:round(size(pa_fintercept,1)*convert_factor_PA),:);
    pa_fintercept = polar2cartImgPAT(pa_fintercept,[1,1 + L_PA],[-1*sc_ang,sc_ang])'; pa_fintercept = imresize(pa_fintercept,[1024,2048]);
    ROI_PA_Fintercept = pa_fintercept(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPPA(3), ROI_MPPA(4)] = FindFWHM(ROI_PA_Fintercept);
    clearvars pa_fintercept ROI_PA_Fintercept
    pa_thb = PA_results_tmp.THb; H_pa = size(pa_thb,1); L_pa = size(pa_thb,2);
    crop_pa = round(H_pa*(1-convert_factor_PA*1.125)); pa_thb = pa_thb(1:round(H_pa*convert_factor_PA*1.125), (crop_pa+1):(L_pa-crop_pa)); pa_thb = imresize(pa_thb,[1024,2048],"nearest");
    pa_thb = imbilatfilt(medfilt2(pa_thb),0.5,1); pa_thb(pa_thb<0) = 0;
    ROI_PA_THB = pa_thb(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    [ROI_MPPA(5), ROI_MPPA(6)] = FindFWHM(ROI_PA_THB);
    clearvars pa_thb
    pa_so2 = PA_results_tmp.SO2;
    pa_so2 = pa_so2(1:round(H_pa*convert_factor_PA*1.125), (crop_pa+1):(L_pa-crop_pa)); pa_so2 = imresize(pa_so2,[1024,2048],"nearest");
    pa_so2 = imbilatfilt(medfilt2(pa_so2),0.5,1); pa_so2(pa_so2<0) = 0;
    ROI_PA_SO2 = pa_so2(roi_boundary(1,2):roi_boundary(2,2) , roi_boundary(1,1):roi_boundary(2,1));
    ROI_PA_SO2MASK = imfill(imclose(imdilate(ROI_PA_THB>1.25,strel('disk',2)),strel('disk',5)),'holes');

    ROI_PA_SO2MASK2 = (ROI_PA_SO2>0.2).*(ROI_PA_SO2<1.01);
    ROI_PA_SO2MASK = ROI_PA_SO2MASK.*ROI_PA_SO2MASK2;
    so2_data = ROI_PA_SO2(ROI_PA_SO2MASK==1);
    so2_data(so2_data == 0) = NaN;
    [ROI_MPPA(7), ROI_MPPA(8)] = FindFWHM(so2_data);
    ROI.PAMP = ROI_MPPA;
    clearvars pa_so2 H_pa L_pa ROI_PA_SO2MASK ROI_PA_SO2MASK2 ROI_PA_THB ROI_MPPA
    ROI_stat_all(idx_scan,1) = timestamp_tmp;
    ROI_stat_all(idx_scan,2:14) = ROI.US_SHAPE;
    ROI_stat_all(idx_scan,15:26) = ROI.USMP;
    ROI_stat_all(idx_scan,27:34) = ROI.PAMP;
    ROI_stat_all(idx_scan,35:75) = ROI.US_RADIOMICS;
    ROI_stat_all(idx_scan,76:116) = ROI.PA_RADIOMICS;

    %% SAVE DATA
    save(strcat(savePath,'.mat'),"ROI","-v7.3","-nocompression")
end