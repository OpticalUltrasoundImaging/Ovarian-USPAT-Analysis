function processOnePatientMPPA(patientID)
warning('off','all')
addpath('C:\Users\Yixiao\Desktop\Eghbal Yun QPAT simulation\codes');        addpath('C:\Users\Yixiao\Desktop\PAT\PAT recon')
load('C:\Users\Yixiao\Desktop\Eghbal Yun QPAT simulation\codes\Sequence.mat', 'Roi','System');
parentFolder = 'C:\Users\Yixiao\Box\MPUS results v4';
patientResultsFolder = strcat(parentFolder,'\',sprintf('%03d',patientID));
areaFiles = dir(patientResultsFolder);
all_area_dirs = {areaFiles([areaFiles.isdir]).name};
all_area_dirs = all_area_dirs(3:end);
N_scan = length(all_area_dirs);
clearvars areaFiles
fprintf('>>>>>>>> Patient %d has %d scans.\n', patientID,N_scan)

RawDataFolderName = strcat('\\ZHU-LAB-PC\Ovarian_PatientData\p',string(patientID),'\Data');
RawAreaFileList = dir(strcat(RawDataFolderName,'\**\*.*'));
all_raw_data_dirs = {RawAreaFileList([RawAreaFileList.isdir]).name};
all_raw_data_dirs = all_raw_data_dirs(3:end);
clearvars RawAreaFileList
timestampPattern = '^\d{8}-\d{6}$';
tmp = regexp(all_raw_data_dirs,timestampPattern,'match','all');
all_raw_data_dirs = tmp(~cellfun('isempty',tmp));
clearvars timestampPattern tmp
fprintf('>>>>>>>> Found %d has %d RAW data scans.\n', patientID, length(all_raw_data_dirs))

%%
dB_PAT = 40;
min_dB = 10^(-dB_PAT/20);
cm_test = hsv(256); cm_test = [zeros(41,3);flipud(cm_test(1:170,:))];
norm_factor = [17,21,28,27];
f = [0 0.05 0.05 1]; m = [ 0 0 1 1];
DC_cancel = fir2(64,f,m);
f = [0 0.6 0.6 1]; m = [1 1 0 0];
HF_cancel = fir2(64,f,m);
clearvars f m
for idx_scan = 1:N_scan
    time_stamp_tmp = all_area_dirs{idx_scan}; time_stamp_tmp = time_stamp_tmp(5:end);
    raw_data_idx = find(cell2mat(cellfun(@(x) matches(time_stamp_tmp, x), all_raw_data_dirs, 'UniformOutput', false)));
    US_results_tmp = load(strcat(patientResultsFolder,'\',all_area_dirs{idx_scan},'\',all_area_dirs{idx_scan},'_MPUS_ROI_results'));
    US_results_tmp = US_results_tmp.MPUS;
    nfocus_tmp = size(US_results_tmp.RF_LOG,1);
    info = loadCSystemParam(Roi,System);
    info.N_ch = 64;
    R_raw_hist = [];
    R_env_raw_hist = [];
    for i_lambda = 1:4
        idx_area_tmp = raw_data_idx + i_lambda - 1;
        area_id = all_raw_data_dirs(idx_area_tmp);
        area_id = area_id{1,1};
        pname = string(fullfile(RawDataFolderName,'\',area_id,'\'));
        if ~exist(pname,"dir")
            pname = string(fullfile(RawDataFolderName,'\selected\',area_id,'\'));
        end
        ChannelFile = [num2str(1),'_layer0_idx',num2str(1),'_CUSTOMDATA_RF.mat'];
        ChannelFile =  fullfile(pname,ChannelFile);
        if ~exist(ChannelFile,'file')
           ChannelFile = [num2str(2),'_layer0_idx',num2str(2),'_CUSTOMDATA_RF.mat'];
           ChannelFile =  fullfile(pname,ChannelFile);
        end
        ChannelData = load (ChannelFile);
        clearvars pname ChannelFile idx_tmp
        [RF_Sum, RF_env] = PAT_das_curve(ChannelData, info);
        R_raw_hist(:,:,i_lambda) = RF_Sum*1e3/norm_factor(i_lambda);
        R_env_raw_hist(:,:,i_lambda) = RF_env*1e3/norm_factor(i_lambda);
        fprintf('%d / %d.\n',i_lambda,4)
    end
    info.Nfocus = size(R_env_raw_hist,1);
    info.d_sample = (0:info.Nfocus-1) * info.pixel_d;
    R_env_clean = mean(R_env_raw_hist,3);
    [~,ZZ] = meshgrid(1:info.Nsc,1:info.Nfocus);
    R_env_clean = R_env_clean.*(1+1e-7*ZZ);
    R_env_clean = R_env_clean/max(max(R_env_clean));
    idx = find(R_env_clean < min_dB);
    RF_log = (20/dB_PAT)*log10(R_env_clean)+1;
    RF_log(idx) = 0;
    clearvars idx ZZ
    sc_ang = info.ScanAngle(1)*pi/180;
    Nt = size(R_env_raw_hist,1);
    convert_factor = 0.8;
    RF_cartesian = polar2cartImgPAT(RF_log(1:round(Nt*convert_factor),:),[info.ROC*100,info.ROC*100 + max(info.d_sample)*convert_factor*100],[-1*abs(sc_ang),abs(sc_ang)])';
    savePath = strcat(patientResultsFolder,'\',all_area_dirs{idx_scan});
    f = figure;
    imshow(RF_cartesian,[0.0,0.9],'border','tight'); colormap('hot')
    set(gca,'xtick','','ytick','')
    imname = strcat(savePath,'\',all_area_dirs{idx_scan},'_PAT.tiff');
    saveas(gcf,imname)
    close(f)
    MPPA.PA_rfdata = R_raw_hist;
    MPPA.PA_envdata = R_env_raw_hist;
    MPPA.PA = mean(R_env_raw_hist,3);
    disp('BEAMFORMING COMPLETE.')

    %% CALCULATE CENTER FREQUENCY
    N_repeat = 3;
    l_pad = 2^9;
    freq_axis = linspace(0,20,l_pad/2)';
    pulse_length = round((1.5*info.c/info.fc)/info.pixel_d);
    R_slope_list = zeros(size(R_raw_hist,1),size(R_raw_hist,2),4);
    R_intercept_list = zeros(size(R_raw_hist,1),size(R_raw_hist,2),4);
    fwb = waitbar(0,'Solving for frequency spectrum ... ');
    for i_wl = 1:size(R_raw_hist,3)
        %imagesc(R_raw_hist(:,:,1),[-2e2,2e2])
        for i_sc = 1:info.Nsc
            rf_tmp = R_raw_hist(:,i_sc,i_wl);
            rf_tmp = conv(conv(rf_tmp, DC_cancel, 'same'),HF_cancel,'same');
            freq_slope_aline = zeros(info.Nfocus,N_repeat); freq_intercept_aline = zeros(info.Nfocus,N_repeat);
            for i_z = 1:info.Nfocus
                for i_rpt = 1:N_repeat
                    h_roi = ((i_rpt)*2+5)*pulse_length;
                    i_z_start = max(1,i_z-floor(h_roi/2));
                    i_z_end = min(i_z+floor(h_roi/2),info.Nfocus);
                    roi_tmp = rf_tmp(i_z_start:i_z_end);
                    h_roi_tmp = length(roi_tmp);
                    roi_tmp_pad = [roi_tmp;zeros(l_pad-h_roi_tmp,1)];
                    roi_tmp_fft = abs(fftshift(fft(roi_tmp_pad,[],1),1));
                    roi_tmp_fft = flip(roi_tmp_fft(1:l_pad/2))+roi_tmp_fft(l_pad/2+1:l_pad);
                    roi_tmp_fft = roi_tmp_fft/max(roi_tmp_fft);
                    X = [ones(l_pad/2,1) freq_axis];
                    b = X(14:204,:)\roi_tmp_fft(14:204);
                    freq_slope_aline(i_z,i_rpt) = min(0,b(2)); freq_intercept_aline(i_z,i_rpt) = max(0,b(1));
                end
            end
            R_slope_list(:,i_sc,i_wl) = mean(freq_slope_aline,2); R_intercept_list(:,i_sc,i_wl) = mean(freq_intercept_aline,2);
            perc = (i_sc+(i_wl-1)*info.Nsc)/info.Nsc/4;
            msg = sprintf('Solving for frequency spectrum ... %.3f',perc);
            waitbar(perc,fwb,msg);
        end
    end
    close(fwb)
    R_slope_list = mean(R_slope_list,3);
    R_intercept_list = mean(R_intercept_list,3);
    MPPA.Fslope = R_slope_list;
    MPPA.Fintercept = R_intercept_list;
    disp('CENTER FREQUENCY CALCULATION COMPLETE.')

    %% CALCULATE SO2
    RF_env_raw_hist_cart = [];
    convert_factor = 1;
    for idx_tmp = 1:4
        RF_cartesian = polar2cartImgPAT(R_env_raw_hist(1:round(Nt*convert_factor),:,idx_tmp),[info.ROC*100,info.ROC*100 + max(info.d_sample)*convert_factor*100],[-1*abs(sc_ang),abs(sc_ang)])';
        RF_env_raw_hist_cart(:,:,idx_tmp) = abs(RF_cartesian);
    end
    ME_HBO2 = [390, 720, 800, 980];
    ME_HBr = [1100, 1120, 920, 840];
    ME_other = [0.0203, 0.0267, 0.0228, 0.0316]/55.56;
    
    ds_ratio = 2;
    Lz = round(size(RF_env_raw_hist_cart,1)/ds_ratio); Lx = round(size(RF_env_raw_hist_cart,2)/ds_ratio);
    R_env_hist = imresize(RF_env_raw_hist_cart,[Lz,Lx],'bilinear');
    HBR_img = zeros(Lz,Lx);
    HBO2_img = zeros(Lz,Lx);
    so2_mask = ones(Lz,Lx);
    blur = 2;
    fwb = waitbar(0,'Solving for Hb conc ... ');
    for j = 1:Lx
        x_range = (j-blur):(j+blur);
        x_range(x_range<1)=1; x_range(x_range>Lx)=Lx;
        for idx_z = 1:Lz
            z_range = (idx_z-blur):(idx_z+blur);
            z_range(z_range<1)=1; z_range(z_range>Lz)=Lz;
            data_tmp = R_env_hist(z_range,x_range,:);
            % check similarity between wavelengths
            data_check = reshape(data_tmp,[],4);
            R_tmp = abs(corrcoef(data_check));
            R_tmp = (sum(R_tmp,1)-1)/3;
            idx = find(R_tmp<0.33);
            if length(idx)<2 || isempty(idx)
                lambda_used = [1,2,3,4]; lambda_used(idx)=[];
                data_tmp(:,:,idx) = [];
                data_tmp = permute(data_tmp,[3,1,2]);
                data_tmp = data_tmp(:);
                abs_map_d = [ME_HBO2(lambda_used); ME_HBr(lambda_used); ME_other(lambda_used)]';
                Aabs_d = repmat(abs_map_d,(2*blur+1)^2,1);
                c_tmp = lsqnonneg(Aabs_d,data_tmp);
                HBR_img(idx_z,j) = c_tmp(2);
                HBO2_img(idx_z,j) = c_tmp(1);
            else
                so2_mask(idx_z,j) = 0;
            end
        end
        perc = j/Lx;
        msg = sprintf('Solving for Hb conc ... %.3f',perc);
        waitbar(perc,fwb,msg);
    end
    close(fwb)
    THB_img = HBO2_img + HBR_img;
    MPPA.THb = THB_img;
    SO2_img = medfilt2(HBO2_img)./medfilt2(THB_img);
    SO2_img(isnan(SO2_img)) = 0; SO2_img(isinf(SO2_img)) = 0; SO2_img(SO2_img>1)=0;
    MPPA.SO2 = SO2_img;
    MPPA.SO2MASK = so2_mask;

    mask2 = mat2gray(THB_img>0.1);
    SO2_img = SO2_img.*so2_mask.*mask2;

    f = figure;
    imshow(SO2_img,[0.0,1.0],'border','tight'); colormap(cm_test)
    set(gca,'xtick','','ytick','')
    imname = strcat(savePath,'\',all_area_dirs{idx_scan},'_SO2.tiff');
    saveas(gcf,imname)
    close(f)
    disp('SO2 CALCULATION COMPLETE.')
    %% SAVE RESULTS
    MPPA.Nfocus = info.Nfocus;
    save(strcat(savePath,'\',all_area_dirs{idx_scan},'_MPPA_ROI_results.mat'),"MPPA","-v7.3","-nocompression")
    clc
    fprintf('>>>>>>>> Patient %d: %d / %d scans processed.\n', patientID, idx_scan, N_scan)
end
end