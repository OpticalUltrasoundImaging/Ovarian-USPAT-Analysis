function [RF_beam, RF_Sum] = US_das_curve(ChannelData, info)
%% REORDER DATA INTO SCANLINES
info.Nfocus = double(ChannelData.AlignedSampleNum);
info.d_sample = 0:info.pixel_d:(info.pixel_d*(info.Nfocus-1));
Reconstruction = zeros(ChannelData.AlignedSampleNum, info.N_ch, info.Nsc);
for sc = 1:info.Nsc
    ti = sprintf('AdcData_scanline%03d_thi0_sa0',sc-1);
    imgtmp = double(ChannelData.(ti));
    %img = [imgtmp(:,:,1);imgtmp(:,:,2)];
    img = reorder_chdata_ch2dly(sc-1, info.N_ele, info.Nsc, info.N_ch, imgtmp);
    Reconstruction(:,:,sc) = img';
end
clearvars imgtmp img fwd
%% REORDER INFORMATION
data_total = floor(ChannelData.AlignedSampleNum);
f = [0 0.1 0.1 1]; m = [ 0 0 1 1];
DC_cancle = fir2(64,f,m);
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
    RF_beam(:,sc) = sum_tmp;                            RF_Sum(:,sc) = abs(hilbert(sum_tmp.*CF));
end
close(fwd)
clearvars fwd msg perc
clearvars idx1 idx2 idx3 idx4 idx_left idx_right h_aper_size fil_tmp ChData CF abs_rf active_ch
clearvars rf_tmp row_abs_sum row_sum_abs Rx_apo_r Rx_apod_idx rx_ch_mtx rx_idx rx_idx_tmp rx_t RxReadPointer
clearvars x z Xr t ti start_right start_left
clearvars tmp TOF tx_t sum_tmp
% %% dynamic ranging
% min_dB = 10^(-dB_PAT/20);
% RF_env = abs(RF_Sum);
% RF_env_raw = RF_env;
% RF_env_norm = RF_env./max(max(RF_env));
% idx = find(RF_env_norm < min_dB);
% RF_log = (20/dB_PAT)*log10(RF_env_norm)+1;
% RF_log(idx) = 0;
%
% Nt = size(RF_log,1);
% convert_factor = (5e-2)/max(info.d_sample);
% RF_cartesian = polar2cartImgPAT(RF_log(1:round(Nt*convert_factor),:),[info.ROC*100,info.ROC*100 + max(info.d_sample)*convert_factor*100],[-1*abs(sc_ang),abs(sc_ang)])';
% figure
% imshow(RF_cartesian,[0.02,1])
% colormap('hot')
% set(gca,'xtick','','ytick','','fontsize',14,'fontweight','bold')
% colorbar()
% title(sprintf('Recon DAS, %ddB',(dB_PAT)))
end
