function [RF_env_slsc, RF_slsc_log] = SLSC(Reconstruction, RxMux, info)
lambda = info.c/info.fc;
DeltaN = round(1.5*lambda/info.pixel_d);
flp = 10;
ftmp = 2*flp*1e6/info.fs;
f = [0 ftmp ftmp 1];
mlo = [1 1 0 0];
lpf = fir2(64,f,mlo);
f = [0 0.03 0.03 1]; m = [0 0 1 1];
DC_cancel = fir2(64,f,m);
RF_SLSC_sum = zeros(info.Nfocus,info.Nsc);
fwd = waitbar(0,'SLSC scanline ... ');
for sc = 1:info.Nsc
    perc = sc/info.Nsc;
    msg = sprintf('SLSC scanline # %0.3f',perc);
    waitbar(perc,fwd,msg);

    ChData = Reconstruction(:,:,sc);
    MuxTbl = RxMux(sc,:);
    sc_ang = info.ScanAngle(sc)*pi/180;

    RF_scanline_tofadjusted = zeros(info.Nfocus,info.N_ch);
    active_ch = [];
    for ch = 1:info.N_ch
        ElePos = MuxTbl(ch);
        if(ElePos > 0 && ElePos <= info.N_ele)
            active_ch = [active_ch,ch];
            tmp = ChData(:,ch);
            fil_tmp = conv(conv(tmp, DC_cancel, 'same'),lpf,'same');

            % CALCULATE DELAY
            ele_ang = info.EleAngle(ElePos)*pi/180;
            ch_ang = abs(sc_ang - ele_ang);

            x = info.ROC*sin(ch_ang);
            z = info.ROC*(1 - cos(ch_ang));

            tx_t= info.d_sample/info.c;
            rx_t = sqrt(x.^2 + (info.d_sample + z).^2)/info.c;
            TOF = tx_t + rx_t;
            RxReadPointer = round(TOF*info.fs);

            idx1 = find(RxReadPointer > info.Nfocus);    RxReadPointer(idx1) = info.Nfocus;
            idx2 = find(RxReadPointer < 1);             RxReadPointer(idx2) = 1;

            rf_tmp = fil_tmp(RxReadPointer);
            rf_tmp(idx1) = 0;                           rf_tmp(idx2) = 0;
            RF_scanline_tofadjusted(:,ch) = rf_tmp';
        end
    end

    % COMPUTE R(m), ASSUME SC = ELEMENT IDX
    %imagesc(RF_scanline_tofadjusted)
    Rm = zeros(info.Nfocus,round(info.N_ch/2));
    for idx_z = 1:info.Nfocus
        idx_z_start = max(1,idx_z - floor(DeltaN/2));
        idx_z_end = min(info.Nfocus,idx_z + floor(DeltaN/2));
        line_tmp = RF_scanline_tofadjusted(idx_z_start:idx_z_end,active_ch);
        for idx_m = 0:round(info.N_ch/2)-1
            corr_sum = 0;
            for idx_i = 1:(length(active_ch) - idx_m)
                signali = line_tmp(:,idx_i);
                signalim = line_tmp(:,idx_i + idx_m);
                num = sum(signali.*signalim);
                denom = sqrt(sum(signali.^2)*sum(signalim.^2));
                corr_sum = corr_sum + num/denom;
            end
            Rm(idx_z,idx_m+1) = 1/(info.N_ch - idx_m)*corr_sum;
        end
    end
    sum_slsc = sum(Rm(:,2:11),2)/length(active_ch); % MAY NEED TO CHANGE SUM NUMBER
    RF_SLSC_sum(:,sc) = sum_slsc;
end
close(fwd)

%% DYNAMIC RANGING AND LOG COMPRESSION
dB_US = 13;
min_dB = 10^(-dB_US/20);
RF_env_slsc = abs(RF_SLSC_sum);
RF_env_slsc_norm = RF_env_slsc./max(max(RF_env_slsc));
idx = find(RF_env_slsc_norm < min_dB);
RF_slsc_log = (20/dB_US)*log10(RF_env_slsc_norm)+1;
RF_slsc_log(idx) = 0;
%RF_slsc_log = imsharpen(RF_slsc_log);
clearvars idx
end