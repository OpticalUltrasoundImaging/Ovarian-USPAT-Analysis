function naka_avg = TSI_CURVE(R_env_raw_hist,info)
pulse_length = (1.5*info.c/info.fc);
N_repeat = 7; N_frame = size(R_env_raw_hist,3);
naka_avg = zeros(size(R_env_raw_hist,1),size(R_env_raw_hist,2),N_repeat*N_frame);
fwb = waitbar(0,'CALCULATING NAKAGAMI PARAMETER ... ');
for i_frame = 1:N_frame
    R_env_raw = R_env_raw_hist(:,:,i_frame);
    for i_repeat = 1:N_repeat
        w = min(1,floor(i_repeat/3))*0.5e-3;    h = (i_repeat+2)*pulse_length;
        w_roi = round(w/info.pitch);        if mod(w_roi,2)==0; w_roi = w_roi+1; end
        h_roi = round(h/info.pixel_d);      if mod(h_roi,2)==0; h_roi = h_roi+1; end
        naka_map = zeros(size(R_env_raw));
        for i_x = 1:info.Nsc
            i_x_start = max(1,i_x-floor(w_roi/2));
            i_x_end = min(i_x+floor(w_roi/2),info.Nsc);
            parfor i_z = 1:info.Nfocus
                i_z_start = max(1,i_z-floor(h_roi/2));
                i_z_end = min(i_z+floor(h_roi/2),info.Nfocus);
                roi_tmp = R_env_raw(i_z_start:i_z_end,i_x_start:i_x_end);
                roi_tmp = roi_tmp(:);
                m2_tmp = moment(roi_tmp(:),2);
                m4_tmp = moment(roi_tmp(:),4);
                naka_map(i_z,i_x) = m2_tmp^2/(m4_tmp-m2_tmp^2);
            end
            naka_avg(:,:,i_repeat+(i_frame-1)*N_frame) = naka_map;
            perc = (i_x+(i_repeat-1+(i_frame-1)*N_repeat)*info.Nsc)/(info.Nsc*N_repeat*N_frame);
            msg = sprintf('CALCULATING NAKAGAMI PARAMETER %0.3f',perc);
            waitbar(perc,fwb,msg);
        end
    end
end
close(fwb)
naka_avg = median(naka_avg,3);
end
