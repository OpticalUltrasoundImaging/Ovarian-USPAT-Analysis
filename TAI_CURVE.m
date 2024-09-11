function [fc_map, tai_avg, fc_map_var] = TAI_CURVE(R_raw_hist,info)
l_pad = 2^9;
x = linspace(0,16,l_pad/2+1); x = x(2:end);
freq_template = chi2pdf(x,4);
[~,max_idx_template] = max(freq_template);
clearvars x
%plot(1:512,freq_template)
pulse_length = (1.5*info.c/info.fc);
N_repeat = 3;     N_frame = size(R_raw_hist,3);
w = 0.5e-3;       w_roi = round(w/info.pitch);        if mod(w_roi,2)==0; w_roi = w_roi+1; end

fc_test = 7.1e6;
fc_scaled_tmp = round(fc_test/info.fs*l_pad);
fc_map = fc_test*ones(size(R_raw_hist,1),size(R_raw_hist,2),N_repeat*N_frame);
fc_map_var = 50*ones(size(R_raw_hist,1),size(R_raw_hist,2),N_repeat*N_frame);
% f = [0 0.1 0.1 1]; m = [1 0.9 0.9 0];
% AC_cancel = fir2(64,f,m);
f = [0 0.6 0.6 1]; m = [1 1 0 0];
LPF = fir2(64,f,m);
f_list = 1e-6*linspace(0.5*fc_test,0.75*fc_test,round(0.25*fc_scaled_tmp)+1)';
Ntap = 11; fcutoff = 0.2;
lof = sin(2*pi*fcutoff*(-floor(Ntap/2):floor(Ntap/2))).*hamming(Ntap); lof = lof/sum(lof); % 11-TAP LOWPASS FILTER BASED ON HAMMING WINDOW WITH 0.2 FREQUENCY CUTOFF

A = [f_list.^2,  -f_list,  ones(size(f_list))]; 
fwb = waitbar(0,'CALCULATING CENTER FREQUENCY ... ');
for i_frame = 1:N_frame
    R_env_raw = R_raw_hist(:,:,i_frame);
    R_env_raw = conv2(conv2(R_env_raw,LPF','same'),lof','same');
    for i_repeat = 1:N_repeat
        i_factor = 9+i_repeat;
        h_roi = round(i_factor*pulse_length/info.pixel_d);      if mod(h_roi,2)==0; h_roi = h_roi+1; end
        for i_x = 1:info.Nsc
            i_x_start = max(1,i_x-floor(w_roi/2));
            i_x_end = min(i_x+floor(w_roi/2),info.Nsc);
            fc_aline = zeros(info.Nfocus,1); fcvar_aline = zeros(info.Nfocus,1);
            parfor i_z = 1:info.Nfocus
                i_z_start = max(1,i_z-floor(h_roi/2));
                i_z_end = min(i_z+floor(h_roi/2),info.Nfocus);
                roi_tmp = R_env_raw(i_z_start:i_z_end,i_x_start:i_x_end);
                [h_roi_tmp,w_roi_tmp] = size(roi_tmp);
                roi_tmp_pad = [roi_tmp;zeros(l_pad-h_roi_tmp,w_roi_tmp)];
                roi_tmp_fft = mean(abs(fftshift(fft(roi_tmp_pad,[],1),1)),2);

                roi_tmp_fft = flip(roi_tmp_fft(1:l_pad/2))+roi_tmp_fft(l_pad/2+1:l_pad);
                roi_tmp_fft = conv(smoothdata(roi_tmp_fft,'lowess',29),lof,'same');
                roi_tmp_fft = roi_tmp_fft/max(roi_tmp_fft(:));
                [r_test,lag_test] = xcorr(roi_tmp_fft,freq_template,l_pad/4);
                [~,max_idx] = max(r_test);
                lag_max = lag_test(max_idx)+max_idx_template;
                fc = lag_max*20/l_pad*2;

                index1 = find(roi_tmp_fft >= 0.4, 1, 'first');  index2 = find(roi_tmp_fft >= 0.4, 1, 'last');
                fwhm = (index2-index1 + 1)*20/512;
                %plot(lag_test,r_test)
                %fft_bw = roi_tmp_fft(floor(0.5*fc_scaled_tmp):ceil(0.75*fc_scaled_tmp));
                %fft_bw = 5 - fft_bw/(max(fft_bw)+1e-8);

%                 x = lsqnonneg(A,sqrt(fft_bw));
%                 fc = x(2)/x(1)/2;
%                 fc_var = 1/x(1);
                %                 figure
                %                 plot(f_list,b)
                %                 hold on
                %                 plot(f_list,x(1)*f_list.^2+x(2)*f_list+x(3))
                %                 hold off
                %fc = sum((linspace(0.425*fc_test,1.675*fc_test,round(1.25*fc_scaled_tmp)+1)'.*fft_bw))/sum(fft_bw);
                fc_aline(i_z) = fc;
                fcvar_aline(i_z) = fwhm;
            end
            fc_map(:,i_x,i_repeat+(i_frame-1)*N_repeat) = fc_aline;
            fc_map_var(:,i_x,i_repeat+(i_frame-1)*N_repeat) = fcvar_aline;
            perc = (i_x+(i_repeat-1+(i_frame-1)*N_repeat)*info.Nsc)/(info.Nsc*N_repeat*N_frame);
            msg = sprintf('CALCULATING CENTER FREQ %0.3f',perc);
            waitbar(perc,fwb,msg);
        end
    end
end
close(fwb)
fc_map_var = median(abs(fc_map_var),3);
fc_map = median(abs(fc_map),3);
%imagesc(fc_map,[0,10])
%figure; subplot(121); imagesc(fc_map); subplot(122); imagesc(fc_map_var)

%%
fc_map_smooth = medfilt2(fillmissing(fc_map,'makima'));
fc_map_smooth(fc_map_smooth<0)=0;
fc_map_smooth = conv2(fc_map_smooth,lof,'same');
%fc_map_smooth = medfilt2(fc_map./sqrt(fc_map_var));
%figure; imagesc(fc_map_smooth)

TAI = zeros(size(fc_map_smooth,1),size(fc_map_smooth,2),7);
fwb = waitbar(0,'CALCULATING ATTENUATION COEFF...');
for i_repeat = 1:5
    i_factor = 27+i_repeat;
    z_kernel = round(i_factor*pulse_length/info.pixel_d);      if mod(z_kernel,2)==0; z_kernel = z_kernel+1; end
    z_half = floor(z_kernel/2);
    z_coord_kernel = 0:info.pixel_d:(z_kernel-1)*info.pixel_d;
    for i_sc = 1:info.Nsc
        parfor i_z = 1:info.Nfocus
            zstrip_tmp = fc_map_smooth(max(1,i_z-z_half):min(info.Nfocus,i_z+z_half),i_sc);
            B = 1e2*[ones(1,length(zstrip_tmp));z_coord_kernel(1:length(zstrip_tmp))]'\zstrip_tmp;
            TAI(i_z,i_sc,i_repeat) = B(2);
        end
        perc = (i_sc+(i_repeat-1)*info.Nsc)/(info.Nsc*7);
        msg = sprintf('CALCULATING ATTENUTATION COEFF %0.3f',perc);
        waitbar(perc,fwb,msg);
    end
    TAI(:,:,i_repeat) = fillmissing(TAI(:,:,i_repeat),"movmean",[15,15]);
end
close(fwb)

%tai_std = std(TAI,0,3);
TAI(TAI>0) = NaN;
tai_avg = mean(TAI,3); tai_avg = -1*tai_avg;

tai_avg = medfilt2(fillmissing(fillmissing(tai_avg,'movmedian',11),'constant',0),[7,3]);
%tai_avg = medfilt2(fillmissing(fillmissing(tai_avg,"movmedian",[13,13]),'movmedian',[1,1],2),[7,3]);
tai_avg = tai_avg./4*2.217;
%tai_avg(tai_avg<0) = 0;
tai_avg = 1.5*tai_avg;
tai_avg(tai_avg<0.025) = 0.025;
%imagesc(tai_avg)

% tai_avg = tai_avg
% figure; subplot(121); imagesc(tai_avg); subplot(122); imagesc(tai_std)
% TAI = -1*TAI; TAI(TAI<-1) = -1; TAI = TAI+1;
% TAI = TAI*2.217/12.5;
%figure;imagesc(TAI)

% dz_filter = designfilt('differentiatorfir','FilterOrder',34, ...
%     'PassbandFrequency',0.25,'StopbandFrequency',2.5, ...
%     'SampleRate',20);
% dz_filter = dz_filter.Coefficients;
% plot(1:35,dz_filter)
% TAI = -1*(conv2(fc_map_smooth,dz_filter', 'same')/(info.pixel_d*1e2)/12.5);
% TAI(TAI<0) = 0;
% TAI = 2.217*imbilatfilt(ordfilt2(TAI,120,ones(125,1)),0.1*max(TAI(:))^2,2);
%figure; imagesc(TAI,[0,1])
end