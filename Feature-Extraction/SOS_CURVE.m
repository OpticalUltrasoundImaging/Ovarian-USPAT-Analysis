function [CAVG,CLOC] = SOS_CURVE(Reconstruction_hist, RxMux, info, lambdaG, lambdaR, N_iter)
fwb = waitbar(0,'CALCULATING AVERAGE SOS ...');
dr = info.c/info.fc/2;
z_start = 0.25e-2; 
if info.Nfocus==2816
    z_end = 5.25e-2; % [M]
elseif info.Nfocus == 2304
    z_end = 4.25e-2;
elseif info.Nfocus == 3328
    z_end = 6.25e-2;
elseif info.Nfocus < 2000
    z_end = 3.25e-2;
elseif info.Nfocus == 5376
    z_end = 9.5e-2;
elseif info.Nfocus == 4864
    z_end = 8.75e-2;
elseif info.Nfocus == 4352
    z_end = 8.25e-2;
else
    z_end = 7.25e-2;
end
zgrid_range = (z_start):dr:(z_end); % [M]
%daq_i = round(z_start/info.pixel_d); daq_f = round(z_end/info.pixel_d);
L_cut = 22;
f = [0 0.6 0.6 1]; m = [1 1 0 0];
LPF = fir2(64,f,m);
clearvars f m
N_frame = size(Reconstruction_hist,4);
CAVG = zeros(length(zgrid_range),info.Nsc,N_frame);
for i_frame = 1:N_frame
    Reconstruction = Reconstruction_hist(:,:,:,i_frame);
    for idx_r = 1:info.Nsc
        sc = idx_r;                         sc_ang = info.ScanAngle(sc)*pi/180;
        sc_data = Reconstruction(:,:,sc);
        sc_elem_lookup = RxMux(sc,:);
        start_idx = find(sc_elem_lookup>0, 1); end_idx = find(sc_elem_lookup<info.N_ele+1, 1, 'last');
        elem_start = sc_elem_lookup(start_idx);
        elem_end = sc_elem_lookup(end_idx);
        L_ch_active = elem_end-elem_start+1;
        sc_data = sc_data(:,start_idx:end_idx);
        sc_data = conv2(sc_data,LPF','same');
        % UPSAMPLE SIGNAL FROM 40MHZ TO 200MHZ
        sc_data_upsmp = zeros(size(sc_data,1)*5,size(sc_data,2));
        for i = 1:L_ch_active
            sc_tmp = sc_data(:,i);
            sc_data_upsmp(:,i) = interp1(0:info.Nfocus-1,sc_tmp,0:0.2:(info.Nfocus-0.2));
        end
        %figure; imagesc(sc_data)
        parfor idx_z = 1:length(zgrid_range)
            N_neighbor = 3;
            L_axialKernel = 381;
            kernel_decay = exp(-2e-2*abs(-floor(L_axialKernel/2):floor(L_axialKernel/2)))';
            dx1x2 = zeros((L_ch_active-N_neighbor)*N_neighbor,1);
            z_pos = zgrid_range(idx_z);
            z_pos_daq = min(round(z_pos/info.pixel_d*5),info.Nfocus*5);
            z_strip = sc_data_upsmp((z_pos_daq-floor(L_axialKernel/2)):(z_pos_daq+floor(L_axialKernel/2)),:);
            % APPLY EXPONENTIAL DECAY TO SIGNAL TO AVOID PEAK HOPPING WHEN COMPUTING CORRELATION
            for i = 1:L_ch_active
                z_strip(:,i) = kernel_decay.*z_strip(:,i);
            end
            %figure; imagesc(z_strip)
            M = zeros((L_ch_active-N_neighbor)*N_neighbor,L_ch_active);
            for x1 = 1:(L_ch_active-N_neighbor)
                sx1 = z_strip(:,x1);
                M(((x1-1)*N_neighbor+1):x1*N_neighbor,x1) = 1;
                for x2 = (x1+1):(x1+N_neighbor)
                    sx2 = z_strip(:,x2);
                    M(((x1-1)*N_neighbor+x2-x1),x2) = -1;
                    r_tau = xcorr(sx1,sx2,50); % MAX 50 PIXELS DELAY BETWEEN SIGNALS FROM 2 ELEMENTS
                    [~,max_idx] = max(r_tau);
                    dx1x2((x1-1)*N_neighbor+x2-x1,1) = max_idx-51;
                end
            end
            t = lsqnonneg(M,dx1x2/200); % [S]
            xe = zeros(L_ch_active,1); ze = xe; xesign = xe;
            for i = 1:L_ch_active
                elem_tmp = elem_start+i-1;
                ele_ang = info.EleAngle(elem_tmp)*pi/180;
                ch_ang = abs(sc_ang - ele_ang);
                xesign(i) = sign(ele_ang - sc_ang);
                xe(i) = info.ROC*sin(ch_ang);
                ze(i) = info.ROC*(1 - cos(ch_ang)) + z_pos;
            end
            elem_pos_tmp = xe.*xe + ze.*ze;
            %elem_pos_tmp = (sqrt(xe.*xe + ze.*ze) + z_pos).^2; % [M^2]
            %figure; yyaxis left; plot(1:L_ch_active,t); yyaxis right; plot(1:L_ch_active,sqrt(elem_pos_tmp))
            [~,min_idx] = min(elem_pos_tmp);
            elem_pos_tmp_fit = elem_pos_tmp(max(1,min_idx-L_cut):min(L_ch_active,min_idx+L_cut+1));
            t_fit = smoothdata(t(max(1,min_idx-L_cut):min(L_ch_active,min_idx+L_cut+1)),'sgolay',5);
            %figure; yyaxis left; plot(1:(2*L_cut+2),t_fit); yyaxis right; plot(1:1:(2*L_cut+2),sqrt(elem_pos_tmp_fit))
            f = fit(elem_pos_tmp_fit,t_fit.*t_fit,'poly1');
            %f = fit(elem_pos_tmp,t.^2,'poly2');
            %figure; plot(f,elem_pos_tmp,t.^2)
            if f.p1<=0
                cavg_tmp = NaN;
            else
                cavg_tmp = 1e3/sqrt(f.p1); % [MM/US]
            end
            CAVG(idx_z,idx_r,i_frame) = cavg_tmp;
        end
        CAVG(:,idx_r,i_frame) = fillmissing(CAVG(:,idx_r,i_frame),'movmedian',[5,5]);
        CAVG(:,idx_r,i_frame) = filloutliers(medfilt1(CAVG(:,idx_r,i_frame),5),"clip","movmedian",11);
        perc = (idx_r + (i_frame-1)*info.Nsc)/(info.Nsc*N_frame);
        msg = sprintf('CALCULATING AVERAGE SOS %d / %d SCANLINES',idx_r + (i_frame-1)*info.Nsc,info.Nsc*N_frame);
        waitbar(perc,fwb,msg);
    end
end
close(fwb)

disp('>>>>>>>> ESTIMATION OF AVERAGE SOS COMPLETE.')
%imagesc((CAVG),[1.3,1.7]); colorbar

test = medfilt2(median(CAVG,3));
test = fillmissing(fillmissing(test,'movmedian',[7,7]),'makima',2);
CAVG = test;
CAVG(isnan(CAVG)) = 0;
CAVG(isinf(CAVG)) = 0;
%imagesc(test/1e4,[0,2])
%%
%fwb = waitbar(0,'CALCULATING LOCAL SOS ...');
Ntap = 11; fcutoff = 0.1;
lof = sin(2*pi*fcutoff*(-floor(Ntap/2):floor(Ntap/2))).*hamming(Ntap); lof = lof/sum(lof); % 11-TAP LOWPASS FILTER BASED ON HAMMING WINDOW WITH 0.2 FREQUENCY CUTOFF
% lambdaG = 2e-3;
% lambdaR = 0.2e-3;
N = size(CAVG,1);
A = tril(ones(N));
norm_factor = 1./(1:N); norm_factor = norm_factor';
A = A.*repmat(norm_factor, [1 N]);
c_init = 0.75; % INITIAL GUESS OF SOS
% N_iter = 1e4;
CLOC = c_init*ones(size(CAVG));
% %loss_hist = zeros(N_iter,1);
% fprintf('Progress:\n');
% fprintf(['\n' repmat('.',1,info.Nsc) '\n\n']);
% mst_drop = 0.2;
% d1 = designfilt("lowpassfir", ...
%     PassbandFrequency=0.005,StopbandFrequency=0.20, ...
%     PassbandRipple=1,StopbandAttenuation=60, ...
%     DesignMethod="equiripple");
% parfor sc = 1:info.Nsc % SOLVE FOR LOCAL SOS BY SCANLINE
%     cavg_tmp = filtfilt(d1,CAVG(:,sc));
% 
%     c_init = mean(cavg_tmp);
%     cloc_tmp = c_init*ones(N,1);
% %     cavg_tmp = ones(462,1); cavg_tmp(225:462,:) = 1+0.005*(1:238); cloc_tmp = 1.5*ones(462,1);
% %     lambdaG = 1e-4;
% %     lambdaR = 2e-5;
% %     N_iter = 4e4;
%     lambdaG_iter = lambdaG; lambdaR_iter = lambdaR;
%     for idx_iter = 1:N_iter
%         %loss_hist(idx_iter) = (conv(cloc_tmp,lof,'same') - cloc_tmp)'*(conv(cloc_tmp,lof,'same') - cloc_tmp);
%         G_iter = 2*A'*(A*cloc_tmp - cavg_tmp); % MODEL FIDELITY LOSS
%         R_iter = (conv(cloc_tmp,lof,'same') - cloc_tmp); % REGULARIZATION LOSS
%         if idx_iter == 30001 || idx_iter == 37501 || idx_iter == 50001
%             lambdaG_iter = mst_drop*lambdaG_iter; lambdaR_iter = 2.5*mst_drop*lambdaR_iter;
%         end
%         cloc_tmp = cloc_tmp - lambdaG*G_iter + lambdaR*R_iter;
%         cloc_tmp(cloc_tmp<0) = 0;
%     end
% %     figure
% %     plot(1:length(zgrid_range),cavg_tmp,'linewidth',2)
% %     hold on
% %     plot(1:length(zgrid_range),cloc_tmp,'linewidth',2)
% %     hold off
%     %figure; plot(1:N_iter,loss_hist)
%     %figure; plot(1:length(zgrid_range),cavg_tmp); hold on; plot(1:length(zgrid_range),cloc_tmp); hold off; legend('c avg','c loc')
%     CLOC(:,sc) = cloc_tmp;
%     %perc = sc/info.Nsc;
%     %msg = sprintf('CALCULATING LOCAL SOS %d / %d SCANLINES',sc,info.Nsc);
%     %waitbar(perc,fwb,msg);
%     fprintf('\b|\n');
% end
% disp('>>>>>>>> CALCULATION OF LOCAL SOS COMPLETE.')
%close(fwb)
% figure
% subplot(121); imagesc(CAVG,[1.3,1.7]); colorbar
% subplot(122); imagesc(CLOC,[1.3,1.7]); colorbar

end
