function info = loadCSystemParam(Roi,System)
info.c = 1540;
info.fs = 40e6;

info.N_ele = double(System.Transducer.elementCnt);
info.pitch = System.Transducer.elementPitchCm/100;
info.fc = (System.Transducer.frequencyMHz)*1e6;
info.ROC = System.Transducer.radiusCm/100;
info.ele_width = System.Transducer.elementWidthCm/100;
info.ele_height = 6e-3;
info.pixel_d = info.c/info.fs;

info.Nsc = 128;
info.N_ch = 128;

Offset = 0;
info.Nfocus = 2048;
info.fc_scaled = info.fc/info.fs*info.Nfocus/2;

info.RxFnum = 1;
info.FOV = Roi{1}.lateralLength;
info.st_ang = Roi{1}.lateralStart;
info.d_theta = info.FOV/(info.Nsc-1);

info.one_angle = info.FOV/(info.N_ele-1);

info.ScanAngle = [info.st_ang:info.d_theta:info.st_ang+(info.Nsc-1)*info.d_theta];
tmp = -info.FOV/2;
info.EleAngle = [tmp:info.one_angle:-tmp];

info.half_rx_ch = info.N_ch*info.pitch*0.5; % total radius of each scan line

n_sample = [0:info.Nfocus-1]' + Offset;
info.d_sample = n_sample * info.pixel_d;
end
