clear
clc
%% LOAD FEATURE DATA
rootpath = '\\10.229.121.108\Workspace\YIXIAO\Eghbal Yun QPAT simulation\';
excelFileName = 'PAT_imaging_record.xlsx';
sheetName = 'ROI STATS V4 (2)';
[data, headers, raw] = xlsread(strcat(rootpath,excelFileName), sheetName);
features = headers(1:2,:);
D = data(:,[1,4:end]);
clearvars headers data
Nroi = size(D,1);
patientID = D(:,1);
ovarySide = raw(3:end,2);
GT = D(:,end);
Ncancer = sum(GT==0); Ncyst = sum(GT==1); Nfibroid = sum(GT==2); Nmetriosis = sum(GT==3); Nteratoma = sum(GT==4); Nnormal = sum(GT==5);
IdxCancer = find(GT==0); IdxCyst = find(GT==1); IdxFibroid = find(GT==2); IdxMetriosis = find(GT==3); IdxTeratoma = find(GT==4); IdxNormal = find(GT==5); 
fprintf('>>>>> Total ROI: %d, %d normal, %d cyst, %d fibroid, %d endometriosis, %d teratoma, %d cancer.\n', Nroi, Nnormal, Ncyst, Nfibroid, Nmetriosis, Nteratoma, Ncancer);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% CLASSIFICATION BY OVARY %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SORT DATA BY OVARY
Nfeatures = 116;
Npatients = unique(patientID,'stable');
Novary = 0;
Nov_normal = 0; Nov_cyst = 0; Nov_fibroid = 0; Nov_teratoma = 0; Nov_endo = 0; Nov_cancer = 0;
ovary_repeat_hist = [];
ovary_gt_hist = [];
for i = 1:length(Npatients)
    patientID_tmp = Npatients(i);
    roi_idx = find(patientID == patientID_tmp);
    side_all = ovarySide(roi_idx,1);
    side_unique = unique(side_all,'stable');
    Novary = Novary + (size(side_unique,1));    
    gt_all = GT(roi_idx);
    n_roi_idx = zeros(size(side_unique,1),1);
    for j = 1:(size(side_unique,1))
        side_tmp = side_unique(j); side_tmp = side_tmp{1,1};
        n_roi_tmp = sum(char(side_all{:,1}) == side_tmp);
        ovary_repeat_hist = [ovary_repeat_hist , n_roi_tmp];
        n_roi_idx(j,1) = n_roi_tmp;
    end
    n_roi_idx = cumsum(n_roi_idx);
    n_roi_idx = [1;n_roi_idx+1];
    for k = 1:(size(side_unique,1))
        gt_tmp = gt_all(n_roi_idx(k));
        ovary_gt_hist = [ovary_gt_hist , gt_tmp];
        if gt_tmp == 0; Nov_cancer = Nov_cancer+1;
        elseif gt_tmp == 1; Nov_cyst = Nov_cyst+1;
        elseif gt_tmp == 2; Nov_fibroid = Nov_fibroid+1;
        elseif gt_tmp == 3; Nov_endo = Nov_endo+1;
        elseif gt_tmp == 4; Nov_teratoma = Nov_teratoma+1;
        else; Nov_normal=Nov_normal+1;
        end
    end
end
fprintf('>>>>> Results include %d patients, %d ovaries.\n', length(Npatients), Novary);
fprintf('>>>>> %d normal ovaries, %d cyst, %d fibroid, %d endometriosis, %d teratoma, %d cancer.\n', Nov_normal, Nov_cyst, Nov_fibroid, Nov_endo, Nov_teratoma, Nov_cancer);

Davg = zeros(Novary,Nfeatures);
ovary_repeat_hist = [ovary_repeat_hist ; cumsum(ovary_repeat_hist)];
ovary_repeat_hist = [[0;0] ovary_repeat_hist];
for i = 1:Novary
    data_tmp = D((ovary_repeat_hist(2,i)+1):ovary_repeat_hist(2,i+1),2:117);
    Davg(i,:) = median(data_tmp,1);
end
IdxCancer = find(ovary_gt_hist==0); IdxCyst = find(ovary_gt_hist==1); IdxFibroid = find(ovary_gt_hist==2); IdxMetriosis = find(ovary_gt_hist==3); IdxTeratoma = find(ovary_gt_hist==4); IdxNormal = find(ovary_gt_hist==5); 

%% NORMALIZE FEATURE DATA
D1 = Davg([IdxCancer,IdxCyst,IdxFibroid,IdxTeratoma,IdxMetriosis,IdxNormal],:);
N1 = size(D1,1); % NUMBER OF ROIS
Nfeatures = size(D1,2);
featureName = features(2,4:end);
GT1 = ovary_gt_hist([IdxCancer,IdxCyst,IdxFibroid,IdxTeratoma,IdxMetriosis,IdxNormal]);
Mu1 = repmat(mean(D1,1),[N1,1]); Std1 = repmat(std(D1,1),[N1,1]);
D1norm = (D1 - Mu1)./Std1; % NORMALIZE VARIABLES

D1cancer = Davg(IdxCancer,:);
D1cyst = Davg(IdxCyst,:);
D1fibroid = Davg(IdxFibroid,:);
D1normal = Davg(IdxNormal,:);
D1endometriosis = Davg(IdxMetriosis,:);
D1teratoma = Davg(IdxTeratoma,:);

%% FIRST PASS: FEATURE SELECTION WITH ANOVA
features_pass1 = zeros(Nfeatures,1);
for idx_feature = 1:Nfeatures
    [p,tbl,stats] = anova1(D1norm(:,idx_feature),GT1,'off');
    [c,m] = multcompare(stats,'Display','on');
    fprintf('>>> Feature %d (%s) PAIRWISE P-VAL: CANCER: [%.5f,%.5f,%.5f] , CYST: [%.5f,%.5f,%.5f] , FIBROID: [%.5f,%.5f,%.5f]\n',...
        idx_feature, featureName{1,idx_feature}, c(1,6),c(2,6),c(3,6),c(1,6),c(4,6),c(5,6),c(2,6),c(4,6),c(6,6))
    x_in = 'a';
    f_box = figure;
    data_tmp = D1(:,idx_feature);
    lb_tmp = min(data_tmp); ub_tmp = max(data_tmp);
    subplot(411)
    h_cancer = histogram(D1cancer(:,idx_feature),'Normalization','pdf'); h_cancer.FaceColor = 'red'; h_cancer.FaceAlpha = 0.5; h_cancer.EdgeColor = 'none'; h_cancer.NumBins = 21; h_cancer.BinLimits = [lb_tmp,ub_tmp];
    subplot(412)
    h_cyst = histogram(D1cyst(:,idx_feature),'Normalization','pdf'); h_cyst.FaceColor = 'blue'; h_cyst.FaceAlpha = 0.5; h_cyst.EdgeColor = 'none'; h_cyst.NumBins = 21; h_cyst.BinLimits = [lb_tmp,ub_tmp];
    subplot(413)
    h_fibroid = histogram(D1fibroid(:,idx_feature),'Normalization','pdf'); h_fibroid.FaceColor = 'yellow'; h_fibroid.FaceAlpha = 0.5; h_fibroid.EdgeColor = 'none'; h_fibroid.NumBins = 21; h_fibroid.BinLimits = [lb_tmp,ub_tmp];
    subplot(414)
    h_normal = histogram(D1normal(:,idx_feature),'Normalization','pdf'); h_normal.FaceColor = 'green'; h_normal.FaceAlpha = 0.5; h_normal.EdgeColor = 'none'; h_normal.NumBins = 21; h_normal.BinLimits = [lb_tmp,ub_tmp];

    while strcmp(x_in,'y') == false && strcmp(x_in,'n') == false && strcmp(x_in,'b') == false
        x_in = input('DO YOU KEEP THIS FEATURE [Y\N]: \n','s');
        if strcmp(x_in,'y') == true
            features_pass1(idx_feature) = 1;
        elseif strcmp(x_in,'b') == true
            break
        end
    end
    close(f_box)
end
Nfeautures_pass1 = sum(features_pass1);
Nfeautures_pass1_shape = sum(features_pass1(1:13));
Nfeautures_pass1_mpus = sum(features_pass1(14:25));
Nfeautures_pass1_mppa = sum(features_pass1(26:33));
Nfeautures_pass1_usrad = sum(features_pass1(34:74));
Nfeautures_pass1_parad = sum(features_pass1(75:115));
fprintf('>>> %d / %d FEATURES REMAINS AFTER 1ST PASS: %d / 13 SHAPE, %d / 12 MPUS, %d / 8 MPPA, %d / 41 USRAD, %d / 41 PARAD\n',...
    sum(features_pass1) , Nfeatures ,  sum(features_pass1(1:13)) , sum(features_pass1(14:25)) , sum(features_pass1(26:33)) , sum(features_pass1(34:74)) , sum(features_pass1(75:115)))

features_pass1 = ones(116,1);
features_pass1_idx = find(features_pass1 == 1);
features_pass1_name = featureName(features_pass1_idx);

%% PLOT ALL 6 CATEGORIES
for idx_plot = 1:length(features_pass1_idx)
feature_name_tmp = features_pass1_name(idx_plot);
feature_name_tmp = feature_name_tmp{1,1};
L_spacing = 1.25;
x0 = L_spacing*ones(1,size(D1cancer,1));    x3 = 4*L_spacing*ones(1,size(D1teratoma,1));    x1 = 2*L_spacing*ones(1,size(D1cyst,1));    x2 = 3*L_spacing*ones(1,size(D1fibroid,1));     x4 = 5*L_spacing*ones(1,size(D1endometriosis,1));       x5 = 6*L_spacing*ones(1,size(D1normal,1));
y0 = D1cancer(:,features_pass1_idx(idx_plot))';     y3 = D1teratoma(:,features_pass1_idx(idx_plot))';   y1 = D1cyst(:,features_pass1_idx(idx_plot))'; y2 = D1fibroid(:,features_pass1_idx(idx_plot))';  y4 = D1endometriosis(:,features_pass1_idx(idx_plot))';  y5 = D1normal(:,features_pass1_idx(idx_plot))';

boxwidth = 0.85;
figure;
swarmchart(x0,y0,'MarkerEdgeColor',"#FF0000",'LineWidth',1,'MarkerFaceColor',"#FF0000",'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.5)
hold on
boxchart(x0,y0,'BoxEdgeColor','k','BoxFaceColor','none','MarkerColor','none','BoxWidth',boxwidth*L_spacing,'LineWidth',2)

swarmchart(x1,y1,'MarkerEdgeColor',"#0072BD",'LineWidth',1,'MarkerFaceColor',"#0072BD",'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.5)
boxchart(x1,y1,'BoxEdgeColor','k','BoxFaceColor','none','MarkerColor','none','BoxWidth',boxwidth*L_spacing,'LineWidth',2)

swarmchart(x2,y2,'MarkerEdgeColor',"#EDB120",'LineWidth',1,'MarkerFaceColor',"#EDB120",'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.5)
boxchart(x2,y2,'BoxEdgeColor','k','BoxFaceColor','none','MarkerColor','none','BoxWidth',boxwidth*L_spacing,'LineWidth',2)

swarmchart(x3,y3,'MarkerEdgeColor',"#7E2F8E",'LineWidth',1,'MarkerFaceColor',"#7E2F8E",'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.5)
boxchart(x3,y3,'BoxEdgeColor','k','BoxFaceColor','none','MarkerColor','none','BoxWidth',boxwidth*L_spacing,'LineWidth',2)

swarmchart(x4,y4,'MarkerEdgeColor',	"#FF00FF",'LineWidth',1,'MarkerFaceColor',"#FF00FF",'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.5)
boxchart(x4,y4,'BoxEdgeColor','k','BoxFaceColor','none','MarkerColor','none','BoxWidth',boxwidth*L_spacing,'LineWidth',2)

swarmchart(x5,y5,'MarkerEdgeColor',"#77AC30",'LineWidth',1,'MarkerFaceColor',"#77AC30",'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.5)
boxchart(x5,y5,'BoxEdgeColor','k','BoxFaceColor','none','MarkerColor','none','BoxWidth',boxwidth*L_spacing,'LineWidth',2)

set(gca,'xtick',L_spacing*[1,2,3,4,5,6],'xticklabel',{'Cancer','Cyst','Solid','Teratoma','Endometriosis','Normal'})
set(gca,'fontweight','bold','fontsize',14)
title(feature_name_tmp)
xlim([0.25*L_spacing,6.75*L_spacing])
%ylim([0.935,1.005])
hold off
end

%% SECOND PASS: FEATURE SELECTION WITH MRMR
feature_idx =  [26:33,75:115];
D1TRAIN_SUB = D1TRAIN(:,feature_idx);
featureName_sub=featureName(1,feature_idx);
[idx,scores] = fscmrmr(D1TRAIN_SUB,GTTRAIN);
N_feature_plot = 20;
bar(scores(idx(1:N_feature_plot)))
set(gca, 'XTick', 1:N_feature_plot); % center x-axis ticks on bins
set(gca, 'XTickLabel', featureName_sub(idx(1:N_feature_plot))); % set x-axis labels

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASS. MODEL TRAINING %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('off','all')
N_repeat = 50;
N_concat = 20;
N_class = 6;
wt_endometriosis = 5;
feature_indices = [9,24,52,30,32,96,102,116];
CM_test_hist_all = zeros(N_class,N_class,N_repeat);
x_roc_coord = linspace(0,1,101);
ROC_test_hist = zeros(N_repeat,101); auc_hist = zeros(N_repeat,1);
calibration_curve = [];
N_validation = [5,8,3,2,2,3];
model_selection = 'KNN';
fwb = waitbar(0,'Training KNN MODEL ... ');
for i_repeat = 1:N_repeat
    %% TRAIN-TEST SPLIT
    CM_test_hist = zeros(N_class,N_class,N_concat);
    p_test_hist = [];
    gt_test_hist = [];
    for j_concat = 1:N_concat
        test_idx_normal = sort(randperm(length(IdxNormal),N_validation(6)));      train_idx_normal = IdxNormal;   train_idx_normal(test_idx_normal) = [];     test_idx_normal = IdxNormal(test_idx_normal);
        test_idx_cyst = sort(randperm(length(IdxCyst),N_validation(2)));          train_idx_cyst = IdxCyst;       train_idx_cyst(test_idx_cyst) = [];         test_idx_cyst = IdxCyst(test_idx_cyst);
        test_idx_fibroid = sort(randperm(length(IdxFibroid),N_validation(3)));    train_idx_fibroid = IdxFibroid; train_idx_fibroid(test_idx_fibroid) = [];   test_idx_fibroid = IdxFibroid(test_idx_fibroid);
        test_idx_cancer = sort(randperm(length(IdxCancer),N_validation(1)));      train_idx_cancer = IdxCancer;   train_idx_cancer(test_idx_cancer) = [];     test_idx_cancer = IdxCancer(test_idx_cancer);
        test_idx_teratoma = sort(randperm(length(IdxTeratoma),N_validation(4)));  train_idx_teratoma = IdxTeratoma;   train_idx_teratoma(test_idx_teratoma) = [];     test_idx_teratoma = IdxTeratoma(test_idx_teratoma);
        test_idx_metriosis = sort(randperm(length(IdxMetriosis),N_validation(5)));  train_idx_metriosis = IdxMetriosis;   train_idx_metriosis(test_idx_metriosis) = [];     test_idx_metriosis = IdxMetriosis(test_idx_metriosis);

        % GENERATE TRAINING DATASET
        Mu1 = repmat(mean(Davg,1),[Novary,1]); Std1 = repmat(std(Davg,1),[Novary,1]);
        D1norm = (Davg - Mu1)./Std1; % NORMALIZE VARIABLES
        % BALANCE CATEGORIES IN THE TRAINING DATASET
        D1normal = D1norm(train_idx_normal,:); idx_tmp = randperm(length(train_idx_normal),1); D1normal = [D1normal;D1normal(idx_tmp,:)]; D1normal = repmat(D1normal,[3,1]);
        D1cyst = D1norm(train_idx_cyst,:); idx_tmp = randperm(length(train_idx_cyst),6); D1cyst = [D1cyst; D1cyst(idx_tmp,:)];
        D1fibroid = D1norm(train_idx_fibroid,:); idx_tmp = randperm(length(train_idx_fibroid),1); D1fibroid = [D1fibroid;D1fibroid(idx_tmp,:)];D1fibroid = repmat(D1fibroid,[3,1]);
        D1cancer = D1norm(train_idx_cancer,:); D1cancer = repmat(D1cancer,[2,1]); 
        D1teratoma = D1norm(train_idx_teratoma,:); D1teratoma = repmat(D1teratoma,[4,1]); idx_tmp = randi([1,length(train_idx_teratoma)*4],1,2); D1teratoma(idx_tmp,:)=[];
        D1metriosis = D1norm(train_idx_metriosis,:); D1metriosis = repmat(D1metriosis,[wt_endometriosis,1]);
        %idx_tmp = randperm(length(train_idx_cancer)*2,2); D1cancer(idx_tmp,:) = []; clearvars idx_tmp
        fprintf('>>>>>>>> Training set: normal %d, cyst %d, fibroid %d, cancer %d\n', size(D1normal,1), size(D1cyst,1), size(D1fibroid,1), size(D1cancer,1))
        fprintf('>>>>>>>> MAKE SURE TRAINING SET IS BALANCED!\n')
        D1TRAIN =   [D1cancer; D1cyst; D1fibroid; D1teratoma; D1metriosis; D1normal];
        GTTRAIN = [zeros(size(D1cancer,1),1); ones(size(D1cyst,1),1); 2*ones(size(D1fibroid,1),1); 3*ones(size(D1teratoma,1),1); 4*ones(size(D1metriosis,1),1); 5*ones(size(D1normal,1),1)];
        % GENERATE TESTING DATASET
        V1normal = D1norm(test_idx_normal,:); V1normal = repmat(V1normal,[5,1]);
        V1cyst = D1norm(test_idx_cyst,:); V1cyst = repmat(V1cyst,[2,1]); idx_tmp = randperm(length(test_idx_cyst)*2,1); V1cyst(idx_tmp,:) = [];
        V1fibroid = D1norm(test_idx_fibroid,:); V1fibroid = repmat(V1fibroid,[5,1]);
        V1cancer = D1norm(test_idx_cancer,:); V1cancer = repmat(V1cancer,[3,1]);
        V1teratoma = D1norm(test_idx_teratoma,:); V1teratoma = repmat(V1teratoma,[7,1]);
        V1metriosis = D1norm(test_idx_metriosis,:); V1metriosis = repmat(V1metriosis,[7,1]);
        D1TEST =   [V1cancer; V1cyst; V1fibroid; V1teratoma; V1metriosis; V1normal];
        GTTEST = [zeros(size(V1cancer,1),1); ones(size(V1cyst,1),1); 2*ones(size(V1fibroid,1),1); 3*ones(size(V1teratoma,1),1); 4*ones(size(V1metriosis,1),1); 5*ones(size(V1normal,1),1)];
        %Nsample_1class = size(V1normal,1);
        %% KNN MODEL
        D1TRAIN_SUB = D1TRAIN(:,feature_indices);
        if strcmp(model_selection, 'KNN') == 1
            testmodel = fitcknn(D1TRAIN_SUB, GTTRAIN,'NumNeighbors',20);
            D1TEST_SUB = D1TEST(:,feature_indices);
            [pred_test,score] = predict(testmodel,D1TEST_SUB);
        elseif strcmp(model_selection, 'SVM') == 1
            testmodel = fitcecoc(D1TRAIN_SUB, GTTRAIN);
            D1TEST_SUB = D1TEST(:,feature_indices);
            [pred_test,score] = predict(testmodel,D1TEST_SUB);
            score_norm = score - repmat(min(score,[],2),[1,N_class]);
            score_norm = score_norm./repmat(sum(score_norm,2),[1,N_class]);
        elseif strcmp(model_selection, 'MNR') == 1
            testmodel = fitmnr(D1TRAIN_SUB, GTTRAIN);     model_coeffs = table2array(testmodel.Coefficients(:,1))/1e2;
            D1TEST_SUB = D1TEST(:,feature_indices);
            p_cancer_tmp = model_coeffs(1:(length(feature_indices)+1),:)'*[ones(size(D1TEST_SUB,1),1) D1TEST_SUB]';
            p_cyst_tmp = model_coeffs((length(feature_indices)+2):2*(length(feature_indices)+1),:)'*[ones(size(D1TEST_SUB,1),1) D1TEST_SUB]';
            p_fibroid_tmp = model_coeffs((2*(length(feature_indices)+1)+1):3*(length(feature_indices)+1),:)'*[ones(size(D1TEST_SUB,1),1) D1TEST_SUB]';
            p_teratoma_tmp = model_coeffs((3*(length(feature_indices)+1)+1):4*(length(feature_indices)+1),:)'*[ones(size(D1TEST_SUB,1),1) D1TEST_SUB]';
            p_metriosis_tmp = model_coeffs((4*(length(feature_indices)+1)+1):5*(length(feature_indices)+1),:)'*[ones(size(D1TEST_SUB,1),1) D1TEST_SUB]';
            norm_factor = exp(p_cancer_tmp) + exp(p_cyst_tmp) + exp(p_fibroid_tmp) + exp(p_teratoma_tmp) + exp(p_metriosis_tmp) + 1;
            p_cancer = exp(p_cancer_tmp)./norm_factor;
            p_cyst = exp(p_cyst_tmp)./norm_factor;
            p_fibroid = exp(p_fibroid_tmp)./norm_factor;
            p_teratoma = exp(p_teratoma_tmp)./norm_factor;
            p_metriosis = exp(p_metriosis_tmp)./norm_factor;
            p_normal = 1./norm_factor;
            score  = [p_cancer' p_cyst' p_fibroid' p_teratoma' p_metriosis' p_normal'];
            [~,pred_test] = max(score,[],2);
            pred_test = pred_test-1;
        elseif strcmp(model_selection, 'MLP') == 1
            testmodel = patternnet([16,24,8]); testmodel.trainParam.showWindow=0;
            GTTRAIN_MLP = zeros(N_class, length(GTTRAIN));
            GTTRAIN_MLP(sub2ind(size(GTTRAIN_MLP), (GTTRAIN+1)', (1:length(GTTRAIN))))=1;
            testmodel = train(testmodel,D1TRAIN_SUB',GTTRAIN_MLP);
            D1TEST_SUB = D1TEST(:,feature_indices);
            score = testmodel(D1TEST_SUB'); score = score';
            [~,pred_test] = max(score,[],2);
            pred_test = pred_test-1;
        end
        CM_test = confusionmat(GTTEST, pred_test);
        CM_test = CM_test./repmat(sum(CM_test,2),[1,N_class]);
        CM_test_hist(:,:,j_concat) = CM_test;
        p_test_hist = [p_test_hist; score];
        gt_test_hist = [gt_test_hist; GTTEST];
    end
    CM_test_hist_all(:,:,i_repeat) = mean(CM_test_hist,3);
    score1vall = [p_test_hist(:,1), max(p_test_hist(:,2:end),[],2)];
    score1vall = p_test_hist(:,1)./sum(p_test_hist,2);
    calibration_curve = [calibration_curve; [score1vall , gt_test_hist==0]];
    %diffscore0 = p_test(:,1) - sum(p_test(:,2:end),2);
    %[X_roc,Y_roc] = perfcurve(GTTEST,diffscore0,0);
    rocObj = rocmetrics(gt_test_hist,p_test_hist,[0,1:N_class-1]);
    roc_Data = table2array(rocObj.Metrics);
    X_roc = roc_Data(1:round(size(roc_Data,1)/4),3);
    Y_roc = roc_Data(1:round(size(roc_Data,1)/4),4);
    [x_roc_unique,ia] = unique(X_roc); y_roc_unique = Y_roc(ia);
    y_roc = interp1(x_roc_unique, y_roc_unique, x_roc_coord);
    ROC_test_hist(i_repeat, :) = y_roc;
    auc_hist(i_repeat) = (sum(y_roc)*2 - y_roc(1) - y_roc(end))*(x_roc_coord(2) - x_roc_coord(1))/2;
    perc = i_repeat/N_repeat;
    msg = sprintf('TRAINING KNN MODEL: %d / %d', i_repeat, N_repeat);
    waitbar(perc,fwb,msg);
    clearvars p_cancer_tmp p_cancer p_cyst_tmp p_cyst p_fibroid_tmp p_fibroid norm_factor
end
close(fwb)
clearvars fwb

% PLOT AVERAGE CONFUSION MATRIX
CM_avg = mean(CM_test_hist_all,3);
wt_acc = (CM_avg(1,1)*(N_class-2) + sum(diag(CM_avg)))/10;
wt_ppv = (diag(CM_avg)'./sum(CM_avg,1)); wt_ppv = (wt_ppv(1)*4 + sum(wt_ppv))/10;
wt_std = reshape(CM_test_hist_all,[],N_repeat); wt_std = reshape(std(wt_std'),N_class, N_class);
wt_std = (sum(diag(wt_std))+wt_std(1,1)*(N_class-2));
figure
cmplot = confusionchart(round(CM_avg*100),{'0 cancer','1 cyst','2 solid','3 teratoma','4 endometriosis','5 normal'});
cmplot.FontSize = 12;
cmplot.RowSummary = 'row-normalized';
cmplot.ColumnSummary = 'column-normalized';
cmplot.DiagonalColor = "#0072BD";
cmplot.OffDiagonalColor = "#0072BD";
% PLOT AVERAGE ROC CURVE
ROC_test_hist(isnan(ROC_test_hist)) = 1;
auc_hist = (sum(ROC_test_hist,2)*2 - ROC_test_hist(:,1) - ROC_test_hist(:,end))*(x_roc_coord(2) - x_roc_coord(1))/2;
roc_avg = mean(medfilt2(ROC_test_hist,[1,1]),1); roc_std = std(medfilt2(ROC_test_hist,[1,1]),1);
auc_avg = mean(auc_hist); auc_std = std(auc_hist);

x2 = [x_roc_coord, fliplr(x_roc_coord)];
inBetween = [roc_avg-roc_std, fliplr(min(roc_avg+roc_std,1))];

figure
plot(x_roc_coord, roc_avg,'r','linewidth',2)
hold on
fill(x2, inBetween, [0.5,0.5,0.5],'FaceAlpha',0.5,'EdgeColor','none');
hold off
title(sprintf('Cancer vs rest: AUC = %.3f ± %.3f',auc_avg, auc_std))
set(gca,'fontweight','bold','fontsize',12)
xlim([-0.02,1.02])
ylim([-0.02,1.02])

roc_tosave = [roc_avg', roc_std'];
fprintf('>>>>>>>> CANCER AUC = %.3f, WEIGHTED ACC = %.3f, WEIGHTED PPV = %.3f, STD = %.4f\n', auc_avg, wt_acc, wt_ppv, wt_std);

%% DECISION CURVE ANALYSIS
thresholds_arr = linspace(0,1,26);
net_benefit = dca(calibration_curve(:,1), calibration_curve(:,2), thresholds_arr);
figure
plot(thresholds_arr , net_benefit)
net_benefit_all = [net_benefit_all, net_benefit'];

%% PLOT CALIBRATION CURVE
N_bin = 10;
bin_edges  = linspace(0,1,N_bin+1);
r_curve = zeros(N_bin,1);
ct_empirical = zeros(N_bin,1);
proba = discretize(calibration_curve(:,1),bin_edges);
for i_tmp =  1:N_bin
    idx_tmp =  find(proba == i_tmp);
    gt_tmp = calibration_curve(idx_tmp,2);
    ct_empirical(i_tmp) = length(idx_tmp);
    r_curve(i_tmp) = sum(gt_tmp)/length(idx_tmp);
end
uncert = 2e-2*log(sum(ct_empirical)./ct_empirical);
figure
plot(bin_edges(2:end)-0.05,r_curve,'x')
hold on
errorbar(bin_edges(2:end)-0.05,bin_edges(2:end)-0.05,uncert,'.')
plot(0:0.01:1,0:0.01:1)
hold off