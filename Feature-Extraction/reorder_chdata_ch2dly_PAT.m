function [reordered_data] = reorder_chdata_ch2dly_PAT(sc_idx, N_ele, N_sc, nNumCh, data_128ch)

data128ch = data_128ch;
MuxRatio = N_ele/N_sc;

sc_idx1 = floor(sc_idx*MuxRatio);
start_idx = sc_idx1-nNumCh/2;
k = rem(N_ele*3 + start_idx, N_ele) + 1;
idx1 = [k:1:128];
idx2 = [1:1:k-1];
reorder_idx= [idx1 idx2];
    
reordered_data(:,:) = data128ch(reorder_idx,:); 

% trim side data
trim_idx = [start_idx : 1 : start_idx+nNumCh-1];
reordered_data(find(trim_idx < 0),:) = 0;
reordered_data(find(trim_idx > N_ele-1),:) = 0;



