function whats_happening_W(ii)
%load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));
n = 374;
S = 100;
orig_FCT = A{ii};
orig_FCRT = A{ii+1};
clear A;

mask_ut = triu(true(n,n),1);
mask_d_3 = logical(repmat(eye(n),1,1,S));

configs.numRegions = n; % aparc2009 parcellation. Number of brain regions
configs.numFCs = 2*S; % number of connectomes (2 FCs per subject in this case)
configs.numEdges = nnz(mask_ut);
configs.numVisits = 2; % 2 visits per subject (test-retest)
configs.max_numPCs = configs.numFCs; % maximum number of PCs == data dimension
Test_index = 1:configs.numVisits:configs.numFCs; % change this and next line if you have a different FCs ordering
Retest_index = 2:configs.numVisits:configs.numFCs;

orig_matrix = zeros(configs.numEdges, configs.numFCs);
orig_W = zeros(n*n, configs.numFCs);

% Building matrix to do PCA decomposition
for iii = 1:S
aux = orig_FCT(:,:,iii);
orig_matrix(:,(2*iii)-1) = aux(mask_ut);
aux = orig_FCRT(:,:,iii);
orig_matrix(:,2*iii) = aux(mask_ut);
end

orig_FCT(orig_FCT<eps) = eps;
orig_FCRT(orig_FCRT<eps) = eps;
orig_FCT(mask_d_3) = 0;
orig_FCRT(mask_d_3) = 0;

orig_WT = zeros(n,n,S);
orig_WRT = zeros(n,n,S);

for iii = 1:S
[SP, ~, ~] = get_shortest_path_lengths(1./orig_FCT(:,:,iii));
[~, MFPT, ~] = f_mfpt(orig_FCT(:,:,iii));
orig_WT(:,:,iii) = MFPT ./ SP;

[SP, ~, ~] = get_shortest_path_lengths(1./orig_FCRT(:,:,iii));
[~, MFPT, ~] = f_mfpt(orig_FCRT(:,:,iii));
orig_WRT(:,:,iii) = MFPT ./ SP;

aux = orig_WT(:,:,iii);
orig_W(:,(2*iii)-1) = aux(:);
aux = orig_WRT(:,:,iii);
orig_W(:,2*iii) = aux(:);
end

orig_matrix_test = orig_matrix(:, Test_index);
orig_matrix_retest = orig_matrix(:, Retest_index);

orig_W_test = orig_W(:, Test_index);
orig_W_retest = orig_W(:, Retest_index);

% Decomposing and also reconstructing
[FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
[W_modes, projected_W_modes] = pca(orig_W, 'NumComponents', configs.max_numPCs);

% v = VideoWriter('emotion_W.avi');
% v.FrameRate = 25; % control frame-rate
% open(v);
% h=figure(2); %creates the figure and the handle to the figure

Idscore_W_recon = zeros(configs.numFCs);
Idscore_W_FC_recon = zeros(configs.numFCs);

for i = 2:configs.numFCs
recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
recon_matrix_test = recon_matrix(:, Test_index);
recon_matrix_retest = recon_matrix(:, Retest_index);

recon_W = projected_W_modes(:,1:i) * W_modes(:,1:i)';
recon_W = bsxfun(@plus, recon_W, mean(orig_W));
recon_W_test = recon_W(:, Test_index);
recon_W_retest = recon_W(:, Retest_index);

aux = zeros(n,n);
recon_FCT = zeros(n,n,S);
recon_FCRT = zeros(n,n,S);

% Threshold
recon_matrix_test(recon_matrix_test < eps) = eps;
recon_matrix_retest(recon_matrix_retest < eps) = eps;

for iii = 1:S
aux = zeros(n,n);
aux(mask_ut) = recon_matrix_test(:,iii);
recon_FCT(:,:,iii) = aux + aux';
aux(mask_ut) = recon_matrix_retest(:,iii);
recon_FCRT(:,:,iii) = aux + aux';

[SP, ~, ~] = get_shortest_path_lengths(1./recon_FCT(:,:,iii));
[~, MFPT, ~] = f_mfpt(recon_FCT(:,:,iii));
recon_FC_WT(:,:,iii) = MFPT ./ SP;

[SP, ~, ~] = get_shortest_path_lengths(1./recon_FCRT(:,:,iii));
[~,MFPT,~] = f_mfpt(recon_FCRT(:,:,iii));
recon_FC_WRT(:,:,iii) = MFPT ./ SP;

aux = zeros(n,n);
aux = recon_FC_WT(:,:,iii);
recon_FC_W_test(:,iii) = aux(:);
aux = recon_FC_WRT(:,:,iii);
recon_FC_W_retest(:,iii) = aux(:);

end

mask_diag = logical(eye(S));

Ident_W_recon = corr(recon_W_retest,recon_W_test);
Iself_W_recon = mean(Ident_W_recon(mask_diag));
Iothers_W_recon = mean(Ident_W_recon(~mask_diag));
Idscore_W_recon(i) = (Iself_W_recon - Iothers_W_recon) * 100;

Ident_W_FC_recon = corr(recon_FC_W_retest,recon_FC_W_test);
Iself_W_FC_recon = mean(Ident_W_FC_recon(mask_diag));
Iothers_W_FC_recon = mean(Ident_W_FC_recon(~mask_diag));
Idscore_W_FC_recon(i) = (Iself_W_FC_recon - Iothers_W_FC_recon) * 100;

% figure(2);
% subplot(2,2,1); imagesc(Ident_W_recon); axis square; colorbar;
% title({['Identifiability of fully recon W'];['Identifiability ',num2str(Idscore_W_recon)];['Iself =',num2str(Iself_W_recon),' Iothers =', num2str(Iothers_W_recon)]});
% subplot(2,2,2); imagesc(Ident_W_FC_recon); axis square; colorbar;
% title({['Identifiability of W on fully recon FC'];['Identifiability ',num2str(Idscore_W_FC_recon)];['Iself =',num2str(Iself_W_FC_recon),' Iothers =', num2str(Iothers_W_FC_recon)]});
% 
% subplot(2,3,5); hold on;
% plot(i,Idscore_W_recon, 'o--r');
% plot(i, Idscore_W_FC_recon, 'o--b');xlim([1 200]); ylim([2 50]);


% frame = getframe(h);
% writeVideo(v,frame);
% finish loop
end
% close(v)
save(sprintf('Ident_W_%d.mat',ii), 'Idscore_W_recon', 'Idscore_W_FC_recon');
