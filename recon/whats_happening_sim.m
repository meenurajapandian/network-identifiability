function whats_happening_sim(ii)
%load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
tic;
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
orig_sosim = zeros(configs.numEdges, configs.numFCs);
orig_tasim = zeros(configs.numEdges, configs.numFCs);

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

orig_SST = zeros(n,n,S);
orig_SSRT = zeros(n,n,S);
orig_TST = zeros(n,n,S);
orig_TSRT = zeros(n,n,S);

for iii = 1:S
[SP, ~, ~] = get_shortest_path_lengths(1./orig_FCT(:,:,iii));
[~, MFPT, ~] = f_mfpt(orig_FCT(:,:,iii));
W = MFPT ./ SP;

sosim = W * W';
snorm = vecnorm(W');
sosim = sosim ./ snorm;
sosim = (sosim' ./ snorm)';
orig_SST(:,:,iii) = sosim;

tasim = W' * W;
tnorm = vecnorm(W);
tasim = tasim ./ tnorm;
tasim = (tasim' ./ tnorm)';
orig_TST(:,:,iii) = tasim;

[SP, ~, ~] = get_shortest_path_lengths(1./orig_FCRT(:,:,iii));
[~, MFPT, ~] = f_mfpt(orig_FCRT(:,:,iii));
W = MFPT ./ SP;

sosim = W * W';
snorm = vecnorm(W');
sosim = sosim ./ snorm;
sosim = (sosim' ./ snorm)';
orig_SSRT(:,:,iii) = sosim;

tasim = W' * W;
tnorm = vecnorm(W);
tasim = tasim ./ tnorm;
tasim = (tasim' ./ tnorm)';
orig_TSRT(:,:,iii) = tasim;

aux = orig_SST(:,:,iii);
orig_sosim(:,(2*iii)-1) = aux(mask_ut);
aux = orig_SSRT(:,:,iii);
orig_sosim(:,2*iii) = aux(mask_ut);

aux = orig_TST(:,:,iii);
orig_tasim(:,(2*iii)-1) = aux(mask_ut);
aux = orig_TSRT(:,:,iii);
orig_tasim(:,2*iii) = aux(mask_ut);
end

orig_matrix_test = orig_matrix(:, Test_index);
orig_matrix_retest = orig_matrix(:, Retest_index);

orig_sosim_test = orig_sosim(:, Test_index);
orig_sosim_retest = orig_sosim(:, Retest_index);
orig_tasim_test = orig_tasim(:, Test_index);
orig_tasim_retest = orig_tasim(:, Retest_index);


% Decomposing and also reconstructing
[FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
[SS_modes, projected_SS_modes] = pca(orig_sosim, 'NumComponents', configs.max_numPCs);
[TS_modes, projected_TS_modes] = pca(orig_tasim, 'NumComponents', configs.max_numPCs);

% v = VideoWriter('emotion_W.avi');
% v.FrameRate = 25; % control frame-rate
% open(v);
% h=figure(2); %creates the figure and the handle to the figure

Idscore_SS_recon = zeros(configs.numFCs);
Idscore_SS_FC_recon = zeros(configs.numFCs);
Idscore_TS_recon = zeros(configs.numFCs);
Idscore_TS_FC_recon = zeros(configs.numFCs);
%
for i = 2:configs.numFCs
recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
recon_matrix_test = recon_matrix(:, Test_index);
recon_matrix_retest = recon_matrix(:, Retest_index);

recon_sosim = projected_SS_modes(:,1:i) * SS_modes(:,1:i)';
recon_sosim = bsxfun(@plus, recon_sosim, mean(orig_sosim));
recon_sosim_test = recon_sosim(:, Test_index);
recon_sosim_retest = recon_sosim(:, Retest_index);

recon_tasim = projected_TS_modes(:,1:i) * TS_modes(:,1:i)';
recon_tasim = bsxfun(@plus, recon_tasim, mean(orig_tasim));
recon_tasim_test = recon_tasim(:, Test_index);
recon_tasim_retest = recon_tasim(:, Retest_index);

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
W = MFPT ./ SP;

sosim = W * W';
snorm = vecnorm(W');
sosim = sosim ./ snorm;
sosim = (sosim' ./ snorm)';
recon_FC_SST(:,:,iii) = sosim;

tasim = W' * W;
tnorm = vecnorm(W);
tasim = tasim ./ tnorm;
tasim = (tasim' ./ tnorm)';
recon_FC_TST(:,:,iii) = tasim;

[SP, ~, ~] = get_shortest_path_lengths(1./recon_FCRT(:,:,iii));
[~,MFPT,~] = f_mfpt(recon_FCRT(:,:,iii));
W = MFPT ./ SP;

sosim = W * W';
snorm = vecnorm(W');
sosim = sosim ./ snorm;
sosim = (sosim' ./ snorm)';
recon_FC_SSRT(:,:,iii) = sosim;

tasim = W' * W;
tnorm = vecnorm(W);
tasim = tasim ./ tnorm;
tasim = (tasim' ./ tnorm)';
recon_FC_TSRT(:,:,iii) = tasim;

aux = zeros(n,n);
aux = recon_FC_SST(:,:,iii);
recon_FC_sosim_test(:,iii) = aux(mask_ut);
aux = recon_FC_SSRT(:,:,iii);
recon_FC_sosim_retest(:,iii) = aux(mask_ut);

aux = recon_FC_TST(:,:,iii);
recon_FC_tasim_test(:,iii) = aux(mask_ut);
aux = recon_FC_TSRT(:,:,iii);
recon_FC_tasim_retest(:,iii) = aux(mask_ut);

end

mask_diag = logical(eye(S));

Ident_SS_recon = corr(recon_sosim_retest,recon_sosim_test);
Iself_SS_recon = mean(Ident_SS_recon(mask_diag));
Iothers_SS_recon = mean(Ident_SS_recon(~mask_diag));
Idscore_SS_recon(i) = (Iself_SS_recon - Iothers_SS_recon) * 100;

Ident_SS_FC_recon = corr(recon_FC_sosim_retest,recon_FC_sosim_test);
Iself_SS_FC_recon = mean(Ident_SS_FC_recon(mask_diag));
Iothers_SS_FC_recon = mean(Ident_SS_FC_recon(~mask_diag));
Idscore_SS_FC_recon(i) = (Iself_SS_FC_recon - Iothers_SS_FC_recon) * 100;

Ident_TS_recon = corr(recon_tasim_retest,recon_tasim_test);
Iself_TS_recon = mean(Ident_TS_recon(mask_diag));
Iothers_TS_recon = mean(Ident_TS_recon(~mask_diag));
Idscore_TS_recon(i) = (Iself_TS_recon - Iothers_TS_recon) * 100;

Ident_TS_FC_recon = corr(recon_FC_tasim_retest,recon_FC_tasim_test);
Iself_TS_FC_recon = mean(Ident_TS_FC_recon(mask_diag));
Iothers_TS_FC_recon = mean(Ident_TS_FC_recon(~mask_diag));Idscore_TS_FC_recon(i) = (Iself_TS_FC_recon - Iothers_TS_FC_recon) * 100;

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
save(sprintf('Ident_SS_%d.mat',ii), 'Idscore_SS_recon', 'Idscore_SS_FC_recon');
save(sprintf('Ident_TS_%d.mat',ii), 'Idscore_TS_recon', 'Idscore_TS_FC_recon');
toc;
