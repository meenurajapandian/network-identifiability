load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
n = 374;
S = 100;
orig_FCT = A{1};
orig_FCRT = A{2};
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
orig_SP = zeros(configs.numEdges, configs.numFCs);
orig_MFPT = zeros(n*n, configs.numFCs);

% Building matrix to do PCA decomposition
for iii = 1:S
aux = orig_FCT(:,:,iii);
orig_matrix(:,(2*iii)-1) = aux(mask_ut);
aux = orig_FCRT(:,:,iii);
orig_matrix(:,2*iii) = aux(mask_ut);
end

orig_FCT(orig_FCT<eps) = eps;
orig_FCRT(orig_FCRT<eps) = eps;
orig_FCT(orig_FCT<eps) = eps;
orig_FCRT(orig_FCRT<eps) = eps;


for iii = 1:S
%[orig_SPT(:,:,iii), ~, ~] = get_shortest_path_lengths(1./orig_FCT(:,:,iii));
%[orig_SPRT(:,:,iii), ~, ~] = get_shortest_path_lengths(1./orig_FCRT(:,:,iii));
[~, orig_MFPTT(:,:,iii), ~] = f_mfpt(orig_FCT(:,:,iii));
[~, orig_MFPTRT(:,:,iii), ~] = f_mfpt(orig_FCRT(:,:,iii));

%aux = orig_SPT(:,:,iii);
%orig_SP(:,(2*iii)-1) = aux(mask_ut);
%aux = orig_SPRT(:,:,iii);
%orig_SP(:,2*iii) = aux(mask_ut);

aux = orig_MFPTT(:,:,iii);
orig_MFPT(:,(2*iii)-1) = aux(:);
aux = orig_MFPTRT(:,:,iii);
orig_MFPT(:,2*iii) = aux(:);
end

orig_matrix_test = orig_matrix(:, Test_index);
orig_matrix_retest = orig_matrix(:, Retest_index);

%orig_SP_test = orig_SP(:, Test_index);
%orig_SP_retest = orig_matrix(:, Retest_index);

orig_MFPT_test = orig_MFPT(:, Test_index);
orig_MFPT_retest = orig_MFPT(:, Retest_index);

% Decomposing and also reconstructing
[FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
%[SP_modes, projected_SP_modes] = pca(orig_SP, 'NumComponents', configs.max_numPCs);
[MFPT_modes, projected_MFPT_modes] = pca(orig_MFPT, 'NumComponents', configs.max_numPCs);

v = VideoWriter('yourFile.avi');
v.FrameRate = 25; % control frame-rate
open(v);
h=figure(1); %creates the figure and the handle to the figure

for i = 2:configs.numFCs
recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
recon_matrix_test = recon_matrix(:, Test_index);
recon_matrix_retest = recon_matrix(:, Retest_index);

%recon_SP = projected_SP_modes(:,1:i) * SP_modes(:,1:i)';
%recon_SP = bsxfun(@plus, recon_SP, mean(orig_SP));
%recon_SP_test = recon_SP(:, Test_index);
%recon_SP_retest = recon_SP(:, Retest_index);

recon_MFPT = projected_MFPT_modes(:,1:i) * MFPT_modes(:,1:i)';
recon_MFPT = bsxfun(@plus, recon_MFPT, mean(orig_MFPT));
recon_MFPT_test = recon_MFPT(:, Test_index);
recon_MFPT_retest = recon_MFPT(:, Retest_index);

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

% [recon_FC_SPT(:,:,iii), ~, ~] = get_shortest_path_lengths(1./recon_FCT(:,:,iii));
% [recon_FC_SPRT(:,:,iii), ~, ~] = get_shortest_path_lengths(1./recon_FCRT(:,:,iii));
[~, recon_FC_MFPTT(:,:,iii), ~] = f_mfpt(recon_FCT(:,:,iii));
[~, recon_FC_MFPTRT(:,:,iii), ~] = f_mfpt(recon_FCRT(:,:,iii));

% aux = recon_FC_SPT(:,:,iii);
% recon_FC_SP_test(:,iii) = aux(mask_ut);
% aux = recon_FC_SPRT(:,:,iii);
% recon_FC_SP_retest(:,iii) = aux(mask_ut);
aux = zeros(n,n);
aux = recon_FC_MFPTT(:,:,iii);
recon_FC_MFPT_test(:,iii) = aux(:);
aux = recon_FC_MFPTRT(:,:,iii);
recon_FC_MFPT_retest(:,iii) = aux(:);
end

mask_diag = logical(eye(S));
%Ident_MFPT_orig = corr(orig_MFPT_retest,orig_MFPT_test);
%Iself_MFPT_orig = mean(Ident_MFPT_orig(mask_diag));
%Iothers_MFPT_orig = mean(Ident_MFPT_orig(~mask_diag));
%Idscore_MFPT_orig = (Iself_MFPT_orig - Iothers_MFPT_orig) * 100;

Ident_MFPT_recon = corr(recon_MFPT_retest,recon_MFPT_test);
Iself_MFPT_recon = mean(Ident_MFPT_recon(mask_diag));
Iothers_MFPT_recon = mean(Ident_MFPT_recon(~mask_diag));
Idscore_MFPT_recon = (Iself_MFPT_recon - Iothers_MFPT_recon) * 100;

Ident_MFPT_FC_recon = corr(recon_FC_MFPT_retest,recon_FC_MFPT_test);
Iself_MFPT_FC_recon = mean(Ident_MFPT_FC_recon(mask_diag));
Iothers_MFPT_FC_recon = mean(Ident_MFPT_FC_recon(~mask_diag));
Idscore_MFPT_FC_recon = (Iself_MFPT_FC_recon - Iothers_MFPT_FC_recon) * 100;

figure(1);
subplot(2,2,1); imagesc(Ident_MFPT_recon); axis square; colorbar;
title({['Identifiability of fully recon MFPT'];['Identifiability ',num2str(Idscore_MFPT_recon)];['Iself =',num2str(Iself_MFPT_recon),' Iothers =', num2str(Iothers_MFPT_recon)]});
subplot(2,2,2); imagesc(Ident_MFPT_FC_recon); axis square; colorbar;
title({['Identifiability of MFPT on fully recon FC'];['Identifiability ',num2str(Idscore_MFPT_FC_recon)];['Iself =',num2str(Iself_MFPT_FC_recon),' Iothers =', num2str(Iothers_MFPT_FC_recon)]});

subplot(2,1,2); hold on;
plot(i,Idscore_MFPT_recon, 'o--r');
plot(i, Idscore_MFPT_FC_recon, 'o--b');xlim([1 200]); ylim([2 90]);


frame = getframe(h);
writeVideo(v,frame);
pause(2);
% finish loop
end
close(v)