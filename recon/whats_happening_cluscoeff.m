function whats_happening_cluscoeff(ii)
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
orig_C = zeros(n, configs.numFCs);

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

orig_C_test = zeros(n,S);
orig_C_retest = zeros(n,S);

for iii = 1:S
orig_C_test(:,iii) = clustering_coef_wu(orig_FCT(:,:,iii));
orig_C_retest(:,iii) = clustering_coef_wu(orig_FCRT(:,:,iii));

orig_C(:,(2*iii)-1) = orig_C_test(:,iii);
orig_C(:,2*iii) = orig_C_retest(:,iii);
end

% Decomposing and also reconstructing
[FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
[C_modes, projected_C_modes] = pca(orig_C, 'NumComponents', configs.max_numPCs);

% v = VideoWriter('emotion_sp.avi');
% v.FrameRate = 25; % control frame-rate
% open(v);
% h=figure(2); %creates the figure and the handle to the figure

Idscore_C_recon = zeros(configs.numFCs,1);
Idscore_C_FC_recon = zeros(configs.numFCs,1);

for i = 2:configs.numFCs
recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
recon_matrix_test = recon_matrix(:, Test_index);
recon_matrix_retest = recon_matrix(:, Retest_index);

recon_C = projected_C_modes(:,1:i) * C_modes(:,1:i)';
recon_C = bsxfun(@plus, recon_C, mean(orig_C));
recon_C_test = recon_C(:, Test_index);
recon_C_retest = recon_C(:, Retest_index);

aux = zeros(n,n);
recon_FCT = zeros(n,n);
recon_FCRT = zeros(n,n);

% Threshold
recon_matrix_test(recon_matrix_test < eps) = eps;
recon_matrix_retest(recon_matrix_retest < eps) = eps;

for iii = 1:S
aux = zeros(n,n);
aux(mask_ut) = recon_matrix_test(:,iii);
recon_FCT = aux + aux';
aux(mask_ut) = recon_matrix_retest(:,iii);
recon_FCRT = aux + aux';

recon_FC_C_test(:,iii) = clustering_coef_wu(recon_FCT);
recon_FC_C_retest(:,iii) = clustering_coef_wu(recon_FCRT);

end

mask_diag = logical(eye(S));

Ident_C_recon = corr(recon_C_test,recon_C_retest);
Iself_C_recon = mean(Ident_C_recon(mask_diag));
Iothers_C_recon = mean(Ident_C_recon(~mask_diag));
Idscore_C_recon(i) = (Iself_C_recon - Iothers_C_recon) * 100;

Ident_C_FC_recon = corr(recon_FC_C_test,recon_FC_C_retest);
Iself_C_FC_recon = mean(Ident_C_FC_recon(mask_diag));
Iothers_C_FC_recon = mean(Ident_C_FC_recon(~mask_diag));
Idscore_C_FC_recon(i) = (Iself_C_FC_recon - Iothers_C_FC_recon) * 100;

% figure(2);
% subplot(2,2,1); imagesc(Ident_SP_recon); axis square; colorbar;
% title({['Identifiability of fully recon SP'];['Identifiability ',num2str(Idscore_SP_recon)];['Iself =',num2str(Iself_SP_recon),' Iothers =', num2str(Iothers_SP_recon)]});
% subplot(2,2,2); imagesc(Ident_SP_FC_recon); axis square; colorbar;
% title({['Identifiability of SP on fully recon FC'];['Identifiability ',num2str(Idscore_SP_FC_recon)];['Iself =',num2str(Iself_SP_FC_recon),' Iothers =', num2str(Iothers_SP_FC_recon)]});
% 
% subplot(2,3,5); hold on;
% plot(i,Idscore_SP_recon, 'o--r');
% plot(i, Idscore_SP_FC_recon, 'o--b');xlim([1 200]); ylim([2 50]);


% frame = getframe(h);
% writeVideo(v,frame);
% finish loop
end
% close(v)
save(sprintf('Ident_C_%d.mat',ii), 'Idscore_C_recon', 'Idscore_C_FC_recon');