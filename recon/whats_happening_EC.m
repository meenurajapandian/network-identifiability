function whats_happening_EC(ii)
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
orig_EC = zeros(configs.numEdges, configs.numFCs);

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

orig_ECT = zeros(n,n,S);
orig_ECRT = zeros(n,n,S);
D = zeros(n,n);
mask_diag = logical(eye(n));

for iii = 1:S
aux = orig_FCT(:,:,iii);
D(mask_diag) = sum(aux,2);
orig_ECT(:,:,iii) = expm((D^-0.5) * aux * (D^ -0.5));
aux = orig_FCRT(:,:,iii);
D(mask_diag) = sum(aux,2);
orig_ECRT(:,:,iii) = expm((D^-0.5) * aux * (D^ -0.5));

aux = orig_ECT(:,:,iii);
orig_EC(:,(2*iii)-1) = aux(mask_ut);
aux = orig_ECRT(:,:,iii);
orig_EC(:,2*iii) = aux(mask_ut);
end

orig_matrix_test = orig_matrix(:, Test_index);
orig_matrix_retest = orig_matrix(:, Retest_index);

orig_EC_test = orig_EC(:, Test_index);
orig_EC_retest = orig_EC(:, Retest_index);

% Decomposing and also reconstructing
[FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
[EC_modes, projected_EC_modes] = pca(orig_EC, 'NumComponents', configs.max_numPCs);

% v = VideoWriter('emotion_sp.avi');
% v.FrameRate = 25; % control frame-rate
% open(v);
% h=figure(2); %creates the figure and the handle to the figure

Idscore_EC_recon = zeros(configs.numFCs);
Idscore_EC_FC_recon = zeros(configs.numFCs);

for i = 2:configs.numFCs
recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
recon_matrix_test = recon_matrix(:, Test_index);
recon_matrix_retest = recon_matrix(:, Retest_index);

recon_EC = projected_EC_modes(:,1:i) * EC_modes(:,1:i)';
recon_EC = bsxfun(@plus, recon_EC, mean(orig_EC));
recon_EC_test = recon_EC(:, Test_index);
recon_EC_retest = recon_EC(:, Retest_index);

aux = zeros(n,n);
recon_FCT = zeros(n,n,S);
recon_FCRT = zeros(n,n,S);

% Threshold
recon_matrix_test(recon_matrix_test < eps) = eps;
recon_matrix_retest(recon_matrix_retest < eps) = eps;

D = zeros(n,n);
mask_diag = logical(eye(n));

for iii = 1:S
aux = zeros(n,n);
aux(mask_ut) = recon_matrix_test(:,iii);
recon_FCT(:,:,iii) = aux + aux';
aux(mask_ut) = recon_matrix_retest(:,iii);
recon_FCRT(:,:,iii) = aux + aux';


D(mask_diag) = sum(recon_FCT(:,:,iii),2);
recon_FC_ECT(:,:,iii) = expm((D^-0.5) * recon_FCT(:,:,iii) * (D^ -0.5));

D(mask_diag) = sum(recon_FCRT(:,:,iii),2);
recon_FC_ECRT(:,:,iii) = expm((D^-0.5) * recon_FCRT(:,:,iii) * (D^ -0.5));

aux = zeros(n,n);
aux = recon_FC_ECT(:,:,iii);
recon_FC_EC_test(:,iii) = aux(mask_ut);
aux = recon_FC_ECRT(:,:,iii);
recon_FC_EC_retest(:,iii) = aux(mask_ut);

end

mask_diag = logical(eye(S));

Ident_EC_recon = corr(recon_EC_retest,recon_EC_test);
Iself_EC_recon = mean(Ident_EC_recon(mask_diag));
Iothers_EC_recon = mean(Ident_EC_recon(~mask_diag));
Idscore_EC_recon(i) = (Iself_EC_recon - Iothers_EC_recon) * 100;

Ident_EC_FC_recon = corr(recon_FC_EC_retest,recon_FC_EC_test);
Iself_EC_FC_recon = mean(Ident_EC_FC_recon(mask_diag));
Iothers_EC_FC_recon = mean(Ident_EC_FC_recon(~mask_diag));
Idscore_EC_FC_recon(i) = (Iself_EC_FC_recon - Iothers_EC_FC_recon) * 100;

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
save(sprintf('Ident_EC_%d.mat',ii), 'Idscore_EC_recon', 'Idscore_EC_FC_recon');