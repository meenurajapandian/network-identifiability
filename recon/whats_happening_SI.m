function whats_happening_SI(ii)
%load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));
orig_FCT = A{ii};
orig_FCRT = A{ii+1};
clear A;
n = 374;
S = 100;

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
orig_SI = zeros(configs.numEdges, configs.numFCs);

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

orig_SIT = zeros(n,n,S);
orig_SIRT = zeros(n,n,S);

for iii = 1:S
[~, SP, SPB] = get_shortest_path_lengths(1./orig_FCT(:,:,iii));
orig_SIT(:,:,iii) = get_information_shortest_paths_wei_und(orig_FCT(:,:,iii),SP,SPB,sum((orig_FCT(:,:,iii))),1);

[~, SP, SPB] = get_shortest_path_lengths(1./orig_FCRT(:,:,iii));
orig_SIRT(:,:,iii) = get_information_shortest_paths_wei_und(orig_FCRT(:,:,iii),SP,SPB,sum((orig_FCRT(:,:,iii))),1);


aux = orig_SIT(:,:,iii);
orig_SI(:,(2*iii)-1) = aux(mask_ut);
aux = orig_SIRT(:,:,iii);
orig_SI(:,2*iii) = aux(mask_ut);
end

orig_matrix_test = orig_matrix(:, Test_index);
orig_matrix_retest = orig_matrix(:, Retest_index);

orig_SI_test = orig_SI(:, Test_index);
orig_SI_retest = orig_SI(:, Retest_index);

% Decomposing and also reconstructing
[FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
[SI_modes, projected_SI_modes] = pca(orig_SI, 'NumComponents', configs.max_numPCs);

v = VideoWriter('emotion_si.avi');
v.FrameRate = 25; % control frame-rate
open(v);
h=figure(2); %creates the figure and the handle to the figure

Idscore_SI_recon = zeros(configs.numFCs);
Idscore_SI_FC_recon = zeros(configs.numFCs);

for i = 2:configs.numFCs
recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
recon_matrix_test = recon_matrix(:, Test_index);
recon_matrix_retest = recon_matrix(:, Retest_index);

recon_SI = projected_SI_modes(:,1:i) * SI_modes(:,1:i)';
recon_SI = bsxfun(@plus, recon_SI, mean(orig_SI));
recon_SI_test = recon_SI(:, Test_index);
recon_SI_retest = recon_SI(:, Retest_index);

aux = zeros(n,n);
recon_FCT = zeros(n,n,S);
recon_FCRT = zeros(n,n,S);

% Threshold
recon_matrix_test(recon_matrix_test < eps) = 10^-10;
recon_matrix_retest(recon_matrix_retest < eps) = 10^-10;

for iii = 1:S
aux = zeros(n,n);
aux(mask_ut) = recon_matrix_test(:,iii);
recon_FCT(:,:,iii) = aux + aux';
aux(mask_ut) = recon_matrix_retest(:,iii);
recon_FCRT(:,:,iii) = aux + aux';

[~, SP, SPB] = get_shortest_path_lengths(1./recon_FCT(:,:,iii));
recon_FC_SIT(:,:,iii) = get_information_shortest_paths_wei_und(recon_FCT(:,:,iii),SP,SPB,sum((recon_FCT(:,:,iii))),1);

[~, SP, SPB] = get_shortest_path_lengths(1./recon_FCRT(:,:,iii));
recon_FC_SIRT(:,:,iii) = get_information_shortest_paths_wei_und(recon_FCRT(:,:,iii),SP,SPB,sum((recon_FCRT(:,:,iii))),1);

aux = zeros(n,n);
aux = recon_FC_SIT(:,:,iii);
recon_FC_SI_test(:,iii) = aux(mask_ut);
aux = recon_FC_SIRT(:,:,iii);
recon_FC_SI_retest(:,iii) = aux(mask_ut);

end

mask_diag = logical(eye(S));

Ident_SI_recon = corr(recon_SI_retest,recon_SI_test);
Iself_SI_recon = mean(Ident_SI_recon(mask_diag));
Iothers_SI_recon = mean(Ident_SI_recon(~mask_diag));
Idscore_SI_recon(i) = (Iself_SI_recon - Iothers_SI_recon) * 100;

Ident_SI_FC_recon = corr(recon_FC_SI_retest,recon_FC_SI_test);
Iself_SI_FC_recon = mean(Ident_SI_FC_recon(mask_diag));
Iothers_SI_FC_recon = mean(Ident_SI_FC_recon(~mask_diag));
Idscore_SI_FC_recon(i) = (Iself_SI_FC_recon - Iothers_SI_FC_recon) * 100;

% figure(2);
% subplot(2,2,1); imagesc(Ident_SI_recon); axis square; colorbar;
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
save(sprintf('Ident_SI_%d.mat',ii), 'Idscore_SI_recon', 'Idscore_SI_FC_recon');
