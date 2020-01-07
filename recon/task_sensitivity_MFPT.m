function task_sensitivity_MFPT

load(fullfile(pwd,'data','optimal_recon.mat'));
load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));

%load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
%load('/home/mrajapan/Documents/MATLAB/data/Identifiability/optimal_recon.mat');

MFPT_opt = optimal.MFPT;
MFPT_FC_opt = optimal.MFPT_FC;

clear optimal;

n = 374;
S = 100;
T = 18;
numTasks = 8;

mask_ut = triu(true(n,n),1);
mask_d_3 = logical(repmat(eye(n),1,1,S));

configs.numRegions = n; % aparc2009 parcellation. Number of brain regions
configs.numFCs = 2*S; % number of connectomes (2 FCs per subject in this case)
configs.numEdges = nnz(mask_ut);
configs.numVisits = 2; % 2 visits per subject (test-retest)
configs.max_numPCs = configs.numFCs; % maximum number of PCs == data dimension
Test_index = 1:configs.numVisits:configs.numFCs; % change this and next line if you have a different FCs ordering
Retest_index = 2:configs.numVisits:configs.numFCs;

orig_MFPTT = cell(numTasks,1);
orig_MFPTRT = cell(numTasks,1);
recon_FC_MFPTT = cell(numTasks,1);
recon_FC_MFPTRT = cell(numTasks,1);
recon_MFPTT = cell(numTasks,1);
recon_MFPTRT = cell(numTasks,1);

for ii = 1:2:(T-2)
    % Reconstruct based on optimal reach
    orig_FCT = A{ii};
    orig_FCRT = A{ii+1};
    
    orig_matrix = zeros(configs.numEdges, configs.numFCs);
    orig_MFPT = zeros(n*n, configs.numFCs);

    % Building FC matrix to do PCA decomposition
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

    for iii = 1:S
    [~, orig_MFPTT{(ii+1)/2}(:,:,iii), ~] = f_mfpt(orig_FCT(:,:,iii));
    [~, orig_MFPTRT{(ii+1)/2}(:,:,iii), ~] = f_mfpt(orig_FCRT(:,:,iii));

    aux = orig_MFPTT{(ii+1)/2}(:,:,iii);
    orig_MFPT(:,(2*iii)-1) = aux(:);
    aux = orig_MFPTRT{(ii+1)/2}(:,:,iii);
    orig_MFPT(:,2*iii) = aux(:);
    end

    clear orig_FCT orig_FCRT;
    
    orig_matrix_test = orig_matrix(:, Test_index);
    orig_matrix_retest = orig_matrix(:, Retest_index);

    orig_MFPT_test = orig_MFPT(:, Test_index);
    orig_MFPT_retest = orig_MFPT(:, Retest_index);

    % Decomposing and also reconstructing
    [FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
    [MFPT_modes, projected_MFPT_modes] = pca(orig_MFPT, 'NumComponents', configs.max_numPCs);

    ifc = MFPT_FC_opt((ii+1)/2,2);
    imfpt = MFPT_opt((ii+1)/2,2);
    % Reconstruct based on optimal for FC
    recon_matrix = projected_FC_modes(:,1:ifc) * FC_modes(:,1:ifc)';
    recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
    recon_matrix_test = recon_matrix(:, Test_index);
    recon_matrix_retest = recon_matrix(:, Retest_index);

    recon_MFPT = projected_MFPT_modes(:,1:imfpt) * MFPT_modes(:,1:imfpt)';
    recon_MFPT = bsxfun(@plus, recon_MFPT, mean(orig_MFPT));
    recon_MFPT_test = recon_MFPT(:, Test_index);
    recon_MFPT_retest = recon_MFPT(:, Retest_index);

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
    aux = zeros(n,n);
    aux(mask_ut) = recon_matrix_retest(:,iii);
    recon_FCRT = aux + aux';

    [~, recon_FC_MFPTT{(ii+1)/2}(:,:,iii), ~] = f_mfpt(recon_FCT);
    [~, recon_FC_MFPTRT{(ii+1)/2}(:,:,iii), ~] = f_mfpt(recon_FCRT);
    
    aux = zeros(n,n);
    aux(:) = recon_MFPT_test(:,iii);
    recon_MFPTT{(ii+1)/2}(:,:,iii) = aux;
    aux = zeros(n,n);
    aux(:) = recon_MFPT_retest(:,iii);
    recon_MFPTRT{(ii+1)/2}(:,:,iii) = aux;
    end
end


% Task sensitivity for both MFPT and MFPT on recon_FC
% recon_FC_MFPTT, recon_FC_MFPTRT, recon_MFPTT, recon_MFPTRT

recon_fc_mfpt_ts = zeros(n,n,S);
recon_mfpt_ts = zeros(n,n,S);
orig_mfpt_ts = zeros(n,n,S);
for iii = 1:S
    for i = 1:n
        for j = 1:n
%            if i < j
                recon_fc_mfpt_M = zeros(numTasks,2);
                recon_mfpt_M = zeros(numTasks,2);
                mfpt_M = zeros(numTasks,2);
                for jj = 1:(numTasks)
                    recon_fc_mfpt_M(jj,1) = recon_FC_MFPTT{jj}(i,j,iii);
                    recon_fc_mfpt_M(jj,2) = recon_FC_MFPTRT{jj}(i,j,iii);
                    recon_mfpt_M(jj,1) = recon_MFPTT{jj}(i,j,iii);
                    recon_mfpt_M(jj,2) = recon_MFPTRT{jj}(i,j,iii);
                    mfpt_M(jj,1) = orig_MFPTT{jj}(i,j,iii);
                    mfpt_M(jj,2) = orig_MFPTRT{jj}(i,j,iii);
                end
                recon_fc_mfpt_ts(i,j,iii) = icc21(recon_fc_mfpt_M);
                recon_mfpt_ts(i,j,iii) = icc21(recon_mfpt_M);
                orig_mfpt_ts(i,j,iii) = icc21(mfpt_M);
%            end
        end
    end
end

save('TS_MFPT.mat', 'recon_fc_mfpt_ts', 'recon_mfpt_ts','orig_mfpt_ts');

