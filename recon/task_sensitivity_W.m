function task_sensitivity_W

load(fullfile(pwd,'data','optimal_recon.mat'));
load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));

%load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
%load('/home/mrajapan/Documents/MATLAB/data/Identifiability/optimal_recon.mat');

W_opt = optimal.W;
W_FC_opt = optimal.W_FC;

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

orig_WT = cell(numTasks,1);
orig_WRT = cell(numTasks,1);
recon_FC_WT = cell(numTasks,1);
recon_FC_WRT = cell(numTasks,1);
recon_WT = cell(numTasks,1);
recon_WRT = cell(numTasks,1);

for ii = 1:2:(T-2)
    % Reconstruct based on optimal reach
    orig_FCT = A{ii};
    orig_FCRT = A{ii+1};
    
    orig_matrix = zeros(configs.numEdges, configs.numFCs);
    orig_W = zeros(n*n, configs.numFCs);

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
    
    [SP, ~, ~] = get_shortest_path_lengths(1./orig_FCT(:,:,iii));
    [~, MFPT, ~] = f_mfpt(orig_FCT(:,:,iii));
    orig_WT{(ii+1)/2}(:,:,iii) = MFPT ./ SP;

    [SP, ~, ~] = get_shortest_path_lengths(1./orig_FCRT(:,:,iii));
    [~, MFPT, ~] = f_mfpt(orig_FCRT(:,:,iii));
    orig_WRT{(ii+1)/2}(:,:,iii) = MFPT ./ SP;

    aux = orig_WT{(ii+1)/2}(:,:,iii);
    orig_W(:,(2*iii)-1) = aux(:);
    aux = orig_WRT{(ii+1)/2}(:,:,iii);
    orig_W(:,2*iii) = aux(:);
    
    end

    clear orig_FCT orig_FCRT;
    
    orig_matrix_test = orig_matrix(:, Test_index);
    orig_matrix_retest = orig_matrix(:, Retest_index);

    orig_W_test = orig_W(:, Test_index);
    orig_W_retest = orig_W(:, Retest_index);

    % Decomposing and also reconstructing
    [FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
    [W_modes, projected_W_modes] = pca(orig_W, 'NumComponents', configs.max_numPCs);

    ifc = W_FC_opt((ii+1)/2,2);
    iw = W_opt((ii+1)/2,2);
    % Reconstruct based on optimal for FC
    recon_matrix = projected_FC_modes(:,1:ifc) * FC_modes(:,1:ifc)';
    recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
    recon_matrix_test = recon_matrix(:, Test_index);
    recon_matrix_retest = recon_matrix(:, Retest_index);

    recon_W = projected_W_modes(:,1:iw) * W_modes(:,1:iw)';
    recon_W = bsxfun(@plus, recon_W, mean(orig_W));
    recon_W_test = recon_W(:, Test_index);
    recon_W_retest = recon_W(:, Retest_index);

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

    [SP, ~, ~] = get_shortest_path_lengths(1./recon_FCT);
    [~, MFPT, ~] = f_mfpt(recon_FCT);
    recon_FC_WT{(ii+1)/2}(:,:,iii) = MFPT ./ SP;

    [SP, ~, ~] = get_shortest_path_lengths(1./recon_FCRT);
    [~,MFPT,~] = f_mfpt(recon_FCRT);
    recon_FC_WRT{(ii+1)/2}(:,:,iii) = MFPT ./ SP;
    
    aux = zeros(n,n);
    aux(:) = recon_W_test(:,iii);
    recon_WT{(ii+1)/2}(:,:,iii) = aux ;
    aux = zeros(n,n);
    aux(:) = recon_W_retest(:,iii);
    recon_WRT{(ii+1)/2}(:,:,iii) = aux;
    end
end


% Task sensitivity for both W and W on recon_FC
% recon_FC_WT, recon_FC_WRT, recon_WT, recon_WRT

recon_fc_w_ts = zeros(n,n,S);
recon_w_ts = zeros(n,n,S);
orig_w_ts = zeros(n,n,S);
for iii = 1:S
    for i = 1:n
        for j = 1:n
%            if i < j
                recon_fc_w_M = zeros(numTasks,2);
                recon_w_M = zeros(numTasks,2);
                w_M = zeros(numTasks,2);
                for jj = 1:(numTasks)
                    recon_fc_w_M(jj,1) = recon_FC_WT{jj}(i,j,iii);
                    recon_fc_w_M(jj,2) = recon_FC_WRT{jj}(i,j,iii);
                    recon_w_M(jj,1) = recon_WT{jj}(i,j,iii);
                    recon_w_M(jj,2) = recon_WRT{jj}(i,j,iii);
                    w_M(jj,1) = orig_WT{jj}(i,j,iii);
                    w_M(jj,2) = orig_WRT{jj}(i,j,iii);
                end
                recon_fc_w_ts(i,j,iii) = icc21(recon_fc_w_M);
                recon_w_ts(i,j,iii) = icc21(recon_w_M);
                orig_w_ts(i,j,iii) = icc21(w_M);
%            end
        end
    end
end

save('TS_W_1.mat', 'recon_fc_w_ts');
save('TS_W_2.mat', 'recon_w_ts');
save('TS_W_3.mat', 'orig_w_ts');


