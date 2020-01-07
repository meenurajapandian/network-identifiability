function task_sensitivity_SP

load(fullfile(pwd,'data','optimal_recon.mat'));
load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));

SP_opt = optimal.SP;
SP_FC_opt = optimal.SP_FC;

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

orig_SPT = cell(numTasks,1);
orig_SPRT = cell(numTasks,1);
recon_FC_SPT = cell(numTasks,1);
recon_FC_SPRT = cell(numTasks,1);
recon_SPT = cell(numTasks,1);
recon_SPRT = cell(numTasks,1);

for ii = 1:2:(T-2)
    % Reconstruct based on optimal reach
    orig_FCT = A{ii};
    orig_FCRT = A{ii+1};
    
    orig_matrix = zeros(configs.numEdges, configs.numFCs);
    orig_SP = zeros(configs.numEdges, configs.numFCs);

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
    [orig_SPT{(ii+1)/2}(:,:,iii), ~, ~] = get_shortest_path_lengths(1./orig_FCT(:,:,iii));
    [orig_SPRT{(ii+1)/2}(:,:,iii), ~, ~] = get_shortest_path_lengths(1./orig_FCRT(:,:,iii));

    aux = orig_SPT{(ii+1)/2}(:,:,iii);
    orig_SP(:,(2*iii)-1) = aux(mask_ut);
    aux = orig_SPRT{(ii+1)/2}(:,:,iii);
    orig_SP(:,2*iii) = aux(mask_ut);
    end

    clear orig_FCT orig_FCRT;
    
    orig_matrix_test = orig_matrix(:, Test_index);
    orig_matrix_retest = orig_matrix(:, Retest_index);

    orig_SP_test = orig_SP(:, Test_index);
    orig_SP_retest = orig_SP(:, Retest_index);

    % Decomposing and also reconstructing
    [FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
    [SP_modes, projected_SP_modes] = pca(orig_SP, 'NumComponents', configs.max_numPCs);

    ifc = SP_FC_opt((ii+1)/2,2);
    isp = SP_opt((ii+1)/2,2);
    % Reconstruct based on optimal for FC
    recon_matrix = projected_FC_modes(:,1:ifc) * FC_modes(:,1:ifc)';
    recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
    recon_matrix_test = recon_matrix(:, Test_index);
    recon_matrix_retest = recon_matrix(:, Retest_index);

    recon_SP = projected_SP_modes(:,1:isp) * SP_modes(:,1:isp)';
    recon_SP = bsxfun(@plus, recon_SP, mean(orig_SP));
    recon_SP_test = recon_SP(:, Test_index);
    recon_SP_retest = recon_SP(:, Retest_index);

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

    [recon_FC_SPT{(ii+1)/2}(:,:,iii), ~, ~] = get_shortest_path_lengths(1./recon_FCT);
    [recon_FC_SPRT{(ii+1)/2}(:,:,iii), ~, ~] = get_shortest_path_lengths(1./recon_FCRT);
    
    aux = zeros(n,n);
    aux(mask_ut) = recon_SP_test(:,iii);
    recon_SPT{(ii+1)/2}(:,:,iii) = aux + aux';
    aux = zeros(n,n);
    aux(mask_ut) = recon_SP_retest(:,iii);
    recon_SPRT{(ii+1)/2}(:,:,iii) = aux + aux';
    end
end


% Task sensitivity for both SP and SP on recon_FC
% recon_FC_SPT, recon_FC_SPRT, recon_SPT, recon_SPRT

recon_fc_sp_ts = zeros(n,n,S);
recon_sp_ts = zeros(n,n,S);
orig_sp_ts = zeros(n,n,S);
for iii = 1:S
    for i = 1:n
        for j = 1:n
            if i < j
                recon_fc_sp_M = zeros(numTasks,2);
                recon_sp_M = zeros(numTasks,2);
                sp_M = zeros(numTasks,2);
                for jj = 1:(numTasks)
                    recon_fc_sp_M(jj,1) = recon_FC_SPT{jj}(i,j,iii);
                    recon_fc_sp_M(jj,2) = recon_FC_SPRT{jj}(i,j,iii);
                    recon_sp_M(jj,1) = recon_SPT{jj}(i,j,iii);
                    recon_sp_M(jj,2) = recon_SPRT{jj}(i,j,iii);
                    sp_M(jj,1) = orig_SPT{jj}(i,j,iii);
                    sp_M(jj,2) = orig_SPRT{jj}(i,j,iii);
                end
                recon_fc_sp_ts(i,j,iii) = icc21(recon_fc_sp_M);
                recon_sp_ts(i,j,iii) = icc21(recon_sp_M);
                orig_sp_ts(i,j,iii) = icc21(sp_M);
            end
        end
    end
end

save('TS_SP.mat', 'recon_fc_sp_ts', 'recon_sp_ts','orig_sp_ts');

