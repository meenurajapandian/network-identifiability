function task_sensitivity_EC

load(fullfile(pwd,'data','optimal_recon.mat'));
load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));

%load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
%load('/home/mrajapan/Documents/MATLAB/data/Identifiability/optimal_recon.mat');

EC_opt = optimal.EC;
EC_FC_opt = optimal.EC_FC;

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

orig_ECT = cell(numTasks,1);
orig_ECRT = cell(numTasks,1);
recon_FC_ECT = cell(numTasks,1);
recon_FC_ECRT = cell(numTasks,1);
recon_ECT = cell(numTasks,1);
recon_ECRT = cell(numTasks,1);

for ii = 1:2:(T-2)
    % Reconstruct based on optimal reach
    orig_FCT = A{ii};
    orig_FCRT = A{ii+1};
    
    orig_matrix = zeros(configs.numEdges, configs.numFCs);
    orig_EC = zeros(configs.numEdges, configs.numFCs);

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

    D = zeros(n,n);
    mask_diag = logical(eye(n));

    for iii = 1:S
    
    aux = orig_FCT(:,:,iii);
    D(mask_diag) = sum(aux,2);
    orig_ECT{(ii+1)/2}(:,:,iii) = expm((D^-0.5) * aux * (D^ -0.5));
    aux = orig_FCRT(:,:,iii);
    D(mask_diag) = sum(aux,2);
    orig_ECRT{(ii+1)/2}(:,:,iii) = expm((D^-0.5) * aux * (D^ -0.5));

    aux = orig_ECT{(ii+1)/2}(:,:,iii);
    orig_EC(:,(2*iii)-1) = aux(mask_ut);
    aux = orig_ECRT{(ii+1)/2}(:,:,iii);
    orig_EC(:,2*iii) = aux(mask_ut);
    end

    clear orig_FCT orig_FCRT;
    
    orig_matrix_test = orig_matrix(:, Test_index);
    orig_matrix_retest = orig_matrix(:, Retest_index);

    orig_EC_test = orig_EC(:, Test_index);
    orig_EC_retest = orig_EC(:, Retest_index);

    % Decomposing and also reconstructing
    [FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
    [EC_modes, projected_EC_modes] = pca(orig_EC, 'NumComponents', configs.max_numPCs);

    ifc = EC_FC_opt((ii+1)/2,2);
    iec = EC_opt((ii+1)/2,2);
    % Reconstruct based on optimal for FC
    recon_matrix = projected_FC_modes(:,1:ifc) * FC_modes(:,1:ifc)';
    recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
    recon_matrix_test = recon_matrix(:, Test_index);
    recon_matrix_retest = recon_matrix(:, Retest_index);

    recon_EC = projected_EC_modes(:,1:iec) * EC_modes(:,1:iec)';
    recon_EC = bsxfun(@plus, recon_EC, mean(orig_EC));
    recon_EC_test = recon_EC(:, Test_index);
    recon_EC_retest = recon_EC(:, Retest_index);

    aux = zeros(n,n);
    recon_FCT = zeros(n,n);
    recon_FCRT = zeros(n,n);

    % Threshold
    recon_matrix_test(recon_matrix_test < eps) = eps;
    recon_matrix_retest(recon_matrix_retest < eps) = eps;

    D = zeros(n,n);
    mask_diag = logical(eye(n));

    for iii = 1:S
    aux = zeros(n,n);
    aux(mask_ut) = recon_matrix_test(:,iii);
    recon_FCT = aux + aux';
    aux(mask_ut) = recon_matrix_retest(:,iii);
    recon_FCRT = aux + aux';

    D(mask_diag) = sum(recon_FCT,2);
    recon_FC_ECT{(ii+1)/2}(:,:,iii) = expm((D^-0.5) * recon_FCT * (D^ -0.5));

    D(mask_diag) = sum(recon_FCRT,2);
    recon_FC_ECRT{(ii+1)/2}(:,:,iii) = expm((D^-0.5) * recon_FCRT * (D^ -0.5));
    
    aux = zeros(n,n);
    aux(mask_ut) = recon_EC_test(:,iii);
    recon_ECT{(ii+1)/2}(:,:,iii) = aux + aux';
    aux(mask_ut) = recon_EC_retest(:,iii);
    recon_ECRT{(ii+1)/2}(:,:,iii) = aux + aux';
    end
end


% Task sensitivity for both EC and EC on recon_FC
% recon_FC_ECT, recon_FC_ECRT, recon_ECT, recon_ECRT

recon_fc_ec_ts = zeros(n,n,S);
recon_ec_ts = zeros(n,n,S);
orig_ec_ts = zeros(n,n,S);
for iii = 1:S
    for i = 1:n
        for j = 1:n
            if i < j
                recon_fc_ec_M = zeros(numTasks,2);
                recon_ec_M = zeros(numTasks,2);
                ec_M = zeros(numTasks,2);
                for jj = 1:(numTasks)
                    recon_fc_ec_M(jj,1) = recon_FC_ECT{jj}(i,j,iii);
                    recon_fc_ec_M(jj,2) = recon_FC_ECRT{jj}(i,j,iii);
                    recon_ec_M(jj,1) = recon_ECT{jj}(i,j,iii);
                    recon_ec_M(jj,2) = recon_ECRT{jj}(i,j,iii);
                    ec_M(jj,1) = orig_ECT{jj}(i,j,iii);
                    ec_M(jj,2) = orig_ECRT{jj}(i,j,iii);
                end
                recon_fc_ec_ts(i,j,iii) = icc21(recon_fc_ec_M);
                recon_ec_ts(i,j,iii) = icc21(recon_ec_M);
                orig_ec_ts(i,j,iii) = icc21(ec_M);
            end
        end
    end
end

save('TS_EC.mat', 'recon_fc_ec_ts', 'recon_ec_ts','orig_ec_ts');

