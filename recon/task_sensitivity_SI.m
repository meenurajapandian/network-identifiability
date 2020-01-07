function task_sensitivity_SI

load(fullfile(pwd,'data','optimal_recon.mat'));
load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));

%load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
%load('/home/mrajapan/Documents/MATLAB/data/Identifiability/optimal_recon.mat');

SI_opt = optimal.SI;
SI_FC_opt = optimal.SI_FC;

clear optimal;

n = 374;
S = 100;
T = 18;
numTasks = 8;

mask_ut = triu(true(n,n),1);
mask_d_3 = logical(repmat(eye(n),1,1,S));

configs.numRegions = n; % aparc2009 parcellation. Number of brain regions
configs.numFCs = 2*S; % number of connectomes (2 FCs per subjsit in this case)
configs.numEdges = nnz(mask_ut);
configs.numVisits = 2; % 2 visits per subjsit (test-retest)
configs.max_numPCs = configs.numFCs; % maximum number of PCs == data dimension
Test_index = 1:configs.numVisits:configs.numFCs; % change this and next line if you have a different FCs ordering
Retest_index = 2:configs.numVisits:configs.numFCs;

orig_SIT = cell(numTasks,1);
orig_SIRT = cell(numTasks,1);
recon_FC_SIT = cell(numTasks,1);
recon_FC_SIRT = cell(numTasks,1);
recon_SIT = cell(numTasks,1);
recon_SIRT = cell(numTasks,1);

for ii = 1:2:(T-2)
    % Reconstruct based on optimal reach
    orig_FCT = A{ii};
    orig_FCRT = A{ii+1};
    
    orig_matrix = zeros(configs.numEdges, configs.numFCs);
    orig_SI = zeros(n*n, configs.numFCs); %Reconstruct entire matrix because of asymmetry

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

    [~, SP, SPB] = get_shortest_path_lengths(1./orig_FCT(:,:,iii));
    orig_SIT{(ii+1)/2}(:,:,iii) = get_information_shortest_paths_wei_und(orig_FCT(:,:,iii),SP,SPB,sum((orig_FCT(:,:,iii))),1);

    [~, SP, SPB] = get_shortest_path_lengths(1./orig_FCRT(:,:,iii));
    orig_SIRT{(ii+1)/2}(:,:,iii) = get_information_shortest_paths_wei_und(orig_FCRT(:,:,iii),SP,SPB,sum((orig_FCRT(:,:,iii))),1);

    aux = orig_SIT{(ii+1)/2}(:,:,iii);
    orig_SI(:,(2*iii)-1) = aux(:);
    aux = orig_SIRT{(ii+1)/2}(:,:,iii);
    orig_SI(:,2*iii) = aux(:);
    end

    clear orig_FCT orig_FCRT;
    
    orig_matrix_test = orig_matrix(:, Test_index);
    orig_matrix_retest = orig_matrix(:, Retest_index);

    orig_SI_test = orig_SI(:, Test_index);
    orig_SI_retest = orig_SI(:, Retest_index);

    % Decomposing and also reconstructing
    [FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);
    [SI_modes, projected_SI_modes] = pca(orig_SI, 'NumComponents', configs.max_numPCs);

    ifc = SI_FC_opt((ii+1)/2,2);
    isi = SI_opt((ii+1)/2,2);
    % Reconstruct based on optimal for FC
    recon_matrix = projected_FC_modes(:,1:ifc) * FC_modes(:,1:ifc)';
    recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
    recon_matrix_test = recon_matrix(:, Test_index);
    recon_matrix_retest = recon_matrix(:, Retest_index);

    recon_SI = projected_SI_modes(:,1:isi) * SI_modes(:,1:isi)';
    recon_SI = bsxfun(@plus, recon_SI, mean(orig_SI));
    recon_SI_test = recon_SI(:, Test_index);
    recon_SI_retest = recon_SI(:, Retest_index);

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
   
    [~, SP, SPB] = get_shortest_path_lengths(1./recon_FCT);
    recon_FC_SIT{(ii+1)/2}(:,:,iii) = get_information_shortest_paths_wei_und(recon_FCT,SP,SPB,sum(recon_FCT),1);

    [~, SP, SPB] = get_shortest_path_lengths(1./recon_FCRT);
    recon_FC_SIRT{(ii+1)/2}(:,:,iii) = get_information_shortest_paths_wei_und(recon_FCRT,SP,SPB,sum(recon_FCRT),1);

    aux = zeros(n,n);
    aux(:) = recon_SI_test(:,iii);
    recon_SIT{(ii+1)/2}(:,:,iii) = aux;
    aux = zeros(n,n);
    aux(:) = recon_SI_retest(:,iii);
    recon_SIRT{(ii+1)/2}(:,:,iii) = aux;
    end
end


% Task sensitivity for both SI and SI on rsion_FC
% recon_FC_SIT, recon_FC_SIRT, recon_SIT, recon_SIRT

recon_fc_si_ts = zeros(n,n,S);
recon_si_ts = zeros(n,n,S);
orig_si_ts = zeros(n,n,S);
for iii = 1:S
    for i = 1:n
        for j = 1:n
%            if i < j
                recon_fc_si_M = zeros(numTasks,2);
                recon_si_M = zeros(numTasks,2);
                si_M = zeros(numTasks,2);
                for jj = 1:(numTasks)
                    recon_fc_si_M(jj,1) = recon_FC_SIT{jj}(i,j,iii);
                    recon_fc_si_M(jj,2) = recon_FC_SIRT{jj}(i,j,iii);
                    recon_si_M(jj,1) = recon_SIT{jj}(i,j,iii);
                    recon_si_M(jj,2) = recon_SIRT{jj}(i,j,iii);
                    si_M(jj,1) = orig_SIT{jj}(i,j,iii);
                    si_M(jj,2) = orig_SIRT{jj}(i,j,iii);
                end
                recon_fc_si_ts(i,j,iii) = icc21(recon_fc_si_M);
                recon_si_ts(i,j,iii) = icc21(recon_si_M);
                orig_si_ts(i,j,iii) = icc21(si_M);
%            end
        end
    end
end

save('TS_SI.mat', 'recon_fc_si_ts', 'recon_si_ts', 'orig_si_ts', '-v7.3');
