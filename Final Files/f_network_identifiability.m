%{
Input:  configs are the configuration details for PCA
        FCorig_test_3D is the functional connectome test data
                       n x n x numSubj file where
                       n is the number of brain region pairs and numSubj is the total number of subjects
        FCorig_retest_3D is the functional connectome retest data
        find_network_property is the function handle that finds a network property of interest

Output: Idscore_FC_recon    Differential Identifiability score at different levels (no. of PCs) of reconstruction of functional connectomes
        Idscore_NP_recon    Differential Identifiability score at different levels (no. of PCs) of reconstruction of network property on original functional connectomes
        Idscore_NP_FC_recon Differential Identifiability score at different levels (no. of PCs) of reconstruction of network property on reconstructed functional connectomes

%}

function [Idscore_FC_recon, Idscore_NP_recon, Idscore_NP_FC_recon] = network_identifiability(FCorig_test_3D, FCorig_retest_3D, configs, network_property)

Test_index = 1:configs.numVisits:configs.numFCs; % change this and next line if you have a different FCs ordering
Retest_index = 2:configs.numVisits:configs.numFCs;
mask_ut = triu(true(configs.numRegions,configs.numRegions),1);

FCorig_test_2D = zeros(configs.numEdges, configs.numSubj); % Origninal data test
FCorig_retest_2D = zeros(configs.numEdges, configs.numSubj); % Origninal data retest


% Building matrix to do PCA decomposition
% Vectorizing the upper triangular data (FC's are symmetrical)

for subj = 1:configs.numSubj
    aux = FCorig_test_3D(:,:,subj);
    FCorig_test_2D(:,subj) = aux(mask_ut);
    aux = FCorig_retest_3D(:,:,subj);
    FCorig_retest_2D(:,subj) = aux(mask_ut);
end

NPorig_FCorig_test_3D = zeros(configs.numRegions,configs.numRegions,configs.numSubj); % Network Property on original FC test
NPorig_FCorig_retest_3D = zeros(configs.numRegions,configs.numRegions,configs.numSubj); % Network Property on original FC retest

if configs.symmetricnp    
    NPorig_FCorig_test_2D = zeros(configs.numEdges, configs.numSubj);
    NPorig_FCorig_retest_2D = zeros(configs.numEdges, configs.numSubj);
    for subj = 1:configs.numSubj
        NPorig_FCorig_test_3D(:,:,subj) = network_property(FCorig_test_3D(:,:,subj));
        NPorig_FCorig_retest_3D(:,:,subj) = network_property(FCorig_retest_3D(:,:,subj));
        % Making a vector and reconstructing
        aux = NPorig_FCorig_test_3D(:,:,subj);
        NPorig_FCorig_test_2D(:,subj) = aux(mask_ut);
        aux = NPorig_FCorig_retest_3D(:,:,subj);
        NPorig_FCorig_retest_2D(:,subj) = aux(mask_ut);
    end
else
    NPorig_FCorig_test_2D = zeros(configs.numRegions*configs.numRegions, configs.numSubj); % Network Property on original data test
    NPorig_FCorig_retest_2D = zeros(configs.numRegions*configs.numRegions, configs.numSubj); % Network Property on original data retest
 
    for subj = 1:configs.numSubj
        NPorig_FCorig_test_3D(:,:,subj) = network_property(FCorig_test_3D(:,:,subj));
        NPorig_FCorig_retest_3D(:,:,subj) = network_property(FCorig_retest_3D(:,:,subj));

        % Making a vector and reconstructing
        aux = NPorig_FCorig_test_3D(:,:,subj);
        NPorig_FCorig_test_2D(:,subj) = aux(:);
        aux = NPorig_FCorig_retest_3D(:,:,subj);
        NPorig_FCorig_retest_2D(:,subj) = aux(:);
    end
end
    
% Interlacing alternate test and retest for the PCA decomposition
FCorig_2D = FCorig_retest_2D(:,[1;1]*(1:size(FCorig_retest_2D,2)));
FCorig_2D(:,1:2:end) = FCorig_test_2D;

% Interlacing alternate test and retest for the PCA decomposition
NPorig_FCorig_2D = NPorig_FCorig_retest_2D(:,[1;1]*(1:size(NPorig_FCorig_retest_2D,2)));
NPorig_FCorig_2D(:,1:2:end) = NPorig_FCorig_test_2D;

disp("Decomposing into components");
% Decomposing the FC and NP into components
[FC_modes, projected_FC_modes] = pca(FCorig_2D, 'NumComponents', configs.max_numPCs);
[NP_modes, projected_NP_modes] = pca(NPorig_FCorig_2D, 'NumComponents', configs.max_numPCs);


PCvector = 2:configs.stepPC:configs.numFCs; %vector with all number of PCs to be evaluated
Idscore_FC_recon = zeros(length(PCvector),1);
Idscore_NP_recon = zeros(length(PCvector),1);
Idscore_NP_FC_recon = zeros(length(PCvector),1);
disp("Reconstructing");
counter = 1;
% Reconstructing from the components
for i = PCvector % For each number of components from 2 to dimensionality of data
    % Reconstructing functional connectomes
    FCrecon_2D = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
    FCrecon_2D = bsxfun(@plus, FCrecon_2D, mean(FCorig_2D));
    FCrecon_test_2D = FCrecon_2D(:, Test_index);
    FCrecon_retest_2D = FCrecon_2D(:, Retest_index);
    
    % Reconstructing network property found on original functional connectomes
    NPrecon_FCorig_2D = projected_NP_modes(:,1:i) * NP_modes(:,1:i)';
    NPrecon_FCorig_2D = bsxfun(@plus, NPrecon_FCorig_2D, mean(NPorig_FCorig_2D));
    NPrecon_FCorig_test_2D = NPrecon_FCorig_2D(:, Test_index);
    NPrecon_FCorig_retest_2D = NPrecon_FCorig_2D(:, Retest_index);
    
    aux = zeros(configs.numRegions,configs.numRegions);
    FCrecon_test_3D = zeros(configs.numRegions,configs.numRegions,configs.numSubj);
    FCrecon_retest_3D = zeros(configs.numRegions,configs.numRegions,configs.numSubj);
    
    NP_FCrecon_test_3D = zeros(configs.numRegions,configs.numRegions,configs.numSubj);
    NP_FCrecon_retest_3D = zeros(configs.numRegions,configs.numRegions,configs.numSubj);
    
    % Finding network properties on reconstructed functional connectomes
    if configs.symmetricnp
        NP_FCrecon_test_2D = zeros(configs.numEdges, configs.numSubj);
        NP_FCrecon_retest_2D = zeros(configs.numEdges, configs.numSubjs);
        for subj = 1:configs.numSubj % For each subject data available
            aux = zeros(configs.numRegions,configs.numRegions);
            aux(mask_ut) = FCrecon_test_2D(:,subj);
            FCrecon_test_3D(:,:,subj) = aux + aux';
            aux(mask_ut) = FCrecon_retest_2D(:,subj);
            FCrecon_retest_3D(:,:,subj) = aux + aux';

            NP_FCrecon_test_3D(:,:,subj) = network_property(FCrecon_test_3D(:,:,subj));
            NP_FCrecon_retest_3D(:,:,subj) = network_property(FCrecon_retest_3D(:,:,subj));

            aux = NP_FCrecon_test_3D(:,:,subj);
            NP_FCrecon_test_2D(:,subj) = aux(mask_ut);
            aux = NP_FCrecon_retest_3D(:,:,subj);
            NP_FCrecon_retest_2D(:,subj) = aux(mask_ut);
        end
    else
        NP_FCrecon_test_2D = zeros(configs.numRegions*configs.numRegions, configs.numSubj);
        NP_FCrecon_retest_2D = zeros(configs.numRegions*configs.numRegions, configs.numSubj);
        for subj = 1:configs.numSubj % For each subject data available
            aux = zeros(configs.numRegions,configs.numRegions);
            aux(mask_ut) = FCrecon_test_2D(:,subj);
            FCrecon_test_3D(:,:,subj) = aux + aux';
            aux(mask_ut) = FCrecon_retest_2D(:,subj);
            FCrecon_retest_3D(:,:,subj) = aux + aux';

            NP_FCrecon_test_3D(:,:,subj) = network_property(FCrecon_test_3D(:,:,subj));
            NP_FCrecon_retest_3D(:,:,subj) = network_property(FCrecon_retest_3D(:,:,subj));

            aux = NP_FCrecon_test_3D(:,:,subj);
            NP_FCrecon_test_2D(:,subj) = aux(:);
            aux = NP_FCrecon_retest_3D(:,:,subj);
            NP_FCrecon_retest_2D(:,subj) = aux(:);
        end
    end
    
    % Compute differential identifiability matrices and scores
    mask_diag = logical(eye(configs.numSubj));
    
    Ident_FC_recon = corr(FCrecon_retest_2D,FCrecon_test_2D);
    Iself_FC_recon = mean(Ident_FC_recon(mask_diag));
    Iothers_FC_recon = mean(Ident_FC_recon(~mask_diag));
    Idscore_FC_recon(counter) = (Iself_FC_recon - Iothers_FC_recon) * 100;
    
    Ident_NP_recon = corr(NPrecon_FCorig_retest_2D,NPrecon_FCorig_test_2D);
    Iself_NP_recon = mean(Ident_NP_recon(mask_diag));
    Iothers_NP_recon = mean(Ident_NP_recon(~mask_diag));
    Idscore_NP_recon(counter) = (Iself_NP_recon - Iothers_NP_recon) * 100;
    
    Ident_NP_FC_recon = corr(NP_FCrecon_retest_2D,NP_FCrecon_test_2D);
    Iself_NP_FC_recon = mean(Ident_NP_FC_recon(mask_diag));
    Iothers_NP_FC_recon = mean(Ident_NP_FC_recon(~mask_diag));
    Idscore_NP_FC_recon(counter) = (Iself_NP_FC_recon - Iothers_NP_FC_recon) * 100;
    counter = counter + 1;
end

