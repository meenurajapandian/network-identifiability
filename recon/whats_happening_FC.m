function whats_happening_FC(ii)
load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
%load(fullfile(pwd,'data','HCP_100S_FC_glasser_April2018.mat'));
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

% Building matrix to do PCA decomposition
for iii = 1:S
aux = orig_FCT(:,:,iii);
orig_matrix(:,(2*iii)-1) = aux(mask_ut);
aux = orig_FCRT(:,:,iii);
orig_matrix(:,2*iii) = aux(mask_ut);
end

orig_matrix_test = orig_matrix(:, Test_index);
orig_matrix_retest = orig_matrix(:, Retest_index);

[~,~,Idiff_recon,~,~,~,~,~]  = f_PCA_identifiability(orig_matrix,Test_index,Retest_index,configs);


% Decomposing and also reconstructing
% v = VideoWriter('emotion_W.avi');
% v.FrameRate = 25; % control frame-rate
% open(v);
% h=figure(2); %creates the figure and the handle to the figure


% close(v)
save(sprintf('Ident_%d.mat',ii), 'Idiff_recon');