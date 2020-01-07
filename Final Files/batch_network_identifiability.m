% changing default fontsize
fontsize = 20;
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultAxesFontSize',fontsize-2);

set(0,'DefaultTextFontname','Times New Roman');
set(0,'DefaultTextFontSize',fontsize);

% Pick the file you want
load('FC_rest.mat');
%load('FC_language,mat');

n = size(FC_test,1); % Brain region pairs
numSubj = size(FC_test,3); % No. of subjects

configs.numRegions = n; % aparc2009 parcellation. Number of brain regions
configs.numSubj = numSubj; % Number of subjects
configs.numFCs = 2*numSubj; % number of connectomes (2 FCs per subject in this case)
configs.numEdges = nnz(mask_ut);
configs.numVisits = 2; % 2 visits per subject (test-retest)
configs.max_numPCs = configs.numFCs; % maximum number of PCs == data dimension
configs.symmetricnp = false; % True if the network property is symmetric i.e NP_ij == NP_ji
configs.stepPC = 5; % Steps of increase for PCA components reconstruction


% FC_test and FC_retest are each an n x n x numSubj data
disp("Search Information..")
network_property = @search_information; % Function handle for the network property that needs to be derived
[Idscore_FC_recon, Idscore_SI_recon, Idscore_SI_FC_recon] = f_network_identifiability(FC_test, FC_retest, configs, network_property);
disp("Mean First Passage Time..")
network_property = @mean_first_passage_time; % Function handle for the network property that needs to be derived
[~, Idscore_MFPT_recon, Idscore_MFPT_FC_recon] = f_network_identifiability(FC_test , FC_retest, configs, network_property);


PCA_comps_range = 1:configs.stepPC:configs.numFCs;

figure;
subplot(1,2,1);
plot(PCA_comps_range, Idscore_FC_recon, 'Color', [0.8392, 0.1529, 0.1569], 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4); hold on;
plot(PCA_comps_range, Idscore_SI_FC_recon, 'Color', [0.5804, 0.4039, 0.7412], 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
plot(PCA_comps_range, Idscore_SI_recon, 'Color', [0.5804, 0.4039, 0.7412], 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
axis square;
legend('Reconstructed FC', 'Search Information, Reconstructed FC', 'Reconstructed Search Information, Origial FC')
xlabel('# Principal Components', 'FontSize', labelsize); ylabel('Idiff (%)', 'FontSize', labelsize);
title('Search Information');

subplot(1,2,2);
plot(PCA_comps_range, Idscore_FC_recon, 'Color', [0.8392, 0.1529, 0.1569], 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4); hold on;
plot(PCA_comps_range, Idscore_MFPT_FC_recon, 'Color', [0.1725, 0.6275, 0.1725], 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
plot(PCA_comps_range, Idscore_MFPT_recon, 'Color', [0.1725, 0.6275, 0.1725], 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
axis square;
legend('Reconstructed FC', 'MFPT, Reconstructed FC', 'Reconstructed MFPT, Origial FC');
xlabel('# Principal Components', 'FontSize', labelsize); ylabel('Idiff (%)', 'FontSize', labelsize);
title('Mean First Passage Time');