PCA_comps_range = 1:200;
T = 18;
S = 100;
tasks_list = {'EMOTION' 'GAMBLING' 'LANGUAGE' 'MOTOR' 'RELATIONAL' 'SOCIAL' 'WM' 'REST1' 'REST2'};
numTasks = length(tasks_list);

path2ident = '/home/mrajapan/Documents/MATLAB/data/Identifiability';
path2data = '/home/mrajapan/Documents/MATLAB/data';

orig_col = [214 39 40; 31 119 180; 44 160 44; 227 119 194; 34 187 186; 0 229 28; 148 103 189; 255 127 14; 255 127 14;];
col = orig_col ./ 255;

labelsize = 20;
ticksize = 18;
titlesize = 22;

%% Task Sensitivity

load(fullfile(path2data,'TS_SP.mat'));
for iii = 1:S
    orig_sp_ts(:,:,iii) = orig_sp_ts(:,:,iii) + orig_sp_ts(:,:,iii)';
    recon_fc_sp_ts(:,:,iii) = recon_fc_sp_ts(:,:,iii) + recon_fc_sp_ts(:,:,iii)';
    recon_sp_ts(:,:,iii) = recon_sp_ts(:,:,iii) + recon_sp_ts(:,:,iii)';
end

m_orig_sp_ts = mean(orig_sp_ts,3);
m_recon_fc_sp_ts = mean(recon_fc_sp_ts,3);
m_recon_sp_ts = mean(recon_sp_ts,3);

figure;
subplot(3,5,1); plot(m_orig_sp_ts(:),m_recon_fc_sp_ts(:),'.', 'Color', col(2,:));
axis square; box on; xlim([0 1]); ylim([0 1]);
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('ICC - SPL(FC)', 'FontSize',labelsize); ylabel('ICC - SPL(PCA(FC))', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,6); plot(m_orig_sp_ts(:),m_recon_sp_ts(:),'.', 'Color', col(2,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('ICC - SPL(FC)', 'FontSize',labelsize); ylabel('ICC - PCA(SPL(FC))', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,11); plot(m_recon_fc_sp_ts(:),m_recon_sp_ts(:),'.', 'Color', col(2,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('ICC - SPL(PCA(FC))', 'FontSize',labelsize); ylabel('ICC - PCA(SPL(FC))', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

clear orig_sp_ts recon_fc_sp_ts recon_sp_ts m_orig_mfpt_ts m_recon_fc_sp_ts m_recon_sp_ts



load(fullfile(path2data,'TS_MFPT.mat'));

m_orig_mfpt_ts = mean(orig_mfpt_ts,3);
m_recon_fc_mfpt_ts = mean(recon_fc_mfpt_ts,3);
m_recon_mfpt_ts = mean(recon_mfpt_ts,3);

subplot(3,5,2); plot(m_orig_mfpt_ts(:),m_recon_fc_mfpt_ts(:),'.', 'Color', col(3,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (MFPT, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (MFPT, PCA on FC)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,7); plot(m_orig_mfpt_ts(:),m_recon_mfpt_ts(:),'.', 'Color', col(3,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (MFPT, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (MFPT, PCA on MFPT)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,12); plot(m_recon_fc_mfpt_ts(:),m_recon_mfpt_ts(:),'.', 'Color', col(3,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (MFPT, PCA on FC)', 'FontSize',labelsize); ylabel('Task ICC (MFPT, PCA on MFPT)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;

clear orig_mfpt_ts recon_fc_mfpt_ts recon_mfpt_ts m_orig_mfpt_ts m_recon_fc_mfpt_ts m_recon_mfpt_ts



load(fullfile(path2data,'TS_SI.mat'));

m_orig_si_ts = mean(orig_si_ts,3);
m_recon_fc_si_ts = mean(recon_fc_si_ts,3);
m_recon_si_ts = mean(recon_si_ts,3);

subplot(3,5,3); plot(m_orig_si_ts(:),m_recon_fc_si_ts(:),'.', 'Color', col(7,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (SI, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (SI, PCA on FC)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,8); plot(m_orig_si_ts(:),m_recon_si_ts(:),'.', 'Color', col(7,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (SI, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (SI, PCA on SI)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,13); plot(m_recon_fc_si_ts(:),m_recon_si_ts(:),'.', 'Color', col(7,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (SI, PCA on FC)', 'FontSize',labelsize); ylabel('Task ICC (SI, PCA on SI)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

clear orig_si_ts recon_fc_si_ts recon_si_ts m_orig_si_ts m_recon_fc_si_ts m_recon_si_ts



load(fullfile(path2data,'TS_W_1.mat'));
load(fullfile(path2data,'TS_W_2.mat'));
load(fullfile(path2data,'TS_W_3.mat'));

m_orig_w_ts = mean(orig_w_ts,3);
m_recon_fc_w_ts = mean(recon_fc_w_ts,3);
m_recon_w_ts = mean(recon_w_ts,3);

subplot(3,5,4); plot(m_orig_w_ts(:),m_recon_fc_w_ts(:),'.', 'Color', col(4,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (W, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (W, PCA on FC)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,9); plot(m_orig_w_ts(:),m_recon_w_ts(:),'.', 'Color', col(4,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (W, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (W, PCA on W)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,14); plot(m_recon_fc_w_ts(:),m_recon_w_ts(:),'.', 'Color', col(4,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (W, PCA on FC)', 'FontSize',labelsize); ylabel('Task ICC (W, PCA on W)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

clear orig_w_ts recon_fc_w_ts recon_w_ts m_orig_w_ts m_recon_fc_w_ts m_recon_w_ts


load(fullfile(path2data,'TS_EC.mat'));
for iii = 1:S
    orig_ec_ts(:,:,iii) = orig_ec_ts(:,:,iii) + orig_ec_ts(:,:,iii)';
    recon_fc_ec_ts(:,:,iii) = recon_fc_ec_ts(:,:,iii) + recon_fc_ec_ts(:,:,iii)';
    recon_ec_ts(:,:,iii) = recon_ec_ts(:,:,iii) + recon_ec_ts(:,:,iii)';
end

m_orig_ec_ts = mean(orig_ec_ts,3);
m_recon_fc_ec_ts = mean(recon_fc_ec_ts,3);
m_recon_ec_ts = mean(recon_ec_ts,3);

subplot(3,5,5); plot(m_orig_ec_ts(:),m_recon_fc_ec_ts(:),'.', 'Color', col(8,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (EC, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (EC, PCA on FC)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,10); plot(m_orig_ec_ts(:),m_recon_ec_ts(:),'.', 'Color', col(8,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (EC, No PCA)', 'FontSize',labelsize); ylabel('Task ICC (EC, PCA on EC)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

subplot(3,5,15); plot(m_recon_fc_ec_ts(:),m_recon_ec_ts(:),'.', 'Color', col(8,:));
axis square; box on; xlim([0 1]); ylim([0 1]); 
handle = gca; handle.LineWidth = 1.7; handle.FontSize = ticksize;
xlabel('Task ICC (EC, PCA on FC)', 'FontSize',labelsize); ylabel('Task ICC (EC, PCA on EC)', 'FontSize',labelsize);
hold on; plot(0:0.5:1,0:0.5:1, 'Color', [0, 0, 0], 'LineWidth', 1.5);
handle = gca; handle.LineWidth = 1.5;
xticks([0:0.2:1]); yticks([0:0.2:1]);

clear orig_ec_ts recon_fc_ec_ts recon_ec_ts m_orig_ec_ts m_recon_fc_ec_ts m_recon_ec_ts

%% Node Pair Properties - only pca on network
figure;
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_%d.mat',(2*jj)-1)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idiff_recon, 'Color', col(1,:), 'LineStyle', '-', 'LineWidth',4,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SP_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_SP_recon, 'Color', col(2,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SI_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_SI_recon, 'Color', col(7,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_W_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_W_recon, 'Color', col(4,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_EC_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_EC_recon, 'Color', col(8,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end



% for jj = 1:(numTasks)
%     load(fullfile(path2ident, sprintf('Ident_SS_%d.mat',jj)));
%     subplot(2,5,jj); hold on;
%     plot(PCA_comps_range, Idscore_SS_FC_recon(:,1), 'Color', col(5,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%     plot(PCA_comps_range, Idscore_SS_recon(:,1), 'Color', col(5,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
%     hold off;
% end
% 
% for jj = 1:(numTasks)
%     load(fullfile(path2ident, sprintf('Ident_TS_%d.mat',jj)));
%     subplot(2,5,jj); hold on;
%     plot(PCA_comps_range, Idscore_TS_FC_recon, 'Color', col(6,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%     plot(PCA_comps_range, Idscore_TS_recon, 'Color', col(6,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
%     hold off;
% end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_MFPT_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_MFPT_recon, 'Color', col(3,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    ylim([0 35]);
    handle = gca;
    handle.LineWidth = 1.7;
    handle.FontSize = ticksize;
    xlabel('# Principal Components', 'FontSize', labelsize); ylabel('Idiff (%)', 'FontSize', labelsize);
    hold off;
    axis square
    box on;
    title(sprintf('%s', tasks_list{jj}), 'FontSize', titlesize);
    if (jj > 5) 
    handle.Position = handle.Position + [0 0.11 0 0];
    end
    xticks([0:50:200]);
end

%hSub = subplot(2,5,10); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');

%legend('Functional Connectivity (FC)', 'Shortest Path (PCA on FC)','Shortest Path (PCA on Network Property)','Driftness (PCA on FC)','Driftness (PCA on Network Property)','Communicability (PCA on FC)','Communicability (PCA on Network Property)','Search Information (PCA on FC)','Search Information (PCA on Network Property)','Source Similarity (PCA on FC)','Source Similarity (PCA on Network Property)','Target Similarity (PCA on FC)','Target Similarity (PCA on Network Property)', 'MFPT (PCA on FC)','MFPT (PCA on Network Property)');

subplot(2,5,10);
plot(NaN, NaN, 'Color', col(1,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4); hold on;

plot(NaN, NaN, 'Color', col(2,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(7,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(3,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(4,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(8,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

 axis square;
axis off;
lgd = legend('Functional Connectivity (FC)', 'Shortest Path Length (SPL)','Search Information (SI)', 'Mean First Passage Time (MFPT)','Driftness (W)','Communicability (C)','Location','northwest');
%lgd = legend('Functional Connectivity (FC)', 'Shortest Path (PCA on FC)','Shortest Path (PCA on Network Property)','Driftness (PCA on FC)','Driftness (PCA on Network Property)','Communicability (PCA on FC)','Communicability (PCA on Network Property)','Search Information (PCA on FC)','Search Information (PCA on Network Property)', 'MFPT (PCA on FC)','MFPT (PCA on Network Property)','Location','northwest');
legend boxoff
%set(lgd,'position',[0.7811,0.392,0.1559,0.1544])
lgd.FontSize = 20;

%% Node Pair Properties - only pca on fc
figure;
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_%d.mat',(2*jj)-1)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idiff_recon, 'Color', col(1,:), 'LineStyle', '-', 'LineWidth',4,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SP_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_SP_FC_recon, 'Color', col(2,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    %plot(PCA_comps_range, Idscore_SP_recon, 'Color', col(2,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SI_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_SI_FC_recon, 'Color', col(7,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    %plot(PCA_comps_range, Idscore_SI_recon, 'Color', col(7,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_W_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_W_FC_recon, 'Color', col(4,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    %plot(PCA_comps_range, Idscore_W_recon, 'Color', col(4,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_EC_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_EC_FC_recon, 'Color', col(8,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    %plot(PCA_comps_range, Idscore_EC_recon, 'Color', col(8,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end



% for jj = 1:(numTasks)
%     load(fullfile(path2ident, sprintf('Ident_SS_%d.mat',jj)));
%     subplot(2,5,jj); hold on;
%     plot(PCA_comps_range, Idscore_SS_FC_recon(:,1), 'Color', col(5,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%     plot(PCA_comps_range, Idscore_SS_recon(:,1), 'Color', col(5,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
%     hold off;
% end
% 
% for jj = 1:(numTasks)
%     load(fullfile(path2ident, sprintf('Ident_TS_%d.mat',jj)));
%     subplot(2,5,jj); hold on;
%     plot(PCA_comps_range, Idscore_TS_FC_recon, 'Color', col(6,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%     plot(PCA_comps_range, Idscore_TS_recon, 'Color', col(6,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
%     hold off;
% end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_MFPT_%d.mat',jj)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_MFPT_FC_recon, 'Color', col(3,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    %plot(PCA_comps_range, Idscore_MFPT_recon, 'Color', col(3,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    ylim([0 35]);
    handle = gca;
    handle.LineWidth = 1.7;
    handle.FontSize = ticksize;
    xlabel('# Principal Components', 'FontSize', labelsize); ylabel('Idiff (%)', 'FontSize', labelsize);
    hold off;
    axis square
    box on;
    title(sprintf('%s', tasks_list{jj}), 'FontSize', titlesize);
    if (jj > 5) 
    handle.Position = handle.Position + [0 0.11 0 0];
    end
    xticks([0:50:200]);
end

%hSub = subplot(2,5,10); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');

%legend('Functional Connectivity (FC)', 'Shortest Path (PCA on FC)','Shortest Path (PCA on Network Property)','Driftness (PCA on FC)','Driftness (PCA on Network Property)','Communicability (PCA on FC)','Communicability (PCA on Network Property)','Search Information (PCA on FC)','Search Information (PCA on Network Property)','Source Similarity (PCA on FC)','Source Similarity (PCA on Network Property)','Target Similarity (PCA on FC)','Target Similarity (PCA on Network Property)', 'MFPT (PCA on FC)','MFPT (PCA on Network Property)');

subplot(2,5,10);
plot(NaN, NaN, 'Color', col(1,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4); hold on;

plot(NaN, NaN, 'Color', col(2,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%plot(NaN, NaN, 'Color', col(2,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(7,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%plot(NaN, NaN, 'Color', col(7,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(3,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%plot(NaN, NaN, 'Color', col(3,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4); axis square;

plot(NaN, NaN, 'Color', col(4,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%plot(NaN, NaN, 'Color', col(4,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(8,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
%plot(NaN, NaN, 'Color', col(8,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

axis off;
lgd = legend('Functional Connectivity (FC)', 'Shortest Path Length (SPL)', 'Search Information (SI)', 'Mean First Passage Time (MFPT)', 'Driftness (W)','Communicability (C)','Location','northwest');
%lgd = legend('Functional Connectivity (FC)', 'Shortest Path (PCA on FC)','Shortest Path (PCA on Network Property)','Driftness (PCA on FC)','Driftness (PCA on Network Property)','Communicability (PCA on FC)','Communicability (PCA on Network Property)','Search Information (PCA on FC)','Search Information (PCA on Network Property)', 'MFPT (PCA on FC)','MFPT (PCA on Network Property)','Location','northwest');
legend boxoff
%set(lgd,'position',[0.7811,0.392,0.1559,0.1544])
lgd.FontSize = 20;
%% Node Properties

figure;
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_%d.mat',(2*jj)-1)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idiff_recon, 'Color', col(1,:), 'LineStyle', '-', 'LineWidth',4,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_D_%d.mat',(2*jj)-1)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_D_FC_recon(:,1), 'Color', col(5,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    plot(PCA_comps_range, Idscore_D_recon(:,1), 'Color', col(5,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_C_%d.mat',(2*jj)-1)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_C_FC_recon(:,1), 'Color', col(6,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    plot(PCA_comps_range, Idscore_C_recon(:,1), 'Color', col(6,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    hold off;
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_B_%d.mat',(2*jj)-1)));
    subplot(2,5,jj); hold on;
    plot(PCA_comps_range, Idscore_B_FC_recon(:,1), 'Color', col(9,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
    plot(PCA_comps_range, Idscore_B_recon(:,1), 'Color', col(9,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);
    ylim([0 37]);
    handle = gca;
    handle.LineWidth = 1.7;
    handle.FontSize = ticksize;
    xlabel('PCA # Components', 'FontSize', labelsize); ylabel('Idiff (%)', 'FontSize', labelsize);
    hold off;
    axis square
    box on;
    title(sprintf('%s', tasks_list{jj}), 'FontSize', titlesize);
    xticks([0:50:200]);
end


%hSub = subplot(2,5,10); plot(1, nan, 1, nan, 'r'); set(hSub, 'Visible', 'off');

%legend('Functional Connectivity (FC)', 'Shortest Path (PCA on FC)','Shortest Path (PCA on Network Property)','Driftness (PCA on FC)','Driftness (PCA on Network Property)','Communicability (PCA on FC)','Communicability (PCA on Network Property)','Search Information (PCA on FC)','Search Information (PCA on Network Property)','Source Similarity (PCA on FC)','Source Similarity (PCA on Network Property)','Target Similarity (PCA on FC)','Target Similarity (PCA on Network Property)', 'MFPT (PCA on FC)','MFPT (PCA on Network Property)');

subplot(2,5,10);
plot(NaN, NaN, 'Color', col(1,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4); hold on;

plot(NaN, NaN, 'Color', col(5,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
plot(NaN, NaN, 'Color', col(5,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(6,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
plot(NaN, NaN, 'Color', col(6,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4);

plot(NaN, NaN, 'Color', col(9,:), 'LineStyle', '-', 'LineWidth',2,'MarkerSize',4);
plot(NaN, NaN, 'Color', col(9,:), 'LineStyle', '--', 'LineWidth',2,'MarkerSize',4); axis square;
axis off;
lgd = legend('Functional Connectivity (FC)', 'Degree (PCA on FC)','Degree (PCA on Network Property)','Clustering Coeff (PCA on FC)','Clustering Coeff (PCA on Network Property)','Betweenness (PCA on FC)','Betweenness (PCA on Network Property)');

legend boxoff
set(lgd,'position',[0.7839, 0.4029, 0.1468, 0.0990])
lgd.FontSize = 20;

%% MFPT AND SI

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_%d.mat',(2*jj)-1)));
    fc(jj,:) = Idiff_recon(2:200);
end

o_m_fc = mean(fc,1);
o_q1_fc = quantile(fc, 0.05, 1);
o_q2_fc = quantile(fc, 0.95, 1);


for jj = 1:(numTasks-1)
    load(fullfile(path2ident, sprintf('Ident_MFPT_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_MFPT_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_MFPT_recon(2:200,:)';
end

m_fc = mean(Id_fc,1);
q1_fc = quantile(Id_fc, 0.05, 1);
q2_fc = quantile(Id_fc, 0.95, 1);
m = mean(Id,1);
q1 = quantile(Id, 0.05, 1);
q2 = quantile(Id, 0.95, 1);

figure;
subplot(1,2,1);
h = fill([2:200,fliplr(2:200)],[q1,fliplr(q2)],col(3,:),'EdgeColor','none'); hold on;
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hp = findobj(h,'type','patch');
hh = hatchfill(hp,'single', 75, 10, col(3,:)); % Variables of hatchfill (A,STYL,ANGLE,SPACING,FACECOL)
set(hh,'color',[1 1 1 0.85],'linewidth',5);
set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h = fill([2:200,fliplr(2:200)],[q1_fc,fliplr(q2_fc)],col(3,:),'EdgeColor','none');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
alpha(0.25); % specifies the transparency of fill
plot(2:200,m_fc,'Color',col(3,:),'LineWidth',3);
plot(2:200,m,'Color',col(3,:),'LineStyle','--','LineWidth',3);
plot(2:200,o_m_fc,'Color',col(1,:),'LineStyle',':','LineWidth',4); hold on;
xlim([2,200]);
ylim([0,35]);
axis square;
hold off;
box on;
title('Mean First Passage Time (MFPT)', 'FontSize', titlesize);
set(gca,'FontSize',ticksize+5);
set(gca,'LineWidth',2);
xlabel('# Principal Components', 'FontSize', labelsize); ylabel('Idiff (%)', 'FontSize', labelsize);
xticks([25:25:200]);


for jj = 1:(numTasks-1)
    load(fullfile(path2ident, sprintf('Ident_SI_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_SI_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_SI_recon(2:200,:)';
end

m_fc = mean(Id_fc,1);
q1_fc = quantile(Id_fc, 0.05, 1);
q2_fc = quantile(Id_fc, 0.95, 1);
m = mean(Id,1);
q1 = quantile(Id, 0.05, 1);
q2 = quantile(Id, 0.95, 1);

%figure;
subplot(1,2,2);
h = fill([2:200,fliplr(2:200)],[q1,fliplr(q2)],col(7,:),'EdgeColor','none'); hold on;
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hp = findobj(h,'type','patch');
hh = hatchfill(hp,'single',50,10, col(7,:)); % Variables of hatchfill (A,STYL,ANGLE,SPACING,FACECOL)
set(hh,'color',[1 1 1 0.85],'linewidth',5);
set(get(get(hh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
h = fill([2:200,fliplr(2:200)],[q1_fc,fliplr(q2_fc)],col(7,:),'EdgeColor','none');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% specifies the transparency of fill
alpha(0.25);
plot(2:200,m_fc,'Color',col(7,:),'LineWidth',3);
plot(2:200,m,'Color',col(7,:),'LineStyle','--','LineWidth',3);
plot(2:200,o_m_fc,'Color',col(1,:),'LineStyle',':','LineWidth',4);
xlim([2,200]);
ylim([0,35]);   
axis square;
hold off;
box on;
title('Search Information (SI)', 'FontSize', titlesize);
set(gca,'FontSize',ticksize+5);
set(gca,'LineWidth',2);
xlabel('# Principal Components', 'FontSize', labelsize); ylabel('Idiff (%)', 'FontSize', labelsize);
xticks([25:25:200]);

%%
figure;
%96 All Language = 5 & 6
load('/home/mrajapan/Documents/MATLAB/data/HCP_100S_FC_glasser_April2018.mat');
ii = 5;
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

for iii = 1:S
aux = orig_FCT(:,:,iii);
orig_matrix(:,(2*iii)-1) = aux(mask_ut);
aux = orig_FCRT(:,:,iii);
orig_matrix(:,2*iii) = aux(mask_ut);
end

[FC_modes, projected_FC_modes] = pca(orig_matrix, 'NumComponents', configs.max_numPCs);

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

orig_SI_test = orig_SI(:, Test_index);
orig_SI_retest = orig_SI(:, Retest_index);

[SI_modes, projected_SI_modes] = pca(orig_SI, 'NumComponents', configs.max_numPCs);

i = 96;
recon_matrix = projected_FC_modes(:,1:i) * FC_modes(:,1:i)';
recon_matrix = bsxfun(@plus, recon_matrix, mean(orig_matrix));
recon_matrix_test = recon_matrix(:, Test_index);
recon_matrix_retest = recon_matrix(:, Retest_index);

recon_SI = projected_SI_modes(:,1:i) * SI_modes(:,1:i)';
recon_SI = bsxfun(@plus, recon_SI, mean(orig_SI));
recon_SI_test = recon_SI(:, Test_index);
recon_SI_retest = recon_SI(:, Retest_index);

Ident_FC_recon = corr(recon_matrix_test, recon_matrix_retest);
Ident_SI_recon = corr(recon_SI_retest,recon_SI_test);

recon_matrix_test(recon_matrix_test < eps) = eps;
recon_matrix_retest(recon_matrix_retest < eps) = eps;

aux = zeros(n,n);
recon_FCT = zeros(n,n,S);
recon_FCRT = zeros(n,n,S);
recon_FC_MFPTT = zeros(n,n,S);
recon_FC_MFPTRT = zeros(n,n,S);

for iii = 1:S
aux = zeros(n,n);
aux(mask_ut) = recon_matrix_test(:,iii);
recon_FCT(:,:,iii) = aux + aux';
aux(mask_ut) = recon_matrix_retest(:,iii);
recon_FCRT(:,:,iii) = aux + aux';

[~, recon_FC_MFPTT(:,:,iii), ~] = f_mfpt(recon_FCT(:,:,iii));
[~, recon_FC_MFPTRT(:,:,iii), ~] = f_mfpt(recon_FCRT(:,:,iii));

aux = zeros(n,n);
aux = recon_FC_MFPTT(:,:,iii);
recon_FC_MFPT_test(:,iii) = aux(:);
aux = recon_FC_MFPTRT(:,:,iii);
recon_FC_MFPT_retest(:,iii) = aux(:);

end

Ident_MFPT_FC_recon = corr(recon_FC_MFPT_retest,recon_FC_MFPT_test);

subplot(1,3,1);
imagesc(Ident_MFPT_FC_recon); axis square;
title('Idiff of Language MFPT(PCA(FC))', 'FontSize', titlesize);
ylabel('Subject Retest', 'FontSize', labelsize+2); xlabel('Subject Test', 'FontSize', labelsize+2);
xticks([]); yticks([]);

subplot(1,3,2);
imagesc(Ident_FC_recon); axis square;
title('Idiff of Language PCA(FC)', 'FontSize', titlesize);
ylabel('Subject Retest', 'FontSize', labelsize+2); xlabel('Subject Test', 'FontSize', labelsize+2);
xticks([]); yticks([]);

subplot(1,3,3);
imagesc(Ident_SI_recon); axis square;
title('Idiff of Language PCA(SI(FC))', 'FontSize', titlesize);
ylabel('Subject Retest', 'FontSize', labelsize+2); xlabel('Subject Test', 'FontSize', labelsize+2);
xticks([]); yticks([]);

%% Heat Maps Horizontal
figure;
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SP_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_SP_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_SP_recon(2:200,:)';
end

subplot(2,5,1); imagesc(Id_fc); axis square; caxis([0, 26]); colormap(gca,flipud(make_color_map(orig_col(2,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Shortest Path', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
%cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(2,5,6); imagesc(Id); axis square; caxis([0, 26]); colormap(gca,flipud(make_color_map(orig_col(2,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
%handle.Position = handle.Position + [-0.3 0 0 0];
%title('Shortest Path', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_MFPT_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_MFPT_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_MFPT_recon(2:200,:)';
end

subplot(2,5,2); imagesc(Id_fc); axis square; caxis([0, 27]); colormap(gca,flipud(make_color_map(orig_col(3,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Mean First Passage Time', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
%cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(2,5,7); imagesc(Id); axis square; caxis([0, 27]); colormap(gca,flipud(make_color_map(orig_col(3,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
%handle.Position = handle.Position + [-0.3 0 0 0];
%title('Mean First Passage Time', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SI_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_SI_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_SI_recon(2:200,:)';
end

subplot(2,5,3); imagesc(Id_fc); axis square; caxis([0, 32]); colormap(gca,flipud(make_color_map(orig_col(7,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Search Information', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
%cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(2,5,8); imagesc(Id); axis square; caxis([0, 32]); colormap(gca,flipud(make_color_map(orig_col(7,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
%handle.Position = handle.Position + [-0.3 0 0 0];
%title('Search Information', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_W_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_W_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_W_recon(2:200,:)';
end

subplot(2,5,4); imagesc(Id_fc); axis square; caxis([0, 25]); colormap(gca,flipud(make_color_map(orig_col(4,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Driftness', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
%cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(2,5,9); imagesc(Id); axis square; caxis([0, 25]); colormap(gca,flipud(make_color_map(orig_col(4,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
%handle.Position = handle.Position + [-0.3 0 0 0];
%title('Driftness', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_EC_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_EC_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_EC_recon(2:200,:)';
end

subplot(2,5,5); imagesc(Id_fc); axis square; caxis([0, 30]); colormap(gca,flipud(make_color_map(orig_col(8,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Communicability', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
xlabel('PCA # Components');
cb=colorbar;
%cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(2,5,10); imagesc(Id); axis square; caxis([0, 30]); colormap(gca,flipud(make_color_map(orig_col(8,:),32))); 
box on;
%title('Communicability', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
xlabel('PCA # Components');
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
%handle.Position = handle.Position + [-0.3 0 0 0];

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

%% Heat Maps Vertical
figure;
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SP_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_SP_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_SP_recon(2:200,:)';
end

subplot(5,2,1); imagesc(Id_fc); axis square; caxis([0, 26]); colormap(gca,flipud(make_color_map(orig_col(2,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Shortest Path', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(5,2,2); imagesc(Id); axis square; caxis([0, 26]); colormap(gca,flipud(make_color_map(orig_col(2,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
handle.Position = handle.Position + [-0.3 0 0 0];
%title('Shortest Path', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_MFPT_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_MFPT_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_MFPT_recon(2:200,:)';
end

subplot(5,2,5); imagesc(Id_fc); axis square; caxis([0, 27]); colormap(gca,flipud(make_color_map(orig_col(3,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Mean First Passage Time', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(5,2,6); imagesc(Id); axis square; caxis([0, 27]); colormap(gca,flipud(make_color_map(orig_col(3,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
handle.Position = handle.Position + [-0.3 0 0 0];
%title('Mean First Passage Time', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SI_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_SI_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_SI_recon(2:200,:)';
end

subplot(5,2,3); imagesc(Id_fc); axis square; caxis([0, 32]); colormap(gca,flipud(make_color_map(orig_col(7,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Search Information', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(5,2,4); imagesc(Id); axis square; caxis([0, 32]); colormap(gca,flipud(make_color_map(orig_col(7,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
handle.Position = handle.Position + [-0.3 0 0 0];
%title('Search Information', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_W_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_W_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_W_recon(2:200,:)';
end

subplot(5,2,7); imagesc(Id_fc); axis square; caxis([0, 25]); colormap(gca,flipud(make_color_map(orig_col(4,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Driftness', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
cb=colorbar;
cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(5,2,8); imagesc(Id); axis square; caxis([0, 25]); colormap(gca,flipud(make_color_map(orig_col(4,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
handle.Position = handle.Position + [-0.3 0 0 0];
%title('Driftness', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_EC_%d.mat',jj)));
    Id_fc(jj,:) = Idscore_EC_FC_recon(2:200,:)';
    Id(jj,:) = Idscore_EC_recon(2:200,:)';
end

subplot(5,2,9); imagesc(Id_fc); axis square; caxis([0, 30]); colormap(gca,flipud(make_color_map(orig_col(8,:),32))); %colorbar;
box on;
handle = gca;
handle.LineWidth = 1.2;
%title('Communicability', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
xlabel('PCA # Components');
cb=colorbar;
cb.Position = handle.Position - [-0.2280 0 0.3197 0];

for jj = 1:numTasks
    [M,I] = max(Id_fc(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

subplot(5,2,10); imagesc(Id); axis square; caxis([0, 30]); colormap(gca,flipud(make_color_map(orig_col(8,:),32))); 
box on;
%title('Communicability', 'FontSize', 15);
yticks([1,2,3,4,5,6,7,8,9]); yticklabels(tasks_list);
xlabel('PCA # Components');
handle = gca;
handle.LineWidth = 1.2;
handle.YAxisLocation = 'right';
handle.Position = handle.Position + [-0.3 0 0 0];

for jj = 1:numTasks
    [M,I] = max(Id(jj,:));
    text(I,jj,sprintf('%d',round(M)), 'FontWeight', 'bold')
end

%% Save Optimal Reconstruction for everything

FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_%d.mat',(2*jj)-1)));
    [M, I] = max(Idiff_recon);
    FC = [FC; M I];
end

SP = [];
SP_FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SP_%d.mat',jj)));
    [M, I] = max(Idscore_SP_FC_recon);
    SP_FC = [SP_FC; M I];
    [M, I] = max(Idscore_SP_recon);
    SP = [SP; M I];
end

W = [];
W_FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_W_%d.mat',jj)));
    [M, I] = max(Idscore_W_FC_recon);
    W_FC = [W_FC; M I];
    [M, I] = max(Idscore_W_recon);
    W = [W; M I];    
end

EC = [];
EC_FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_EC_%d.mat',jj)));
    [M, I] = max(Idscore_EC_FC_recon);
    EC_FC = [EC_FC; M I];
    [M, I] = max(Idscore_EC_recon);
    EC = [EC; M I];
end

SI = [];
SI_FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SI_%d.mat',jj)));
    [M, I] = max(Idscore_SI_FC_recon);
    SI_FC = [SI_FC; M I];
    [M, I] = max(Idscore_SI_recon);
    SI = [SI; M I];
end

SS = [];
SS_FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_SS_%d.mat',jj)));
    [M, I] = max(Idscore_SS_FC_recon);
    SS_FC = [SS_FC; M I];
    [M, I] = max(Idscore_SS_recon);
    SS = [SS; M I];
end

TS = [];
TS_FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_TS_%d.mat',jj)));
    [M, I] = max(Idscore_TS_FC_recon);
    TS_FC = [TS_FC; M I];
    [M, I] = max(Idscore_TS_recon);
    TS = [TS; M I];
end

MFPT = [];
MFPT_FC = [];
for jj = 1:(numTasks)
    load(fullfile(path2ident, sprintf('Ident_MFPT_%d.mat',jj)));
    [M, I] = max(Idscore_MFPT_FC_recon);
    MFPT_FC = [MFPT_FC; M I];
    [M, I] = max(Idscore_MFPT_recon);
    MFPT = [MFPT; M I];
end

optimal = table(FC, SP_FC, SP, SI_FC, SI, MFPT_FC, MFPT, W_FC, W, EC_FC, EC, SS_FC, SS, TS_FC, TS);

save(fullfile(path2ident, 'optimal_recon.mat'), 'optimal');

