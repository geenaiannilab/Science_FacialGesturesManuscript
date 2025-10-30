function examineDecodedBhvPCs(bhvTest, bhvTestPred1, times2plot, accuracyVal, st, vt, pixHeight,pixWid,evOnsets)


% eigenvalues 
    latent = diag(st);
    figure; subplot(2,1,1);bar(latent); title('eigenvalues');
    subplot(2,1,2); plot(cumsum(latent/sum(latent)), 'o') 
    title('elbow plot'); xlabel('PC#'); ylabel('cumulative variance')

    % visualize each spatial map, & time component, colored by what bhv taking
    %   place during that time period
    mean_vt = mean(vt(:));
    std_vt = std(vt(:));
    vt_min = mean_vt - 3*std_vt;
    vt_max = mean_vt + 3*std_vt;
    clims = [vt_min vt_max];

    for facePC = 1:size(bhvTest,2)
        figure; ax1 = subplot(211);

        vt_PC = vt(:,facePC); % PC direction/axis, or spatial map 
        img = reshape(vt_PC,pixHeight,pixWid);
        imagesc(img, clims);
        colorbar
        colormap(ax1,redblue)
        title(['face PC' num2str(facePC) ' CC= ' num2str(accuracyVal(facePC))]); 

        ax2 = subplot(212); % score, or temporal component scaled by eigenvalue
        plot(times2plot, bhvTest(:,facePC),'Color',[0.8500 0.3250 0.0980]); hold on; plot(times2plot, bhvTestPred1(:,facePC),'Color',[0 0.4470 0.7410]);
        plot([evOnsets; evOnsets], repmat(ylim',1,size(evOnsets,2)), '--k')
        xlim([times2plot(1) times2plot(end)])
        xlabel('time, s')
        hold off
        legend('true bhv','decoded')

    end
end
