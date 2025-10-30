
%%% after running PSID_trials_DLC, run this to compare results of decoding models across regions
%%%
close all; clear all; 
subject = 'Barney';
date = '210704';
workDir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/PSID/']);

trainData = 'allTrialsPCs';
bhvs2plot = [1 2 4];
regions = {'M1','PMv','M3','PMV+M3','M1+PMv','M1+M3','M1+PMv+M3'};

%
if strcmpi(trainData,'LSTrialPCs')
    for rr = 1:length(regions)  
        disp(regions{rr})
        data = load([workDir 'results_' regions{rr} '_' trainData '.mat']);
        
        bhvTestTrialAvgFoldAvg = mean(data.bhvTestTrialAvg,4);    % nTimepts x nPredictors x nBhvs
        bhvPredTrialAvgFoldAvg = mean(data.bhvPredTrialAvg,4);   % nTimepts x nPredictors x nBhvs
        latentPredTrialAvgFoldAvg = mean(data.latentPredTrialAvg,4); % nTimepts x nLatent x nBhvs
        
        CCtrialAvgFoldAvg = mean(data.CCtrialbased,3); % average CC per predictor (x-pos AND y-pos)
        CCtrialAvgSEM = (std(data.CCtrialbased,[],3)' ./ (sqrt(data.nFolds)) ) ;
        CCtrialbasedMarkerAvg = mean(reshape(CCtrialAvgFoldAvg, 2, [])); % avg CC per marker (over x&y)
        
        figure; 
        c = distinguishable_colors(data.nx);
        for i = 1:data.nx
            subplot(2,2,i)
            plot(data.taxis,latentPredTrialAvgFoldAvg(:,i),'Color', c(i,:),'LineWidth',4)
            xline(0,'--r')
            xlabel('time relative to moveOnset, s'); ylabel('au')
            title(['Latent Dimension #' num2str(i)])
        end
        sgtitle(['Behaviorally Relevent Dimensions in  ' regions{rr} ' Lipsmack'])
        
        
        figure;
        for i = 1:data.nx
            thisDim = squeeze(data.latentPredTrialAvg(:,i,1,:))';
            shadedErrorBar(data.taxis, thisDim,...
                {@(m) (mean(thisDim)),...
        @(x) (std(thisDim)) ./ sqrt(size(thisDim,2))},'lineprops',{'Color',c(i,:)});
        
        end
        xline(0,'--k')
        xlabel('time relative to moveOnset, s'); ylabel('au');
        hold off;
        legend('BhvRelv Dim1', 'BhvRelv Dim2','BhvRelv Dim3','Irrelv');
        sgtitle(['Behaviorally Relevent Dimensions in  ' regions{rr} ' :Lipsmack'])
        
        figure; 
             x = smoothdata(latentPredTrialAvgFoldAvg(:,1),'gaussian',10)';
            y = smoothdata(latentPredTrialAvgFoldAvg(:,2),'gaussian',10)';
            z = smoothdata(latentPredTrialAvgFoldAvg(:,3),'gaussian',10)';
            c = data.taxis;
            scatter3(x,y,z,40,c,'linew',8)
            view(3); colorbar
            title(['Neural Traj in  ' regions{rr} ' bhv = Lipsmack' ]);
            

        %  plot CC (trial & fold-avg) per marker (x and Y position) 
        tmp(:,rr) = CCtrialAvgFoldAvg;
        tmpSEM(:,rr) = CCtrialAvgSEM;

        figure; bar(CCtrialAvgFoldAvg); 
        axes = gca;
        axes.XTick = data.bhvPCs;
        axes.XTickLabel = data.bhvPCs;
        xlabel('faceMarker'); ylabel('CC'); ylim([-0.2 0.8])
        title(['CrossValidated, CC of Predicted vs. Actual Face Markers; ' regions{rr}])
        %subtitle(['Grand Mean = ' mean(data.CCMean)])
    
    end
    
    
    compareDecoding = tmp; 
    compareDecodingErrorBars = tmpSEM;

    figure; 
    b = bar(compareDecoding, 'grouped');
    hold on;
    
    % adding error bars 
    [ngroups,nbars] = size(compareDecoding);
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',compareDecoding,compareDecodingErrorBars,'k','linestyle','none', 'linew',1);
    hold off

    axes = gca; 
    axes.XTick = data.bhvPCs; axes.XTickLabel = data.bhvPCs;
    xlabel('Face PC','FontSize',16); ylabel('Decoding, CC','FontSize',16); ylim([-0.1 1]); legend(regions)
    title('Decoding Performance Across All Models','FontSize',24)

else

    for rr = 1:length(regions)  
        disp(regions{rr})
        data = load([workDir '/results_' regions{rr} '_' trainData '.mat']);
        
        for bhv = 1:length(bhvs2plot)
            bhvTestTrialAvgFoldAvg = mean(data.bhvTestTrialAvg(:,:,bhv,:),4);    % nTimepts x nPredictors x nBhvs
            bhvPredTrialAvgFoldAvg = mean(data.bhvPredTrialAvg(:,:,bhv,:),4);   % nTimepts x nPredictors x nBhvs
            latentPredTrialAvgFoldAvg = mean(data.latentPredTrialAvg(:,:,bhv),4); % nTimepts x nLatent x nBhvs
            
            % CCtrialbased is nPredictors x nBhvs x nFolds 
            CCtrialAvgFoldAvg = mean(data.CCtrialbased(:,bhv,:),3); % average CC per predictor (x-pos AND y-pos)
            CCtrialAvgSEM = (std(data.CCtrialbased(:,bhv,:),[],3)' ./ (sqrt(data.nFolds)) ) ;
            CCtrialbasedMarkerAvg = mean(reshape(CCtrialAvgFoldAvg, 2, [])); % avg CC per marker (over x&y)
        
            figure; 
            c = distinguishable_colors(data.nx);
            for i = 1:data.nx
                subplot(2,2,i)
                plot(data.taxis,latentPredTrialAvgFoldAvg(:,i),'Color', c(i,:),'LineWidth',4)
                xline(0,'--r')
                xlabel('time relative to moveOnset, s'); ylabel('au')
                title(['Latent Dimension #' num2str(i)])
            end
            sgtitle(['Behaviorally Relevent Dimensions in  ' regions{rr} ' bhv = ' num2str(bhvs2plot(bhv))]);
            
            
            figure;
            for i = 1:data.nx
                thisDim = squeeze(data.latentPredTrialAvg(:,i,bhv,:))';
                shadedErrorBar(data.taxis, thisDim,...
                    {@(m) (mean(thisDim)),...
            @(x) (std(thisDim)) ./ sqrt(size(thisDim,2))},'lineprops',{'Color',c(i,:)});
            
            end
            xline(0,'--k')
            xlabel('time relative to moveOnset, s'); ylabel('au');
            hold off;
            legend('BhvRelv Dim1', 'BhvRelv Dim2','BhvRelv Dim3','Irrelv');
            sgtitle(['Behaviorally Relevent Dimensions in  ' regions{rr} ' bhv = ' num2str(bhvs2plot(bhv))]);
         

         %  plot CC (trial & fold-avg) per marker (x and Y position) 
            tmp(:,bhv, rr) = CCtrialAvgFoldAvg;
            tmpSEM(:,bhv,rr) = CCtrialAvgSEM;
    
        
            
        end
    end

    for bhv = 1:length(bhvs2plot)
        compareDecoding =   squeeze(tmp(:,bhv,:)); 
        compareDecodingErrorBars = squeeze(tmpSEM(:,bhv,:));

        figure; 
        b = bar(compareDecoding, 'grouped');
        hold on;

         % adding error bars 
        [ngroups,nbars] = size(compareDecoding);
        x = nan(nbars, ngroups);
        for i = 1:nbars
            x(i,:) = b(i).XEndPoints;
        end
        % Plot the errorbars
        errorbar(x',compareDecoding,compareDecodingErrorBars,'k','linestyle','none', 'linew',1);
        hold off
    
        axes = gca; 
        axes.XTick = data.bhvPCs; axes.XTickLabel = data.bhvPCs;
        xlabel('Face PC','FontSize',16); ylabel('Decoding, CC','FontSize',16); ylim([-0.1 1]); legend(regions)
        title('Decoding Performance Across All Models','FontSize',24)

                
    end
        
end