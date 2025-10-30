%%%%%%%  Creates data in Figure 2B 
%%%%%%%%  Single-Cell Activity and Selectivity in Cortical Face-Motor Regions
%%%%%%%%
%%%%%%%% Sample population activity of all simultaneously recorded cells 
%%%%%%%% (top, sorted by area) and facial marker positions (bottom) 
%%%%%%%%% GRI 09/11/2025 

close all; clear all;
set(0,'defaultAxesFontSize',24);
set(0,'defaultAxesFontWeight','bold');


date = '210704';
subject ='Barney';
chls = 1:240;
win = 0.05; %bin size, sec
markers2use =  'all'; excludeTongue = 1; 

bhvColorMap = [1 1 1 ; 1 0 0 ;0 0 1; 1 1 1; 0 1 0];
markersColorMap = 'hsv';

% camera variables different across subjects
if strcmpi(subject, 'Barney') || strcmpi(subject, 'Dexter')
    faceCameraStream = 'FT1dat';
    bhvSR = 70;

elseif strcmpi(subject, 'Thor')
    faceCameraStream = 'FT2dat';
    bhvSR = 120;
end

workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis/']);

flist = natsortfiles(dir([workdir '/bhvspikes_sub*.mat']));
flist = fullfile({flist.folder}, {flist.name})'; 

for session = 6%:length(flist)
    load([flist{session}]);
    
    if isempty(obhv.evScore)
        disp('no scored bhv, skipping')
        continue
    end
    
    tmin = obhv.(faceCameraStream).time(1);
    pulseLength = length(obhv.(faceCameraStream).time);
    nFramesExported = size(obhv.(faceCameraStream).Markers,1);

    tmax = obhv.(faceCameraStream).time(nFramesExported);

    newTimeAxis = tmin:win:tmax;

    rowIndx = 1;
    raster = [];
    spUnitID = [];
    bhv2plot = [];
    
    for ch = chls
        
        units = 1:size(obhv.spikes.el(ch).sp,2);
        
        for uindx = units
            if isnan(obhv.spikes.el(ch).sp(uindx).times)
                continue 
            else
                
            spikes = obhv.spikes.el(ch).sp(uindx).times;
            spUnitID(rowIndx,:) = [ch; uindx];
            raster(rowIndx,:) = histcounts(spikes,newTimeAxis);
            rowIndx = rowIndx +1;
            
            end
        end
        
    end
    
    raster = raster ./win;

    bhv2plot(:,1) = cell2mat(obhv.evScore(:,1)); % onset frame
    bhv2plot(:,2) = cell2mat(obhv.evScore(:,2)); % offset frame
    bhv2plot(:,3) = cell2mat(obhv.evScore(:,5)); % bhv
    
    % remove any events after DLC cut off 
    bhv2plot = bhv2plot(~any(bhv2plot >= nFramesExported,2),:); 
    allBhv = binBhv( bhv2plot, obhv.(faceCameraStream).time, newTimeAxis);
    
    %%%%%%% extract DLC marker positions & labels from obhv
    %%%%%
    % 
    DLCin = obhv.(faceCameraStream);
    DLCtimeBase = DLCin.time(1:nFramesExported);
    
    [DLCmarkers, markerLabels] = processDLCMarkers('DLCin', DLCin, 'markers2use', markers2use,'excludeTongue', true,'plotFlag',false);
    DLCOut = rebinBhv(DLCmarkers, DLCtimeBase, win);    
  
    %%%%%
    %%%%% PLOTTING
    fig1 = figure('units','normalized','outerposition',[0 0 1 1]); 
    fig1.Position =[ 0.3690    0.0855    0.6138    0.7953];

    h1 = axes(fig1,'position',[0.13 0.35 0.75, 0.6]);
    h2 = axes(fig1,'position',[0.13 0.16 0.75, 0.18]);
    h3 = axes(fig1,'position',[0.13 0.08 0.75, 0.07]);
    
    % raster 
    axes(h1);
    imagesc(newTimeAxis, 1:size(raster,1), raster); colormap(h1,(1-gray)); cb = colorbar; 
    ylab = ylabel('Cell ID');

    h1.YLabel
    h1.XTick = [];
    h1.YTick = [];
    h1.YTickLabel = [];
    cb.Position = [0.89 0.35 0.0106    0.6001];
    ylabel(cb,'Spikes/second','FontSize',18,'Rotation',270);
    cb.Label.Units = 'normalized';cb.Label.Position = [6.7 0.5 0 ];

    caxis([0 100])
   
    % face marker traces
    axes(h2)
    plot(newTimeAxis,DLCOut(:,50:76),'linew',2); colormap(h2, hsv);
    ylab1 = ylabel('Marker Position','FontSize',16,'Rotation',45);
    h2.YAxisLocation = 'left';
    h2.YTick = [];  h2.XTick = [];
    

    % bhv scoring 
    axes(h3) 
    imagesc(newTimeAxis, 1:4,allBhv(2:5,:)); colormap(1-gray);%colormap(h3,bhvColorMap);
    h3.YTick = '';
    xlabel('Time, s')
    ylab2 = ylabel('Gesture','FontSize',18,'Rotation',45);

    linkaxes([h1 h2 h3],'x')
    
    axes(h1);title([extractAfter(obhv.Info.blockname,'Barney-')])
    ylab1.Position = [328 400.0002   -1.0000]; 
    axes(h2); ylab2.Position= [328    2.5000    1.0000];
    xlim([330 360]);
   
    % without the bhv scoring 
    fig2 = figure('units','normalized','outerposition',[0 0 1 1]); 
    fig2.Position =[ 0.3690    0.0855    0.6138    0.7953];
 
    h4 = axes(fig2,'position',[0.13 0.3 0.75, 0.6]);
    h5 = axes(fig2,'position',[0.13 0.1 0.75, 0.19]);
    
    % raster 
    axes(h4);
    imagesc(newTimeAxis, 1:size(raster,1), raster); colormap(h4,(1-gray)); cb = colorbar; ylabel(cb,'Spikes/second','FontSize',18,'Rotation',270)
    ylabel('Cell ID')
    h4.YLabel
    h4.XTick = [];
    h4.YTick = [];
    h4.YTickLabel = [];
    cb.Position = [0.89 0.3 0.0106    0.6001];
    cb.Label.Units = 'normalized';cb.Label.Position = [6.7 0.5 0 ];

    caxis([0 100])
   
    % face marker traces
    axes(h5)
    plot(newTimeAxis,DLCOut(:,50:76),'linew',2); colormap(h2, hsv);
    ylabel('Marker Position','FontSize',18)
    h5.YAxisLocation = 'left';
    h5.YTick = [];
   
    linkaxes([h4 h5],'x')
    xlim([330 360])
    axes(h4);title([extractAfter(obhv.Info.blockname,'Barney-')])
    %subtitle([obhv.Info.sessType])
end
fig1; xlim([330 360]);


function allBhv = binBhv( scoring, FT1ticks,timeAxis)
    
allBhv = zeros(7, length(timeAxis)-1);
win = timeAxis(2) - timeAxis(1);

for i=1:size(scoring,1)
   
    indx1 = round(FT1ticks(scoring(i,1))/win);
    indx2 = round(FT1ticks(scoring(i,2))/win);
    stype = scoring(i,3);
    if ~isnan(stype)
        allBhv(stype,indx1:indx2) = stype;
    end

end
end

function bhvOut = rebinBhv(bhvIn, timeBaseOld, newBinsize)
    
downsampledTimeBaseVector = timeBaseOld(1):newBinsize:timeBaseOld(end);
bhvOut = zeros(length(downsampledTimeBaseVector), size(bhvIn,2));

for j = 1:size(bhvIn,2)
    bhvOut(:,j) = interp1(timeBaseOld, bhvIn(:,j), downsampledTimeBaseVector);
end

end

