function [DLCout, markerLabels] = processDLCMarkers(varargin)

p = inputParser;

%p.addParameter('cameraSync',cameraSync, @isdouble);
p.addParameter('markers2use','all');
p.addParameter('DLCin','DLCin', @isstruct);

p.addParameter('plotFlag',false,@islogical);
p.addParameter('excludeTongue',true, @islogical);
p.addParameter('velocityFlag', false, @islogical);
p.addParameter('downsampleFlag',false,@islogical)
p.addParameter('downsampleTarget',0.05)
p.addParameter('smoothFlag',true);
p.addParameter('smoothSize',20);

% parse
p.parse(varargin{:})


DLCin = p.Results.DLCin;
markers2use = p.Results.markers2use;
plotFlag = p.Results.plotFlag;
velocityFlag = p.Results.velocityFlag;
excludeTongue = p.Results.excludeTongue;
smoothFlag = p.Results.smoothFlag;
smoothSize = p.Results.smoothSize;
downsampleFlag = p.Results.downsampleFlag;
downsampleTarget = p.Results.downsampleTarget;
markerTrajs = DLCin.Markers;
markerLabels = DLCin.MarkerLabels;

bhvMeasure = table2array(markerTrajs);

if strcmpi(markers2use,'all')

    colToExclude = 3:3:size(bhvMeasure,2);
    datCols2keep = setdiff(1:size(bhvMeasure,2),colToExclude);  
    
    markerCols2keep = 1:3:length(markerLabels);
    DLCout = bhvMeasure(:,datCols2keep);
    markerLabels = markerLabels(markerCols2keep);
    
    confidence = bhvMeasure(:,3:3:end);
    prctReplaced = sum(confidence <= 0.5,'all') ./ numel(confidence);
    disp([num2str(prctReplaced*100) ' % datapts are < 0.5 confidence; replacing them!']);

    for marker = 1:length(markerCols2keep)
        lowC = confidence(:,marker) <= 0.5;
        if ~isempty(find(lowC))
            DLCout(lowC,(marker*2)-1) = mean(DLCout(:,(marker*2)-1)); % xpos
            DLCout(lowC,marker*2) = mean(DLCout(:,marker*2)); % y pos
        end
    end


    if excludeTongue
        DLCout = DLCout(:,1:end-2); % last two columns are tongue 
        markerLabels = markerLabels(1:end-1);

    end

else
    idx = find(ismember(markerLabels,markers2use));
    idxy = sort(cat(2, idx, idx +1))';
    DLCout = bhvMeasure(:,idxy);
    markerLabels = markers2use;
    
    if excludeTongue
        DLCout = DLCout(:,1:end-2); % last two columns are tongue 
        markerLabels = markerLabels(1:end-1);

    end
end

if smoothFlag
    DLCout = smoothdata(DLCout,'gaussian', smoothSize);

end

if downsampleFlag
    
    timeBaseOld = DLCin.time(1:size(DLCout,1));
    downsampledTimeBaseVector = 0:downsampleTarget:timeBaseOld(end);

    DLCoutDS = zeros(length(downsampledTimeBaseVector),size(DLCout,2));

    for j = 1:size(DLCout,2)
        DLCoutDS(:,j) = interp1(timeBaseOld, DLCout(:,j), downsampledTimeBaseVector);
    end

    DLCout = DLCoutDS;

end

if velocityFlag
    DLCout = gradient(DLCout);
end

if plotFlag
        figure; 
        for i = 1:length(markerLabels)
            xpos = (i*2) - 1;
            ypos = (i*2);
            x = DLCout(:,xpos);
            y = -DLCout(:,ypos);
            plot(x,y); hold on 
        
        end
        hold off
end

end
