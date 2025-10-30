function[STMAout] = applyFaceROIs(varargin)


p = inputParser;

% Optional params
p.addParameter('workDir','D:\Dropbox\FaceExpressionAnalysis\STMA2021\',@ischar);
p.addParameter('subject','Barney', @ischar);
p.addParameter('date','210704',@ischar);
p.addParameter('TDTsessID','110137', @ischar);
p.addParameter('STMA', @isdouble);
p.addParameter('ROIs',{'mouth','Lcheek','Rcheek','nose','eyes','Lear','Rear'}, @iscell);
p.addParameter('verbose',0, @isscalar);

% parse
p.parse(varargin{:})

close all; 

date = p.Results.date;
subject = p.Results.subject;
TDTsessID = p.Results.TDTsessID;
workDir = [p.Results.workDir '/' subject '_' date '/' subject '-' date '-' TDTsessID];
ROIs = p.Results.ROIs;
STMAin = p.Results.STMA;
verbose = p.Results.verbose;

% apply the predefined crop regions 
for ROI = 1:length(ROIs) 
        thisROI = load([workDir '/faceROI_' ROIs{ROI} '.mat']);
        if verbose
            disp(['ROI= ' ROIs{ROI} 'size is ' num2str(thisROI.position2save)])
        end
        
        for ii= 1:size(STMAin,3)   
            STMAout.(ROIs{ROI})(:,:,ii) = imcrop(STMAin(:,:,ii),thisROI.position2save);
        end

end

end