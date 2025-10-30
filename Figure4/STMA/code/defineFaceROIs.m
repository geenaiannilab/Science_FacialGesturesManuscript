
function defineFaceROIs(varargin)

p = inputParser;

% Optional params
p.addParameter('workDir','D:\Dropbox\FaceExpressionAnalysis\STMA2021\',@ischar);
p.addParameter('subject','Barney', @ischar);
p.addParameter('date','210704',@ischar);
p.addParameter('TDTsessID','110137', @ischar);
p.addParameter('ROIs',{'mouth','Lcheek','Rcheek','nose','eyes','Lear','Rear'}, @iscell);

% parse
p.parse(varargin{:})

close all; 

date = p.Results.date;
subject = p.Results.subject;
TDTsessID = p.Results.TDTsessID;
workDir = [p.Results.workDir subject '_' date '/' subject '-' date '-' TDTsessID];
ROIs = p.Results.ROIs;
image2load = dir([workDir '/cropped_*']); image2load = ([image2load.folder '/' image2load.name]);

for ROI = 1:length(ROIs) 
    
    % for the first ROI, define mouth & save (option to keep aspect ratio &
    %               size for later ROIs) 
    if ROI == 1
        image = matfile(image2load);
        image.Properties.Writable = true; 
        
        if strcmp(subject,'Thor')
            firstFrame = rot90(image.catFrames(:,:,1));
        else
            firstFrame = image.catFrames(:,:,1);
        end
        
        figure; imagesc(firstFrame(:,:,1))  

        % Create an ROI in the current axis 
        roi = images.roi.Rectangle(gca);

        % Interactively draw a rectangular ROI 
        rect = drawrectangle('LineWidth',3);

        % Get position of the rectangle 
        position2save = floor(rect.Position);
        
        save([workDir '/faceROI_' ROIs{ROI} '.mat'],'position2save');
        save([workDir '/faceROI_template.mat'], 'rect');
    
    % for all other ROIs (besides the 1st), option to resize or keep the
    %           same size ROI, just translate it on image 
    else
        image = matfile(image2load);
        image.Properties.Writable = true; 
        
        if strcmp(subject,'Thor')
            firstFrame = rot90(image.catFrames(:,:,1));
        else
            firstFrame = image.catFrames(:,:,1);
        end
        
        figure; imagesc(firstFrame(:,:,1)); title(['define: ' ROIs{ROI}]);  

        ROItmp = load([workDir '/faceROI_template.mat']);

        prompt = (['ROI is ' ROIs{ROI} ' do you want to resize?  [y/n]: ']);
        answer = input(prompt,'s');

        if strcmp(answer,'n')
            rect = drawrectangle('LineWidth',3,'AspectRatio',ROItmp.rect.AspectRatio,...
            'Position', ROItmp.rect.Position,'InteractionsAllowed','translate');

            position2save = customWait(rect);

            save([workDir '/faceROI_' ROIs{ROI} '.mat'],'position2save');

        elseif strcmp(answer,'y')
            rect = drawrectangle('LineWidth',3,'AspectRatio',ROItmp.rect.AspectRatio,...
            'Position', ROItmp.rect.Position,'InteractionsAllowed','all');

            position2save = customWait(rect);

            save([workDir '/faceROI_' ROIs{ROI} '.mat'],'position2save');
        else
            disp('pls respond y/n')
            return
        end

        close all; 

    end
end

