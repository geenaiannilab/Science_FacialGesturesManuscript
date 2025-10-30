function obhvOut = spacedBhvsPSID(obhv,minRest,varargin)

%minRest = 1; % in seconds, min desired rest period between events

if length(varargin) == 1
    overwriteFlag = varargin{1};
    verbose = 0;
elseif length(varargin) == 2
    overwriteFlag = varargin{1};
    verbose = varargin{2};
else
    overwriteFlag = 0;
    verbose = 0;
end

if contains(obhv.Info.blockname,'Thor')
    srBhv = 1 / median(diff(obhv.FT2dat.time)); % Hz camera
elseif contains(obhv.Info.blockname,'Barney') || contains(obhv.Info.blockname,'Dexter') 
    srBhv = 1 / median(diff(obhv.FT1dat.time)); % Hz camera
end
minFrames = ceil(minRest * srBhv);

ev = 1;
evOut = 1;
removedEvs = 0;
evSpaced = obhv.evScore(1,:);

while ev < size(obhv.evScore,1)
    thisEv = cell2mat(obhv.evScore(ev,1:2));
    nextEv = cell2mat(obhv.evScore(ev+1,1:2));
    if ~(nextEv(1) - thisEv(2) >= minFrames)
        ev = ev + 1;
        removedEvs = removedEvs +1;
        continue
    else
        evSpaced(evOut+1,:) = obhv.evScore(ev+1,:);
        evOut = evOut + 1;
    end
    
    ev = ev + 1;

end

if verbose
    disp([num2str(size(obhv.evScore,1) - removedEvs) ' of ' num2str(size(obhv.evScore,1)) '  events remain'])
end
if overwriteFlag
     obhv.evScore = evSpaced;
end

obhvOut = obhv;

end