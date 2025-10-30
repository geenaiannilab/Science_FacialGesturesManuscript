function[regions] = getChannel2CorticalRegionMapping(subject, collapseFlag)

%%%%%% this is necessary since channel #s were different in each animal
%%%%%% (Barney = 8 arrays, Thor = 6 arrays) 

if strcmpi(subject, 'Barney') || strcmpi(subject, 'Dexter')
    regions{1}.label = 'S1';
    regions{1}.channels = 1:32;
    
    if collapseFlag % collapse M1m, M1lat, and F4; collapse PMvm & Pmvlat, collapse M3 & ACC
        regions{2}.label = 'M1';
        regions{2}.channels = 33:128;
        regions{3}.label = 'PMv';
        regions{3}.channels = 129:192;
        regions{4}.label = 'M3';
        regions{4}.channels = 193:240;
        regions{5}.label = 'All';
        regions{5}.channels = 1:240;
    else
        regions{2}.label = 'M1m';
        regions{2}.channels = 33:64;
        regions{3}.label = 'M1lat';
        regions{3}.channels = 65:96;
        regions{4}.label = 'F4';
        regions{4}.channels = 97:128;
        regions{5}.label = 'PMvm';
        regions{5}.channels = 129:160;
        regions{6}.label = 'PMvlat';
        regions{6}.channels = 161:192;
        regions{7}.label = 'M3';
        regions{7}.channels = 193:224;
        regions{8}.label = 'ACC';
        regions{8}.channels = 225:240;
        regions{9}.label = 'All';
        regions{9}.channels = 1:240;
    end
    
elseif strcmpi(subject,'Thor')
    regions{1}.label = 'S1';
    regions{1}.channels = 1:32;
    
    if collapseFlag % collapse M1m, M1lat, collapse M3 & Acc
        regions{2}.label = 'M1';
        regions{2}.channels = 33:96;
        regions{3}.label = 'PMv';
        regions{3}.channels = 97:128;
        regions{4}.label = 'M3';
        regions{4}.channels = 129:192;
        regions{5}.label = 'All';
        regions{5}.channels = 1:192;
    else
        regions{2}.label = 'M1m';
        regions{2}.channels = 33:64;
        regions{3}.label = 'M1lat';
        regions{3}.channels = 65:96;
        regions{4}.label = 'PMv';
        regions{4}.channels = 97:128;
        regions{5}.label = 'M3';
        regions{5}.channels = 129:160;
        regions{6}.label = 'ACC';
        regions{6}.channels = 161:192;
        regions{7}.label = 'All';
        regions{7}.channels = 1:192;
    end
else
    disp('unrecognized subject')
end

end
