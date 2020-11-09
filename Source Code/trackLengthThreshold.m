function [keptTracks,noTracksFound] = trackLengthThreshold(track,threshold,specificMito)
trackNum = 0;
noTracksFound = 0;
for i = 1:length(track)
    if nargin>2
        if ~ismember(i,specificMito)
            continue
        end
    end
    if length(track(i).frame)>threshold
        trackNum = trackNum + 1;
        tempTrack = track(i);
        tempTrack.OGid = i;
        keptTracks(trackNum) = tempTrack;
    end
end
if ~trackNum
    keptTracks = 0;
    noTracksFound = 1;
    return
end
for i = 1:length(keptTracks)
    if keptTracks(i).fusion(end)
        if ~sum(keptTracks(i).fusion(end)==[keptTracks.OGid])
            keptTracks(i).fusion(end) = 0;
        else
            keptTracks(i).fusion(end) = find(keptTracks(i).fusion(end)==[keptTracks.OGid]);
        end
    end
    if keptTracks(i).fission(1)
        if ~sum(keptTracks(i).fission(1)==[keptTracks.OGid])
            keptTracks(i).fission(1) = 0;
        else
            keptTracks(i).fission(1) = find(keptTracks(i).fission(1)==[keptTracks.OGid]);
        end
    end

end

end