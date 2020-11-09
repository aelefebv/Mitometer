function deltaMatrix = getDeltaMatrix(mito,track,field)

numMito = length(mito);
numTracks = length(track);

%initialize
deltaMatrix = zeros(numMito,numTracks);

deltaMito = [mito.(field)];

%we simply take the difference in the feature’s value between the current
%time frame’s mitochondria to those of existing tracks (e.g. the difference
%in surface area) divided by the existing track’s value, producing a
%normalized difference.
for trackNum = 1:numTracks
    deltaMatrix(:,trackNum) = (deltaMito'-track(trackNum).(field)(end))./(track(trackNum).(field)(end));
end

end