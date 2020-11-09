function argStats = getDynamicStats(arg)

numTracks = length(arg);
mean = nan(1,numTracks);
std = nan(1,numTracks);
CoV = nan(1,numTracks);

for trackNum = 1:numTracks
    mean(trackNum) = nanmean(arg{trackNum});
    std(trackNum) = nanstd(arg{trackNum});
    CoV(trackNum) = std(trackNum)/mean(trackNum);
end

argStats.mean = mean;
argStats.std = std;
argStats.CoV = CoV;

end