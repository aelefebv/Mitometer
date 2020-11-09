function [mitoDynamics,CoMDynamics] = getDynamics(micronPerPixel,secondsPerFrame,track,mitoCoM)
%track should be a structure of mitochondrial tracks, and mitoCoM should be the
%CoM stats for which you want to get the dynamics. MitoCoM can be left
%blank.

numTracks = length(track);
numFrames = max([track.frame]);

distance = cell(1,numTracks);
displacement = cell(1,numTracks);

speed = cell(1,numTracks);
velocity = cell(1,numTracks);

travelAngle = cell(1,numTracks);

MSD = nan(length(track),numFrames);
MSDstd = nan(length(track),numFrames);

    
%Track Dynamics
for trackNum = 1:numTracks
    
    [MSD(trackNum,1:(track(trackNum).frame(end)-track(trackNum).frame(1))),MSDstd(trackNum,1:(track(trackNum).frame(end)-track(trackNum).frame(1)))] = getMSD(track(trackNum),micronPerPixel);
    
    speed{trackNum}(1) = nan;
    velocity{trackNum}(1) = nan;
    distance{trackNum}(1) = 0;
    displacement{trackNum}(1) = 0;
    
    travelAngle{trackNum}(1) = nan;
    
    coordsO = [track(trackNum).WeightedCentroid(1),track(trackNum).WeightedCentroid(2)];
        
%(-atan2d((trackToAdd.Col(end)-oldTrack.Col(end)),(trackToAdd.Row(end)-oldTrack.Row(end))))
%migrationAngle(frameNum) = -atan2((coordscurr(1)-coordslast(1)),(coordscurr(2)-coordslast(2)));
%track(i).angle = angle0to180(90-track(i).angle);
       
    for idx = 2:length(track(trackNum).frame) 
        coordsF = [track(trackNum).WeightedCentroid(idx*2-1),track(trackNum).WeightedCentroid(idx*2)];
        coordsL = [track(trackNum).WeightedCentroid((idx-1)*2-1),track(trackNum).WeightedCentroid((idx-1)*2)];

        distance{trackNum}(idx) = norm(coordsF-coordsL)*micronPerPixel;
        displacement{trackNum}(idx) = norm(coordsF-coordsO)*micronPerPixel;
       
        speed{trackNum}(idx) = abs(distance{trackNum}(idx))/(track(trackNum).frame(idx)-track(trackNum).frame(idx-1))/secondsPerFrame;
        velocity{trackNum}(idx) = (displacement{trackNum}(idx))/(track(trackNum).frame(idx)-track(trackNum).frame(idx-1))/secondsPerFrame;
        
        travelAngle{trackNum}(idx) = angle0to180(90-(-atan2d((coordsF(1)-coordsL(1)),(coordsF(2)-coordsL(2)))));
    end
end

%Center of Mass Dynamics
if nargin == 4
    
    CoMdistance = zeros(1,numFrames);
    CoMdisplacement = zeros(1,numFrames);

    CoMspeed = zeros(1,numFrames);
    CoMvelocity = zeros(1,numFrames);

    CoMcoordsO = mitoCoM(1).centroid;
    
    CoMtravelAngle = zeros(1,numFrames);
    
    CoMspeed(1) = nan;
    CoMdistance(1) = 0;
    CoMdisplacement(1) = 0;
    
    CoMvelocity(1) = nan;

    CoMtravelAngle(1) = nan;
    
    for idx = 2:length(mitoCoM)
        
        CoMcoordsF = mitoCoM(idx).centroid;
        CoMcoordsL = mitoCoM(idx-1).centroid;
%         coordsF = [track(trackNum).WeightedCentroid(idx*2-1),track(trackNum).WeightedCentroid(idx*2)];

        CoMdistance(idx) = norm(CoMcoordsF-CoMcoordsL)*micronPerPixel;
        CoMdisplacement(idx) = norm(CoMcoordsF-CoMcoordsO)*micronPerPixel;

        CoMspeed(idx) = abs(CoMdistance(idx))/secondsPerFrame;
        CoMvelocity(idx) = (CoMdisplacement(idx))/secondsPerFrame;

        CoMtravelAngle(idx) = angle0to180(90-(-atan2d((CoMcoordsF(1)-CoMcoordsL(1)),(CoMcoordsF(2)-CoMcoordsL(2)))));
    end
end

mitoDynamics.distance = distance;
mitoDynamics.displacement = displacement;
mitoDynamics.speed = speed;
mitoDynamics.velocity = velocity;
mitoDynamics.MSD = MSD;
mitoDynamics.MSDstd = MSDstd;
mitoDynamics.travelAngle = travelAngle;

if nargin == 4
    CoMDynamics.distance = CoMdistance;
    CoMDynamics.displacement = CoMdisplacement;
    CoMDynamics.speed = CoMspeed;
    CoMDynamics.velocity = CoMvelocity;
    CoMDynamics.travelAngle = CoMtravelAngle;
else
    CoMDynamics = struct();
end

end