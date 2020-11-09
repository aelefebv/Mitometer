function [mitoDynamics,CoMDynamics] = getDynamics3d(micronPerPixel,secondsPerFrame,track,zDistanceMicrons,mitoCoM)
%track should be a structure of mitochondrial tracks, and mitoCoM should be the
%CoM stats for which you want to get the dynamics. MitoCoM can be left
%blank.

zRatio = zDistanceMicrons/micronPerPixel;

numTracks = length(track);
numFrames = max([track.frame]);

distance = cell(1,numTracks);
displacement = cell(1,numTracks);

speed = cell(1,numTracks);
velocity = cell(1,numTracks);

xyAngle = cell(1,numTracks);
xzAngle = cell(1,numTracks);
yzAngle = cell(1,numTracks);


MSD = nan(length(track),numFrames);
    
%Track Dynamics
for trackNum = 1:numTracks
    
    MSD(trackNum,1:(track(trackNum).frame(end)-track(trackNum).frame(1))) = getMSD(track(trackNum),micronPerPixel);
    
    speed{trackNum}(1) = nan;
    velocity{trackNum}(1) = nan;
    distance{trackNum}(1) = 0;
    displacement{trackNum}(1) = 0;
    
    xyAngle{trackNum}(1) = nan;
    xzAngle{trackNum}(1) = nan;
    yzAngle{trackNum}(1) = nan;

    coordsO = [track(trackNum).WeightedCentroid(1),track(trackNum).WeightedCentroid(2),track(trackNum).WeightedCentroid(3)*zRatio];
        
%(-atan2d((trackToAdd.Col(end)-oldTrack.Col(end)),(trackToAdd.Row(end)-oldTrack.Row(end))))
%migrationAngle(frameNum) = -atan2((coordscurr(1)-coordslast(1)),(coordscurr(2)-coordslast(2)));
%track(i).angle = angle0to180(90-track(i).angle);
       
    for idx = 2:length(track(trackNum).frame) 
        coordsF = [track(trackNum).WeightedCentroid(idx*3-2),track(trackNum).WeightedCentroid(idx*3-1),track(trackNum).WeightedCentroid(idx*3)*zRatio];
        coordsL = [track(trackNum).WeightedCentroid((idx-1)*3-2),track(trackNum).WeightedCentroid((idx-1)*3-1),track(trackNum).WeightedCentroid((idx-1)*3)*zRatio];

        distance{trackNum}(idx) = norm(coordsF-coordsL)*micronPerPixel;
        displacement{trackNum}(idx) = norm(coordsF-coordsO)*micronPerPixel;
       
        speed{trackNum}(idx) = abs(distance{trackNum}(idx))/(track(trackNum).frame(idx)-track(trackNum).frame(idx-1))/secondsPerFrame;
        velocity{trackNum}(idx) = (displacement{trackNum}(idx))/(track(trackNum).frame(idx)-track(trackNum).frame(idx-1))/secondsPerFrame;
        
%         travelAngle{trackNum}(idx) = angle0to180(90-(-atan2d((coordsF(1)-coordsL(1)),(coordsF(2)-coordsL(2)))));
        xyAngle{trackNum}(idx) = angle0to180(90-(-atan2d((coordsF(1)-coordsL(1)),(coordsF(2)-coordsL(2)))));
        xzAngle{trackNum}(idx) = angle0to180(90-(-atan2d((coordsF(1)-coordsL(1)),(coordsF(3)-coordsL(3)))));
        yzAngle{trackNum}(idx) = angle0to180(90-(-atan2d((coordsF(2)-coordsL(2)),(coordsF(3)-coordsL(3)))));

    end
end

%Center of Mass Dynamics
if nargin > 4
    
    CoMdistance = zeros(1,numFrames);
    CoMdisplacement = zeros(1,numFrames);

    CoMspeed = zeros(1,numFrames);
    CoMvelocity = zeros(1,numFrames);

    CoMcoordsO = mitoCoM(1).centroid;
    
    CoMxyAngle = zeros(1,numFrames);
    CoMxzAngle = zeros(1,numFrames);
    CoMyzAngle = zeros(1,numFrames);
    
    CoMspeed(1) = nan;
    CoMdistance(1) = 0;
    CoMdisplacement(1) = 0;
    
    CoMvelocity(1) = nan;

    CoMxyAngle(1) = nan;
    CoMxzAngle(1) = nan;
    CoMyzAngle(1) = nan;

    for idx = 2:length(mitoCoM)
        
        CoMcoordsF = [mitoCoM(idx).centroid(1),mitoCoM(idx).centroid(2),mitoCoM(idx).centroid(3)*zRatio];
        CoMcoordsL = [mitoCoM(idx-1).centroid(1),mitoCoM(idx-1).centroid(2),mitoCoM(idx-1).centroid(3)*zRatio];
%         coordsF = [track(trackNum).WeightedCentroid(idx*2-1),track(trackNum).WeightedCentroid(idx*2)];

        CoMdistance(idx) = norm(CoMcoordsF-CoMcoordsL)*micronPerPixel;
        CoMdisplacement(idx) = norm(CoMcoordsF-CoMcoordsO)*micronPerPixel;

        CoMspeed(idx) = abs(CoMdistance(idx))/secondsPerFrame;
        CoMvelocity(idx) = (CoMdisplacement(idx))/secondsPerFrame;

        CoMxyAngle(idx) = angle0to180(90-(-atan2d((CoMcoordsF(1)-CoMcoordsL(1)),(CoMcoordsF(2)-CoMcoordsL(2)))));
        CoMxzAngle(idx) = angle0to180(90-(-atan2d((CoMcoordsF(1)-CoMcoordsL(1)),(CoMcoordsF(3)-CoMcoordsL(3)))));
        CoMyzAngle(idx) = angle0to180(90-(-atan2d((CoMcoordsF(2)-CoMcoordsL(2)),(CoMcoordsF(3)-CoMcoordsL(3)))));
    end
end

mitoDynamics.distance = distance;
mitoDynamics.displacement = displacement;
mitoDynamics.speed = speed;
mitoDynamics.velocity = velocity;
mitoDynamics.MSD = MSD;
mitoDynamics.xyAngle = xyAngle;
mitoDynamics.xzAngle = xzAngle;
mitoDynamics.yzAngle = yzAngle;

if nargin > 4
    CoMDynamics.distance = CoMdistance;
    CoMDynamics.displacement = CoMdisplacement;
    CoMDynamics.speed = CoMspeed;
    CoMDynamics.velocity = CoMvelocity;
    CoMDynamics.xyAngle = CoMxyAngle;
    CoMDynamics.xzAngle = CoMxzAngle;
    CoMDynamics.yzAngle = CoMyzAngle;
else
    CoMDynamics = struct();
end

end