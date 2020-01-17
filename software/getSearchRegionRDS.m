function [longVecS,longVecE,shortVecS,shortVecE,shortVecS3D,shortVecE3D,...
    longVecSMS,longVecEMS,shortVecSMS,shortVecEMS,shortVecS3DMS,shortVecE3DMS,...
    longRedVecS,longRedVecE,longRedVecSMS,longRedVecEMS] = ...
    getSearchRegionRDS(xyzVelS,xyzVelE,brownStd,...
    trackType,undetBrownStd,timeWindow,brownStdMult,linStdMult,...
    timeReachConfB,timeReachConfL,minSearchRadius,maxSearchRadius,...
    useLocalDensity,closestDistScale,maxStdMult,nnDistLinkedFeat,...
    nnWindow,trackStartTime,trackEndTime,probDim,resLimit,brownScaling,...
    linScaling,linearMotion)
%GETSEARCHREGIONRDS determines the search regions for particles undergoing free diffusion with or without drift with or without direction reversal
%
%SYNOPSIS [longVecS,longVecE,shortVecS,shortVecE,shortVecS3D,shortVecE3D,...
%    longVecSMS,longVecEMS,shortVecSMS,shortVecEMS,shortVecS3DMS,shortVecE3DMS,...
%    longRedVecS,longRedVecE,longRedVecSMS,longRedVecEMS] = ...
%    getSearchRegionRDS(xyzVel,brownStd,trackType,...
%    undetBrownStd,timeWindow,brownStdMult,linStdMult,timeReachConfB,timeReachConfL,...
%    minSearchRadius,maxSearchRadius,useLocalDensity,closestDistScale,...
%    maxStdMult,nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime,...
%    probDim,resLimit,brownScaling,linScaling,linearMotion)
%
%INPUT  xyzVelS        : Velocity in x, y and z (if 3D) at track starts.
%       xyzVelE        : Velocity in x, y and z (if 3D) at track ends.
%       brownStd       : Standard deviation of Brownian motion steps.
%       trackType      : Type of track. 1 for linear, 0 for Brownian, NaN for undetermined.
%       undetBrownStd  : Standard deviation of Brownian motion steps to be used
%                        for undetermined tracks.
%       timeWindow     : Maximum gap size.
%       brownStdMult   : Multiplication factor to go from average Brownian
%                        displacement to search radius.
%       linStdMult     : Multiplication factor to go from average linear
%                        displacement to search radius.
%       timeReachConfB : Time gap for Brownian motion to reach confinement.
%       timeReachConfL : Time gap for linear motion to reach confinement.
%       minSearchRadius: Minimum allowed search radius.
%       maxSearchRadius: Maximum allowed search radius for linking between
%                        two consecutive frames. It will be expanded for
%                        different gap lengths based on the time scaling of
%                        Brownian motion.
%       useLocalDensity: 1 if local density of features is used to expand 
%                        their search radius if possible, 0 otherwise.
%       closestDistScale:Scaling factor of nearest neighbor distance.
%       maxStdMult     : Maximum value of factor multiplying std to get
%                        search radius.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       nnWindow       : Time window to be used in estimating the
%                        nearest neighbor distance of a track at its start
%                        and end.
%       trackStartTime : Starting time of all tracks.            gapCloseParam.timeWindow = 2;
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%       trackEndTime   : Ending time of all tracks.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%       resLimit       : Resolution limit, in whatever space units are
%                        being used.
%       brownScaling   : Power by which Brownian part of search radius
%                        inceases with time.
%       linScaling     : Power by which linear part of search radius
%                        increases with time.
%       linearMotion   : Motion model flag.
%                        0 for only random motion;
%                        1 for random + directed motion;
%                        2 for random + directed motion with the
%                        possibility of instantaneous switching to opposite
%                        direction (but same speed), i.e. something like 1D
%                        diffusion.
%
%OUTPUT longVecS  : Vector defining long axis of search region at track starts.
%       longVecE  : Vector defining long axis of search region at track ends.
%       shortVecS : Vector defining short axis of search region at track starts.
%       shortVecE : Vector defining short axis of search region at track ends.
%       shortVecS3D:Vector defining 2nd short axis of search region at
%                   track starts in case of 3D.
%       shortVecE3D:Vector defining 2nd short axis of search region at
%                   track ends in case of 3D.
%
%       longVecSMS, longVecEMS, shortVecSMS, shortVecEMS, shortVecS3DMS,
%       shortVecE3DMS: Same as above, but for merging ang splitting.
%       
%       longRedVecS,longRedVecE,longRedVecSMS,longRedVecEMS: Same as above,
%                   but for links in opposite direction of velocity in case
%                   linearMotion = 1 (not 2).
%
%Khuloud Jaqaman, September 2011

%% Output

longVecS  = [];
longVecE  = [];
longRedVecS  = [];
longRedVecE  = [];
shortVecS = [];
shortVecE = [];
shortVecS3D = [];
shortVecE3D = [];
longVecSMS  = [];
longVecEMS  = [];
longRedVecSMS  = [];
longRedVecEMS  = [];
shortVecSMS = [];
shortVecEMS = [];
shortVecS3DMS = [];
shortVecE3DMS = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= 24
    disp('--getSearchRegionRDS: Incorrect number of input arguments!');
    return
end

%% Determine search region

%determine number of tracks
numTracks = size(brownStd,1);

%reserve memory for output
longVecS  = zeros(probDim,timeWindow,numTracks);
longVecE  = zeros(probDim,timeWindow,numTracks);
longRedVecS  = zeros(probDim,timeWindow,numTracks);
longRedVecE  = zeros(probDim,timeWindow,numTracks);
shortVecS = zeros(probDim,timeWindow,numTracks);
shortVecE = zeros(probDim,timeWindow,numTracks);
longVecSMS  = zeros(probDim,timeWindow,numTracks);
longVecEMS  = zeros(probDim,timeWindow,numTracks);
longRedVecSMS  = zeros(probDim,timeWindow,numTracks);
longRedVecEMS  = zeros(probDim,timeWindow,numTracks);
shortVecSMS = zeros(probDim,timeWindow,numTracks);
shortVecEMS = zeros(probDim,timeWindow,numTracks);
if probDim == 3
    shortVecS3D = zeros(probDim,timeWindow,numTracks);
    shortVecE3D = zeros(probDim,timeWindow,numTracks);
    shortVecS3DMS = zeros(probDim,timeWindow,numTracks);
    shortVecE3DMS = zeros(probDim,timeWindow,numTracks);
end

%define square root of "problem dimension" to avoid calculating it many times
sqrtDim = sqrt(probDim);

%put time scaling of forward linear motion in a vector
timeScalingLin = [(1:timeReachConfL).^linScaling(1) ...
    (timeReachConfL)^linScaling(1) * (2:timeWindow-timeReachConfL+1).^linScaling(2)];

%put time scaling of Brownian motion in a vector
timeScalingBrown = [(1:timeReachConfB).^brownScaling(1) ...
    (timeReachConfB)^brownScaling(1) * (2:timeWindow-timeReachConfB+1).^brownScaling(2)];

%scale maxSearchRadius like Brownian motion (it's only imposed on the
%Brownian aspect of tracks)
maxSearchRadius = maxSearchRadius * timeScalingBrown;

%calculate minimum and maximum search radii for merging and splitting,
%taking into account the resolution limit
minSearchRadiusMS = max(minSearchRadius,resLimit);
maxSearchRadiusMS = max(maxSearchRadius,resLimit);

%determine the nearest neighbor distances of tracks at their starts and ends
windowLimS = min([trackStartTime+nnWindow trackEndTime],[],2);
windowLimE = max([trackEndTime-nnWindow trackStartTime],[],2);
nnDistTracksS = zeros(numTracks,1);
nnDistTracksE = zeros(numTracks,1);
for iTrack = 1 : numTracks
    nnDistTracksS(iTrack) = min(nnDistLinkedFeat(iTrack,...
        trackStartTime(iTrack):windowLimS(iTrack)));
    nnDistTracksE(iTrack) = min(nnDistLinkedFeat(iTrack,...
        windowLimE(iTrack):trackEndTime(iTrack)));
end

for iTrack = 1 : numTracks
    
    switch trackType(iTrack)

        case 1
            
            %get velocity, its magnitude and "direction of motion"
            %at track start
            velDriftS = xyzVelS(iTrack,:)';
            velMagS = sqrt(velDriftS' * velDriftS);
            directionMotionS = velDriftS / velMagS;
            %at track end
            velDriftE = xyzVelE(iTrack,:)';
            velMagE = sqrt(velDriftE' * velDriftE);
            directionMotionE = velDriftE / velMagE;
            
            %obtain vector(s) perpendicular to direction of motion
            if probDim == 2 %in 2D case, 1 vector needed
                %at track start
                perpendicularS = [-directionMotionS(2) directionMotionS(1)]';
                %at track end
                perpendicularE = [-directionMotionE(2) directionMotionE(1)]';
            else %in 3D case, 2 vectors needed
                %at track start
                perpendicularS = [-directionMotionS(2) directionMotionS(1) 0]';
                perpendicularS = perpendicularS / (sqrt(perpendicularS'*perpendicularS));
                perpendicular3DS = cross(directionMotionS,perpendicularS);
                %at track end
                perpendicularE = [-directionMotionE(2) directionMotionE(1) 0]';
                perpendicularE = perpendicularE / (sqrt(perpendicularE'*perpendicularE));
                perpendicular3DE = cross(directionMotionE,perpendicularE);
            end

            %calculate the expected displacement due to drift for all time
            %gaps
            %at track start
            dispDrift1FS = velMagS * timeScalingLin;
            %at track end
            dispDrift1FE = velMagE * timeScalingLin;
            
            %calculate the expected displacement along x (= along y, [z]) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;
            
            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand Brownian search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance at its start
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track start
                %if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand search radius multiplication factor at track end
                %if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the "long vectors" of the search rectangles for all time
            %gaps when direction of motion is continued
            %at track start
            longVec1FS = ...
                (repmat((linStdMult' .* dispDrift1FS),probDim,1) + ...
                repmat((brownStdMult' .* dispBrown1 * sqrtDim),probDim,1)) .* ...
                repmat(directionMotionS,1,timeWindow);
            longVecFSMag = sqrt((diag(longVec1FS' * longVec1FS))');  %magnitude
            longVecFSDir = longVec1FS ./ repmat(longVecFSMag,probDim,1); %direction
            %at track end
            longVec1FE = ...
                (repmat((linStdMult' .* dispDrift1FE),probDim,1) + ...
                repmat((brownStdMult' .* dispBrown1 * sqrtDim),probDim,1)) .* ...
                repmat(directionMotionE,1,timeWindow);
            longVecFEMag = sqrt((diag(longVec1FE' * longVec1FE))');  %magnitude
            longVecFEDir = longVec1FE ./ repmat(longVecFEMag,probDim,1); %direction
            
            %determine the "long vectors" of the search rectangles for all time
            %gaps when direction of motion is reversed
            %only needed when linearMotion = 1 (when linearMotion = 2, 
            %there is no forward or backward direction)
            if linearMotion == 1
                %at start
                longVec1BS = ...
                    repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                    repmat(directionMotionS,1,timeWindow);
                longVecBSMag = sqrt((diag(longVec1BS' * longVec1BS))');  %magnitude
                longVecBSDir = longVec1BS ./ repmat(longVecBSMag,probDim,1); %direction
                %at end
                longVec1BE = ...
                    repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                    repmat(directionMotionE,1,timeWindow);
                longVecBEMag = sqrt((diag(longVec1BE' * longVec1BE))');  %magnitude
                longVecBEDir = longVec1BE ./ repmat(longVecBEMag,probDim,1); %direction
            end
            
            %determine the "short vectors"
            %at track starts
            shortVecS1 = ...
                repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicularS,1,timeWindow);
            shortVecSMag = sqrt((diag(shortVecS1' * shortVecS1))');  %magnitude
            shortVecSDir = shortVecS1 ./ repmat(shortVecSMag,probDim,1); %direction
            %at track ends
            shortVecE1 = ...
                repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicularE,1,timeWindow);
            shortVecEMag = sqrt((diag(shortVecE1' * shortVecE1))');  %magnitude
            shortVecEDir = shortVecE1 ./ repmat(shortVecEMag,probDim,1); %direction
            
            %make sure that "long vectors" are longer than minimum allowed
            %at start
            longVecSMagTmp = max([longVecFSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            longVec1FS = repmat(longVecSMagTmp,probDim,1) .* longVecFSDir; %new long vector
            %do the same for merging and splitting
            longVecSMagTmp = max([longVecFSMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            longVec1MSFS = repmat(longVecSMagTmp,probDim,1) .* longVecFSDir;
            %at end
            longVecEMagTmp = max([longVecFEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            longVec1FE = repmat(longVecEMagTmp,probDim,1) .* longVecFEDir; %new long vector
            %do the same for merging and splitting
            longVecEMagTmp = max([longVecFEMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            longVec1MSFE = repmat(longVecEMagTmp,probDim,1) .* longVecFEDir;
            
            %if linearMotion=1, make sure that backwards "long vectors" are
            %within allowed range
            if linearMotion == 1
                %at start
                longVecSMagTmp = max([longVecBSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
                longVecSMagTmp = min([longVecSMagTmp;maxSearchRadius]); %compare to maximum
                longVec1BS = repmat(longVecSMagTmp,probDim,1) .* longVecBSDir; %new long vector
                %do the same for merging and splitting
                longVecSMagTmp = max([longVecBSMag;repmat(minSearchRadiusMS,1,timeWindow)]);
                longVecSMagTmp = min([longVecSMagTmp;maxSearchRadiusMS]);
                longVec1MSBS = repmat(longVecSMagTmp,probDim,1) .* longVecBSDir;
                %at end
                longVecEMagTmp = max([longVecBEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
                longVecEMagTmp = min([longVecEMagTmp;maxSearchRadius]); %compare to maximum
                longVec1BE = repmat(longVecEMagTmp,probDim,1) .* longVecBEDir; %new long vector
                %do the same for merging and splitting
                longVecEMagTmp = max([longVecBEMag;repmat(minSearchRadiusMS,1,timeWindow)]);
                longVecEMagTmp = min([longVecEMagTmp;maxSearchRadiusMS]);
                longVec1MSBE = repmat(longVecEMagTmp,probDim,1) .* longVecBEDir;
            end
            
            %make sure that "short vectors" at track starts are within
            %allowed range
            shortVecSMagTmp = max([shortVecSMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVecSMagTmp = min([shortVecSMagTmp;maxSearchRadius]); %compare to maximum
            shortVecS1 = repmat(shortVecSMagTmp,probDim,1) .* shortVecSDir; %new short vector
            %do the same for merging and splitting
            shortVecSMagTmpMS = max([shortVecSMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            shortVecSMagTmpMS = min([shortVecSMagTmpMS;maxSearchRadiusMS]);
            shortVecS1MS = repmat(shortVecSMagTmpMS,probDim,1) .* shortVecSDir;
            
            %make sure that "short vectors" at track ends are within allowed
            %range
            shortVecEMagTmp = max([shortVecEMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            shortVecEMagTmp = min([shortVecEMagTmp;maxSearchRadius]); %compare to maximum
            shortVecE1 = repmat(shortVecEMagTmp,probDim,1) .* shortVecEDir; %new short vector
            %do the same for merging and splitting
            shortVecEMagTmpMS = max([shortVecEMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            shortVecEMagTmpMS = min([shortVecEMagTmpMS;maxSearchRadiusMS]);
            shortVecE1MS = repmat(shortVecEMagTmpMS,probDim,1) .* shortVecEDir;
            
            %save values for this track
            longVecS(:,:,iTrack) = longVec1FS;
            longVecE(:,:,iTrack) = longVec1FE;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;
            longVecSMS(:,:,iTrack) = longVec1MSFS;
            longVecEMS(:,:,iTrack) = longVec1MSFE;
            shortVecSMS(:,:,iTrack) = shortVecS1MS;
            shortVecEMS(:,:,iTrack) = shortVecE1MS;
            switch linearMotion
                case 1
                    longRedVecS(:,:,iTrack) = longVec1BS;
                    longRedVecE(:,:,iTrack) = longVec1BE;
                    longRedVecSMS(:,:,iTrack) = longVec1MSBS;
                    longRedVecEMS(:,:,iTrack) = longVec1MSBE;
                case 2
                    longRedVecS(:,:,iTrack) = longVec1FS;
                    longRedVecE(:,:,iTrack) = longVec1FE;
                    longRedVecSMS(:,:,iTrack) = longVec1MSFS;
                    longRedVecEMS(:,:,iTrack) = longVec1MSFE;
            end
            
            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecS13D = repmat(shortVecSMagTmp,probDim,1) .* repmat(perpendicular3DS,1,timeWindow);
                shortVecE13D = repmat(shortVecEMagTmp,probDim,1) .* repmat(perpendicular3DE,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D; %#ok<AGROW>
                shortVecE3D(:,:,iTrack) = shortVecE13D; %#ok<AGROW>
                %do the same for merging and splitting
                shortVecS13DMS = repmat(shortVecSMagTmpMS,probDim,1) .* repmat(perpendicular3DS,1,timeWindow);
                shortVecE13DMS = repmat(shortVecEMagTmpMS,probDim,1) .* repmat(perpendicular3DE,1,timeWindow);
                shortVecS3DMS(:,:,iTrack) = shortVecS13DMS; %#ok<AGROW>
                shortVecE3DMS(:,:,iTrack) = shortVecE13DMS; %#ok<AGROW>
            end

        case 0
            
            %take direction of motion to be along x and construct
            %perpendicular(s)
            if probDim == 2
                directionMotion = [1 0]';
                perpendicular = [0 1]';
            else
                directionMotion = [1 0 0]';
                perpendicular = [0 1 0]';
                perpendicular3D = [0 0 1]';
            end

            %calculate the expected displacement along x (= along y, [z]) due to
            %brownian motion for all time gaps
            dispBrown1 = brownStd(iTrack) * timeScalingBrown;

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end
            
            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand start's search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance at its end
                %/closestDistScale by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand end's search radius multiplication factor if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end

            %determine the long vectors of the search ellipses at track
            %starts for all time gaps
            longVecS1 = ...
                repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vectors of the search ellipses at track
            %ends for all time gaps
            longVecE1 = ...
                repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vectors at track starts
            shortVecS1 = ...
                repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %determine the short vectors at track ends
            shortVecE1 = ...
                repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is within allowed range
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and splitting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]);
            
            %re-calculate both vectors based on modified magnitudes            
            longVecS1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecS1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir;
            shortVecS1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir;

            %construct additional short vectors for 3D problems and save
            %them
            if probDim == 3
                shortVecS13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D; %#ok<AGROW>
                %repeat for merging and splitting
                shortVecS13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3DMS(:,:,iTrack) = shortVecS13DMS; %#ok<AGROW>
            end

            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))');  %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is within allowed range
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and splitting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]);
            
            %re-calculate both vectors based on modified magnitudes            
            longVecE1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecE1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir;
            shortVecE1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir;

            %construct additional short vectors for 3D problems and save
            %them
            if probDim == 3
                shortVecE13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3D(:,:,iTrack) = shortVecE13D; %#ok<AGROW>
                %repeat for merging and splitting
                shortVecE13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3DMS(:,:,iTrack) = shortVecE13DMS; %#ok<AGROW>
            end

            %save values for this track
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            longRedVecS(:,:,iTrack) = longVecS1;
            longRedVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;
            longVecSMS(:,:,iTrack) = longVecS1MS;
            longVecEMS(:,:,iTrack) = longVecE1MS;
            longRedVecSMS(:,:,iTrack) = longVecS1MS;
            longRedVecEMS(:,:,iTrack) = longVecE1MS;
            shortVecSMS(:,:,iTrack) = shortVecS1MS;
            shortVecEMS(:,:,iTrack) = shortVecE1MS;

        otherwise

            %take direction of motion to be along x and construct
            %perpendicular(s)
            if probDim == 2
                directionMotion = [1 0]';
                perpendicular = [0 1]';
            else
                directionMotion = [1 0 0]';
                perpendicular = [0 1 0]';
                perpendicular3D = [0 0 1]';
            end

            %calculate the expected displacement along x (= along y) due to
            %brownian motion for all time gaps
            if brownStd(iTrack)==1
                dispBrown1 = undetBrownStd * timeScalingBrown;
            else
                dispBrown1 = brownStd(iTrack) * timeScalingBrown;
            end

            %copy brownStdMult into vector that might be modified using
            %local density
            brownStdMultModS = brownStdMult'; %for track start
            brownStdMultModE = brownStdMult'; %for track end

            %if local density information is used to expand search radius ...
            if useLocalDensity

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksS(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand start's search radius multiplication factor if possible
                brownStdMultModS = max([brownStdMultModS; ratioDist2Std]);

                %divide the track's nearest neighbor distance/closestDistScale
                %by expected Brownian displacement
                ratioDist2Std = repmat(nnDistTracksE(iTrack)/closestDistScale,...
                    1,timeWindow) ./ dispBrown1;

                %make ratios larger than maxStdMult equal to maxStdMult
                ratioDist2Std(ratioDist2Std > maxStdMult) = maxStdMult;

                %expand end's search radius multiplication factor if possible
                brownStdMultModE = max([brownStdMultModE; ratioDist2Std]);

            end
            
            %determine the long vector of the search ellipse at track
            %starts for all time gaps
            longVecS1 = ...
                repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the long vector of the search ellipse at track
            %ends for all time gaps
            longVecE1 = ...
                repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(directionMotion,1,timeWindow);

            %determine the short vector at track starts
            shortVecS1 = ...
                repmat((brownStdMultModS .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %determine the short vector at track ends
            shortVecE1 = ...
                repmat((brownStdMultModE .* dispBrown1 * sqrtDim),probDim,1) .* ...
                repmat(perpendicular,1,timeWindow);

            %get magnitude and direction of both vectors at track starts
            vecMag = sqrt((diag(longVecS1' * longVecS1))'); %magnitude of both vectors
            longVecDir = longVecS1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecS1 ./ repmat(vecMag,probDim,1); %direction of short vector

            %make sure that magnitude is within allowed range
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and spltting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]); %compare to minimum
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]); %compare to maximum

            %re-calculate both vectors based on modified magnitudes
            longVecS1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecS1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecS1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir;
            shortVecS1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir;

            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecS13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3D(:,:,iTrack) = shortVecS13D; %#ok<AGROW>
                %repeat for merging and splitting
                shortVecS13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecS3DMS(:,:,iTrack) = shortVecS13DMS; %#ok<AGROW>
            end

            %get magnitude and direction of both vectors at track ends
            vecMag = sqrt((diag(longVecE1' * longVecE1))'); %magnitude of both vectors
            longVecDir = longVecE1 ./ repmat(vecMag,probDim,1);   %direction of long vector
            shortVecDir = shortVecE1 ./ repmat(vecMag,probDim,1); %direction of short vector
            
            %make sure that magnitude is within allowed range
            vecMagTmp = max([vecMag;repmat(minSearchRadius,1,timeWindow)]); %compare to minimum
            vecMagTmp = min([vecMagTmp;maxSearchRadius]); %compare to maximum
            %repeat for merging and spltting
            vecMagTmpMS = max([vecMag;repmat(minSearchRadiusMS,1,timeWindow)]);
            vecMagTmpMS = min([vecMagTmpMS;maxSearchRadiusMS]);
            
            %re-calculate both vectors based on modified magnitudes
            longVecE1 = repmat(vecMagTmp,probDim,1) .* longVecDir; %new long vector
            shortVecE1 = repmat(vecMagTmp,probDim,1) .* shortVecDir; %new short vector
            %repeat for merging and splitting
            longVecE1MS = repmat(vecMagTmpMS,probDim,1) .* longVecDir; %new long vector
            shortVecE1MS = repmat(vecMagTmpMS,probDim,1) .* shortVecDir; %new short vector

            %construct additional short vectors for 3D problems
            if probDim == 3
                shortVecE13D = repmat(vecMagTmp,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3D(:,:,iTrack) = shortVecE13D; %#ok<AGROW>
                %repeat for merging and splitting
                shortVecE13DMS = repmat(vecMagTmpMS,probDim,1) .* repmat(perpendicular3D,1,timeWindow);
                shortVecE3DMS(:,:,iTrack) = shortVecE13DMS; %#ok<AGROW>
            end
            
            %save values for this track
            longVecS(:,:,iTrack) = longVecS1;
            longVecE(:,:,iTrack) = longVecE1;
            longRedVecS(:,:,iTrack) = longVecS1;
            longRedVecE(:,:,iTrack) = longVecE1;
            shortVecS(:,:,iTrack) = shortVecS1;
            shortVecE(:,:,iTrack) = shortVecE1;
            longVecSMS(:,:,iTrack) = longVecS1MS;
            longVecEMS(:,:,iTrack) = longVecE1MS;
            longRedVecSMS(:,:,iTrack) = longVecS1MS;
            longRedVecEMS(:,:,iTrack) = longVecE1MS;
            shortVecSMS(:,:,iTrack) = shortVecS1MS;
            shortVecEMS(:,:,iTrack) = shortVecE1MS;

    end %(switch trackType)

end %(for iTrack = 1 : numTracks)

% longVecS = 10*longVecS;
% longVecE = 10*longVecE;
% shortVecS = 10*shortVecS;
% shortVecE = 10*shortVecE;
% shortVecS3D = 10*shortVecS3D;
% shortVecE3D = 10*shortVecE3D;
% longVecSMS = 10*longVecSMS;
% longVecEMS = 10*longVecEMS;
% shortVecSMS = 10*shortVecSMS;
% shortVecEMS = 10*shortVecEMS;
% shortVecS3DMS = 10*shortVecS3DMS;
% shortVecE3DMS = 10*shortVecE3DMS;
% longRedVecS = 10*longRedVecS;
% longRedVecE = 10*longRedVecE;
% longRedVecSMS = 10*longRedVecSMS;
% longRedVecEMS = 10*longRedVecEMS;


%%%%% ~~ the end ~~ %%%%%

