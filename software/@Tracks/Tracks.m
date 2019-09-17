classdef (Abstract = true) Tracks < handle  & matlab.mixin.Copyable
% Tracks is a class that encapsulates the tracksFinal output of
% TrackingProcess. Each individual Tracks object consists of multiple
% segments that may merge and split. An array of Tracks objects represents
% individual tracks in the same dataset.
%
% This class extends the struct that is the output of trackCloseGapsKalman.
%
% See also trackObj , plotTracks2D,
% common/trackWithGapClosing/postprocessing/misc
% common/trackWithGapClosing/postprocessing/motionAnalysis
% common/trackWithGapClosing/postprocessing/HuetImplementation
%
% Properties
%
% .tracksFeatIndxCG  : Indices of the detections in a movieInfo
%                      structure
%
% .tracksCoordAmpCG3D: The positions and amplitudes of the tracked
%                      features, after gap closing.
%
%                      This 3D matrix has dimensions as follows:
%                      numSegments x 8 x numFrames
%                      * numSegments (rows) = number of track segments in compound
%                        track.
%                      * 8 corresponds to [x, y, z, a, dx, dy, dz, da]
%                      * numFrames = the number of time points in each
%                        segment. Column number is relative to
%                        .seqOfEvents(1) or .startFrame
%
%                      NaN indicates frames where track segments do not exist.
%
%                      Due to linear indexing this is mostly backwards
%                      compatible with the tracksCoordAmpCG property below
%                      and in the original structure when indexed.
%                      e.g. tracksCoordAmpCG3D(:,1:8:end) still accesses the
%                      x coordinate
%                      The incompatibility is when this is referenced
%                      without indexing.
%
% .seqOfEvents      : Matrix with number of rows equal to number
%                     of events happening in a track and 4
%                     columns:
%                     1st: Frame where event happens;
%                     2nd: 1 - start of track, 2 - end of track;
%                     3rd: Index of track segment that ends or starts;
%                     4th: NaN - start is a birth and end is a death,
%                          number - start is due to a split,
%                                   end is due to a merge
%                                   number is the index of track segment
%                                   for the merge/split.
% Dependent Properties
%
% X, Y, Z, A, dX, dY, dZ, dA : xyz coordinates, amplitude and uncertanties
% segmentStartFrame          : Absolute frame where each segment starts
% segmentEndFrame            : Absolute frame where each segment ends
% parentSegment              : Segment from which each segment split.
%                              NaN if the segment originated independently
% spouseSegment              : Segment into which the segment merged
%                              NaN if the segment never merged
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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

% Mark Kittisopikul, January 2015
    properties (Abstract = true)
        % numSegments x numFrames matrix of indices. See class description.
        tracksFeatIndxCG
        % numSegments x 8 x numFrames matrix of coordinates and amplitudes. See class description.
        tracksCoordAmpCG3D
        % numEvents x 4 matrix. See class description.
        seqOfEvents
        % 2D matrix, corresponds to tracksCoordAmpCG3D(:,:)
        tracksCoordAmpCG
        % x coordinates as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,1:8:end)
        % Column is relative to the startFrame
        x
        % Y coordinates as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,2:8:end)
        % Column is relative to the startFrame
        y
        % Z coordinates as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,3:8:end)
        % Column is relative to the startFrame
        z
        % Amplitude as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,4:8:end)
        % Column is relative to the startFrame
        A
        % X uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,5:8:end)
        % Column is relative to the startFrame
        dx
        % Y uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,6:8:end)
        % Column is relative to the startFrame
        dy
        % Z uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,7:8:end)
        % Column is relative to the startFrame
        dz
        % A uncertainty as nSeg by nFrame matrix = tracksCoordAmpCG3D(:,8:8:end)
        % Column is relative to the startFrame
        dA
        % Column vector for the absolute time point each segment started
        segmentStartFrame
        % Column vector for the absolute time point each segmented ended
        segmentEndFrame
        % Column vector for the segment from which the segment originated
        parentSegment
        % Column vector for the segment into which the segment merged
        spouseSegment
        % Scalar absolute frame at which the compound track starts
        startFrame
        % Scalar absolute at which the compound track ends
        endFrame
        % Number of segments in each compound track
        % see also getNumSegments
        numSegments
        % Number of frames in which each compound track exists
        numFrames
    end
    properties
        % Index assigned on construction of object or by reindex()
        index
        % Time, defaults to f-1
        t
    end
    properties ( Dependent = true )
        % Length of sequence from startFrame to endFrame
        lifetime
        % Sequence from startFrame to endFrame
        f
    end
    properties ( Dependent = true , Hidden)
        % Alias for numSegments, deprecated
        nSeg
        % Alias for startFrame, deprecated
        start
    end
    methods
        function l = get.lifetime(obj)
            l = obj.endFrame - obj.startFrame + 1;
        end
        function n = get.nSeg(obj)
            n = obj.numSegments;
        end
        function f = get.f(obj)
            f = obj.startFrame : obj.endFrame;
        end
        function t = get.t(obj)
            if(isempty(obj.t))
                % default to f -1 if empty, backwards-compatible
                t = obj.f - 1;
            else
                % otherwiser return the stored value
                t = obj.t;
            end
        end
        function set.t(obj,t)
            assert(numel(t) == numel(obj.t) || isempty(t), ...
                'Tracks:set_t:timePointLength', ...
                ['The new number of timepoints must be' ...
                 'the same as the number of frames, unless empty.']);
            obj.t = t(:).';
        end
        function s = get.start(obj)
            s = obj.startFrame;
            warning(['Tracks.end and Tracks.start are deprecated.' ...
                    'Use Tracks.endFrame and Tracks.startFrame instead']);
        end
        function varargout = end(obj,k,n)
            if(nargin == 1)
                warning(['Tracks.end and Tracks.start are deprecated.' ...
                    'Use Tracks.endFrame and Tracks.startFrame instead']);
                [varargout{1:length(obj)}] = deal(obj.endFrame);
            else
                varargout{1} = builtin('end',obj,k,n);
            end
        end
        out = textGraph(obj)
        disp(obj)
        plot(obj,varargin)
        m = getMatrix(obj,varargin)
        s = getSparse(obj,varargin)
        n = numTimePoints(obj);
        t = totalSegments(obj)
        b = isstruct(obj)
        s = getStruct(obj)
        seqM = getSeqOfEventsMatrix(obj)
        msM = getMergeSplitMatrix(obj)
        idx = getMergeIdx(obj,msM)
        idx = getSplitIdx(obj,msM)
        oldIdx = reindex(obj,newIdx);
        [merge,split] = getMergeSplitXY(obj,matrix);
        [ obj ] = normalizeSeqOfEvents( obj );
        [ obj ] = setFeatFromIdx(obj,movieInfo);
        [ obj ] = setMovieInfo(obj,movieInfo);
        [ obj ] = addCoord(obj,tracks);
        [ obj ] = addOffset(obj,X,Y,Z);
        [ newTracks ] = getAddCoord(obj,tracks);
        [ obj ] = overlapping(obj,tracks);
        [ newTracks ] = getOverlapping(obj,tracks);
        [ obj ] = multCoord(obj,scalar);
        [ newTracks ] = getMultCoord(obj,scalar);
        [ movieInfo ] = getMovieInfo(obj);
        %mapFramesToTime map frames to timepoints for the collection of Tracks
        [ obj ] = mapFramesToTimepoints( obj, timepoints, frames );
        [ combinedTrack ] = combine( obj );
        [vert,edges,frames,edgesLabel]=getGraph(obj);
    end
end
