classdef TracksStruct < Tracks & dynamicprops
% TracksStruct is a Tracks implementation that is a close extension of the
% original struct implementation for backwards compatibility while also
% providing the properties and methods of the Tracks interface class.
%
% Emphasis: Backwards-compatability, fast-instantiation
%
% See also Tracks
%
% Mark Kittisopikul, January 2015
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
    properties
        % numSegments x numFrames matrix of indices. See class description.
        tracksFeatIndxCG
        % numSegments x 8 x numFrames matrix of coordinates and amplitudes. See class description.
        tracksCoordAmpCG3D
        % numEvents x 4 matrix. See class description.
        seqOfEvents
    end
    properties (Dependent = true)
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
    methods
        function obj = TracksStruct(tracks,movieInfo)
            % Takes a tracksFinal structure from trackCloseGapsKalman
            if(nargin ~= 0)
                if(~isstruct(tracks))
                    tracks = convertMat2Struct2(tracks);
                end
                obj(numel(tracks)) = TracksStruct();
                obj = reshape(obj,size(tracks));
                if(nargin > 1)
                    % If movieInfo is given, then
                    % extract new coordinates and ignore previous ones.
                    tracksCoordAmpCG = getFeatFromIdx(tracks,movieInfo);
                    [obj.tracksCoordAmpCG] = deal(tracksCoordAmpCG{:});
                else
                    [obj.tracksCoordAmpCG] = deal(tracks.tracksCoordAmpCG);
                end
                [obj.tracksFeatIndxCG] = deal(tracks.tracksFeatIndxCG);
                if(~isfield(tracks,'seqOfEvents'))
                    [obj.seqOfEvents] = deal([]);
                else
                    [obj.seqOfEvents] = deal(tracks.seqOfEvents);
                end
                obj.normalizeSeqOfEvents();
                obj.reindex();
            end
        end
        function tracksCoordAmpCG = get.tracksCoordAmpCG(obj)
            tracksCoordAmpCG = obj.tracksCoordAmpCG3D(:,:);
        end
        function set.tracksCoordAmpCG(obj,tracksCoordAmpCG)
            obj.tracksCoordAmpCG3D = ... 
                reshape(tracksCoordAmpCG,size(tracksCoordAmpCG,1),8,[]);
        end
        function x = get.x(obj)
                x = obj.tracksCoordAmpCG3D(:,1,:);
                x = x(:,:);
        end
        function y = get.y(obj)
                y = obj.tracksCoordAmpCG3D(:,2,:);
                y = y(:,:);
        end
        function z = get.z(obj)
                z = obj.tracksCoordAmpCG3D(:,3,:);
                z = z(:,:);
        end
        function A = get.A(obj)
                A = obj.tracksCoordAmpCG3D(:,4,:);
                A = A(:,:);
        end
        function dx = get.dx(obj)
                dx = obj.tracksCoordAmpCG3D(:,5,:);
                dx = dx(:,:);
        end
        function dy = get.dy(obj)
                dy = obj.tracksCoordAmpCG3D(:,6,:);
                dy = dy(:,:);
        end
        function dz = get.dz(obj)
                dz = obj.tracksCoordAmpCG3D(:,7,:);
                dz = dz(:,:);
        end
        function dA = get.dA(obj)
                dA = obj.tracksCoordAmpCG3D(:,8,:);
                dA = dA(:,:);
        end
        function S = get.segmentStartFrame(obj)
            S = zeros(obj.numSegments,1);
            startIdx = obj.seqOfEvents(:,2) == 1;
            S(obj.seqOfEvents(startIdx,3)) = obj.seqOfEvents(startIdx,1);
        end
        function E = get.segmentEndFrame(obj)
            E = zeros(obj.numSegments,1);
            endIdx = obj.seqOfEvents(:,2) == 2;
            E(obj.seqOfEvents(endIdx,3)) = obj.seqOfEvents(endIdx,1);
            endIdx = endIdx & ~isnan(obj.seqOfEvents(:,4));
            E(obj.seqOfEvents(endIdx,3)) = obj.seqOfEvents(endIdx,1) - 1;
        end
        function O = get.parentSegment(obj)
            O = zeros(obj.numSegments,1);
            idx = obj.seqOfEvents(:,2) == 1;
            O(obj.seqOfEvents(idx,3)) = obj.seqOfEvents(idx,4);
        end
        function D = get.spouseSegment(obj)
            D = zeros(obj.numSegments,1);
            idx = obj.seqOfEvents(:,2) == 2;
            D(obj.seqOfEvents(idx,3)) = obj.seqOfEvents(idx,4);
        end
        function startTime = get.startFrame(obj)
            startTime = obj.seqOfEvents(1,1);
        end
        function set.startFrame(obj,f)
            if(isempty(obj.seqOfEvents))
                obj.seqOfEvents = ...
                    [ 1 1 1 NaN; ...
                      2 2 1 NaN];
            end
            obj.seqOfEvents(1,1) = f;
        end
        function endTime = get.endFrame(obj)
            endTime = obj.seqOfEvents(end,1);
        end
        function set.endFrame(obj,f)
            if(isempty(obj.seqOfEvents))
                obj.seqOfEvents = ...
                    [ 1 1 1 NaN; ...
                      2 2 1 NaN];
            end
            obj.seqOfEvents(end,1) = f;
        end
        function N = get.numSegments(obj)
            N = size(obj.tracksCoordAmpCG3D,1);
        end
        function N = get.numFrames(obj)
            N = size(obj.tracksCoordAmpCG3D,3);
        end
        function P = addprop(obj,propName)
            if(~isscalar(obj))
                P = arrayfun(@(x) addprop(x,propName),obj,'UniformOutput',false);
                P = [P{:}];
                P = reshape(P,size(obj));
            else
                P = addprop@dynamicprops(obj,propName);
            end
        end
    end
end