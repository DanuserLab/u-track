classdef TracksHandle < Tracks & dynamicprops
% TracksHandle is a Tracks implementation optimized for serving the new
% properties such as X, Y, Z while also providing backwards-compatability
% with the tracksFinal struct
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
        % numSegments x numFrames matrix of indices. See class description.
        tracksFeatIndxCG
    end
    properties (Dependent = true)
        % numSegments x 8 x numFrames matrix of coordinates and amplitudes. See class description.
        tracksCoordAmpCG3D
        % 2D matrix, corresponds to tracksCoordAmpCG3D(:,:)
        tracksCoordAmpCG
        % numEvents x 4 matrix. See class description.
        seqOfEvents
        % Number of segments in each compound track
        % see also getNumSegments
        numSegments
        % Number of frames in which each compound track exists
        numFrames
    end
    properties (Access = protected, Transient)
        cache
    end
    methods
        function obj = TracksHandle(tracks,movieInfo)
            % Takes a tracksFinal structure from trackCloseGapsKalman
            if(nargin ~= 0)
                if(~isstruct(tracks))
                    tracks = convertMat2Struct2(tracks);
                end
                obj(numel(tracks)) = TracksHandle();
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
        function set.tracksFeatIndxCG(obj,tracksFeatIndxCG)
            obj.tracksFeatIndxCG = tracksFeatIndxCG;
        end
        function set.tracksCoordAmpCG(obj,tracksCoordAmpCG)
            threeD = reshape(tracksCoordAmpCG,size(tracksCoordAmpCG,1),8,[]);
            obj.tracksCoordAmpCG3D = threeD;
        end
        function tracksCoordAmpCG = get.tracksCoordAmpCG(obj)
            tracksCoordAmpCG = obj.tracksCoordAmpCG3D(:,:);
        end
        function set.tracksCoordAmpCG3D(obj,tracksCoordAmpCG3D)
            obj.x  = tracksCoordAmpCG3D(:,1,:);
            obj.y  = tracksCoordAmpCG3D(:,2,:);
            obj.z  = tracksCoordAmpCG3D(:,3,:);
            obj.A  = tracksCoordAmpCG3D(:,4,:);
            obj.dx = tracksCoordAmpCG3D(:,5,:);
            obj.dy = tracksCoordAmpCG3D(:,6,:);
            obj.dz = tracksCoordAmpCG3D(:,7,:);
            obj.dA = tracksCoordAmpCG3D(:,8,:);
            
            obj.x  = obj.x(:,:);
            obj.y  = obj.y(:,:);
            obj.z  = obj.z(:,:);
            obj.A  = obj.A(:,:);
            obj.dx = obj.dx(:,:);
            obj.dy = obj.dy(:,:);
            obj.dz = obj.dz(:,:);
            obj.dA = obj.dA(:,:);
        end
        function tracksCoordAmpCG3D = get.tracksCoordAmpCG3D(obj)
            tracksCoordAmpCG3D = zeros(obj.numSegments,8,obj.numFrames);
            tracksCoordAmpCG3D(:,1,:) = obj.x;
            tracksCoordAmpCG3D(:,2,:) = obj.y;
            tracksCoordAmpCG3D(:,3,:) = obj.z;
            tracksCoordAmpCG3D(:,4,:) = obj.A;
            tracksCoordAmpCG3D(:,5,:) = obj.dx;
            tracksCoordAmpCG3D(:,6,:) = obj.dy;
            tracksCoordAmpCG3D(:,7,:) = obj.dz;
            tracksCoordAmpCG3D(:,8,:) = obj.dA;
        end
        function set.seqOfEvents(obj,seqOfEvents)
            init = zeros(obj.numSegments,1);
            
            obj.segmentStartFrame = init;
            obj.segmentEndFrame = init;
            obj.parentSegment = init;
            obj.spouseSegment = init;
            
            if(isempty(seqOfEvents))
                % Reset to start at frame 1
                obj.segmentStartFrame(:) = 1;
                obj.segmentEndFrame(:) = obj.numFrames;
                obj.parentSegment(:) = NaN;
                obj.spouseSegment(:) = NaN;
                return;
            end
            
            startIdx = seqOfEvents(:,2) == 1;
            obj.segmentStartFrame(seqOfEvents(startIdx,3)) = seqOfEvents(startIdx,1);
            
            
            endIdx = seqOfEvents(:,2) == 2;
            obj.segmentEndFrame(seqOfEvents(endIdx,3)) = seqOfEvents(endIdx,1);
            endIdx = endIdx & ~isnan(seqOfEvents(:,4));
            obj.segmentEndFrame(seqOfEvents(endIdx,3)) = seqOfEvents(endIdx,1) - 1;
            
            
            idx = seqOfEvents(:,2) == 1;
            obj.parentSegment(seqOfEvents(idx,3)) = seqOfEvents(idx,4);
            
            
            idx = seqOfEvents(:,2) == 2;
            obj.spouseSegment(seqOfEvents(idx,3)) = seqOfEvents(idx,4);
            
            obj.startFrame = seqOfEvents(1,1);
            obj.endFrame = seqOfEvents(end,1);
        end
        function seqOfEvents = get.seqOfEvents(obj)
            seqOfEvents = zeros(obj.numSegments*2,4);
            
            startIdx = 1:obj.numSegments;
            seqOfEvents(startIdx,1) = obj.segmentStartFrame;
            seqOfEvents(startIdx,2) = 1;
            seqOfEvents(startIdx,3) = 1:obj.numSegments;
            seqOfEvents(startIdx,4) = obj.parentSegment;
            
            % For segments that merge, their true endFrames are offset by 1
            E = obj.segmentEndFrame;
            segmentMerged = ~isnan(obj.spouseSegment);
            E(segmentMerged) = E(segmentMerged) + 1;
            
            endIdx = (1:obj.numSegments) + obj.numSegments;
            seqOfEvents(endIdx,1) = E;
            seqOfEvents(endIdx,2) = 2;
            seqOfEvents(endIdx,3) = 1:obj.numSegments;
            seqOfEvents(endIdx,4) = obj.spouseSegment;
            
            seqOfEvents = sortrows(seqOfEvents);
        end
        function resetCache(obj)
            c = struct();
            c.tracksCoordAmpCG = [];
            obj.cache = c;
        end
        function N = get.numSegments(obj)
            N = size(obj.x,1);
        end
        function N = get.numFrames(obj)
            N = size(obj.x,2);
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
%         function out = subsasgn(A,S,B,varargin)
%             try
%                 if(isempty(A))
%                     A = TracksHandle.empty;
%                 end
%                 out = builtin('subsasgn',A,S,B,varargin{:});
%             catch err
%                 switch(err.identifier)
%                     case 'MATLAB:noPublicFieldForClass'
%                         if(~all(isprop(A,S(1).subs)))
%                             rethrow(err);
%                         end
%                         if(nargin < 4)
%                             % Allow for [tracks.prop] = 5; for dynamic
%                             % properties
%                             out = arrayfun(@(t) subsasgn(t,S,B),A,'UniformOutput',false);
%                             out = [out{:}];
%                             out = reshape(out,size(A));
%                         else
%                             % Allow for
%                             % test = {1,2,3,4,5}
%                             % [tracks.prop] = test{:}
%                             % for dynamic properties
%                             out = arrayfun(@(t,b) subsasgn(t,S,b{1}),A,[{B} varargin],'UniformOutput',false);
%                             out = [out{:}];
%                             out = reshape(out,size(A));
%                         end
%                     otherwise
%                         rethrow(err)
%                 end
%             end
%         end
%         function varargout = subsref(A,S)
%             try
%                 [varargout{1:nargout}] = builtin('subsref',A,S);
%             catch err
%                 switch(err.identifier)
%                     case 'MATLAB:noSuchMethodOrField'
%                         if(all(isprop(A,S(1).subs)))
%                             % Allow for tracks.prop where prop is a dynamic
%                             % property
%                             varargout = arrayfun(@(t) subsref(t,S),A,'Unif',false);
%                         else
%                             rethrow(err);
%                         end
%                     otherwise
%                         rethrow(err);
%                 end
%             end
%         end
    end
end
