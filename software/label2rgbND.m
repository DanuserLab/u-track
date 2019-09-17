function [imLabelRGB, varargout] = label2rgbND( imLabel, labelColorMap )

    imLabel = double(imLabel);
    
    if ~exist( 'labelColorMap', 'var' )

        % assign colors randomly
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
        numLabels = double(max(imLabel(:)));
        imLabel( imLabel == 0 ) = numLabels + 1;
        
        cmap = random('Uniform', 0.00001, 1, [numLabels, 3]);
        cmap = [ cmap ; 0 0 0 ];
        
    else
        assert( all(size(labelColorMap, 2) == 3) );
        assert( all(size(labelColorMap, 1) >= max(imLabel(:))+1) );                
        
        numLabels = size(labelColorMap, 1) - 1; 
        imLabel( imLabel == 0 ) = numLabels + 1;
        
        cmap = labelColorMap;
    end
    
    imLabelRGB = [];
    for i = 1:3
        imLabelRGB = cat( ndims(imLabel) + 1, imLabelRGB, reshape( cmap( imLabel(:), i ), size( imLabel ) ) );
    end
    
    if nargout > 1
        varargout{1} = cmap;
    end
    
end