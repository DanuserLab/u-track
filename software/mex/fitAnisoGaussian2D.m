%FITANISOGAUSSIAN2D Fit a 2-D Anisotropic Gaussian function to data in an image window.
%    [prmVect prmStd C res J] = fitGaussian2D(data, prmVect, mode, options)
%
%    Symbols: xp : x-position
%             yp : y-position
%              A : amplitude
%             sx : standard deviation along the rotated x axis
%             sy : standard deviation along the rotated y axis
%              t : angle
%              c : background
%
%    The origin is defined at the center of the input window.
%
%    Inputs:     data : 2-D image array
%             prmVect : parameter vector with order: [xp, yp, A, sx, sy, t, c]
%                mode : string that defines parameters to be optimized; any among 'xyarstc' (r = sx, s = sy)
%           {options} : vector [MaxIter TolFun TolX]; max. iterations, precision on f(x), precision on x
%
%    Outputs: prmVect : parameter vector
%              prmStd : parameter standard deviations
%                   C : covariance matrix
%                 res : structure with fields
%		    	          .data : residuals
%		    	          .pval : p value of KS test (normally distributed)
%		    	          .mean : mean of residuals
%		    	          .std  : standard deviation of residuals
%		    	          .RSS  : residual sum-of-squares
%                 J : Jacobian
%
% Axis conventions: image processing, see meshgrid
%
% Example: [prmVect prmStd C res J] = fitAnisoGaussian2D(data, [0 0 max(data(:)) 1.5 1.5 pi/6 min(data(:))], 'xyarstc');
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

% (c) Sylvain Berlemont, 2011
