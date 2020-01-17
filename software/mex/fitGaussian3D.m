%FITGAUSSIAN3D Fits a 3D Gaussian function to the input volume.
%    [prmVect, prmStd, C, res, J] = fitGaussian3D(volume, prmVect, mode)
%
%    Notation: x : x-position
%              y : y-position
%              z : z-position 
%              A : amplitude
%              s : x,y standard deviation
%              r : z standard deviation
%              c : background
%
%    The origin [0 0 0] is defined as the first index of the input volume.
%
%    Inputs:   
%              volume : 3D image volume
%             prmVect : parameter vector [x y z A s c] -or-
%                       [x y z A s r c] if the z-sampling is anisotropic
%                mode : string that defines parameters to be optimized;
%                       any among 'xyzasrc'
%
%    Optional inputs (after mode):
%     optim. settings : parameter vector with optimization settings: [maxIter eAbs eRel]
%                       max. iterations, tolerances. See GSL documentation.
%
%    Outputs: prmVect : parameter vector
%              prmStd : parameter standard deviations
%                   C : covariance matrix
%                 res : structure with fields:
%                         .data : residuals
%                         .pval : p-value of the Kolmogorov-Smirnov test (normal dist.)
%                         .mean : mean of the residuals
%                         .std  : standard deviation of the residuals
%                         .RSS  : residual sum of squares
%                   J : Jacobian
%
% Note: Pixels in 'data' that are set to NaN are ignored
%
% Axis conventions: image processing, see meshgrid
% For Gaussian mixture fitting, use fitGaussianMixture3D()
%
% Example: [prmVect, prmStd, C, res, J] = fitGaussian3D(vol, [0 0 0 max(vol(:)) 1.5 min(vol(:))], 'xyzAsc');
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

% Francois Aguet, 2013 (last modified Sep 27, 2013)