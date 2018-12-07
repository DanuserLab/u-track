function pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree)
%MMFDISTPV calculates the p-values of inter-particle distances using a t-test
%
%SYNOPSIS pValue = mmfDistPV(maximaPos,varCovMat,numMaxima,numDegFree)
%
%INPUT  maximaPos  : Particle positions.
%       varCovMat  : Variance-covariance matrix from fit.
%       numMaxima  : Number of particles.
%       numDegFree : Number of degrees of freedom for t-test.
%       PLEASE SEE detectSubResFeatures2D_V2 FOR PROPER CONTEXT.
%
%OUTPUT p-value    : Matrix of inter-particle distance p-values.
%
%REMARKS This function in its current format is not really written for
%general use but as a function to be called by detectSubResFeatures2D_V2.
%
%Khuloud Jaqaman, August 2011
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
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

%reserve memory for output variable
pValue = zeros(numMaxima);

%go over all pairs of particles
for k = 1 : numMaxima-1
    for j = k+1 : numMaxima
        
        %calculate distance between the 2 particles
        x1_x2 = maximaPos(j,1) - maximaPos(k,1);
        y1_y2 = maximaPos(j,2) - maximaPos(k,2);
        distance = sqrt(x1_x2^2+y1_y2^2);
        
        %get the standard deviation of the distance
        j1 = 3*(j-1)+1;
        k1 = 3*(k-1)+1;
        stdDist = x1_x2^2*(varCovMat(j1,j1) + ...
            varCovMat(k1,k1) - 2*varCovMat(j1,k1)) ...
            + y1_y2^2*(varCovMat(j1+1,j1+1) + ...
            varCovMat(k1+1,k1+1) - 2*varCovMat(j1+1,k1+1)) ...
            + 2*x1_x2*y1_y2*(varCovMat(j1,j1+1) - ...
            varCovMat(j1,k1+1) - varCovMat(j1+1,k1) + ...
            varCovMat(k1,k1+1));
        stdDist = sqrt(stdDist)/distance;
        
        %1-sided t-test: H0: T=0, H1: T>0
        %calculate test statistic (t-distributed)
        testStat = distance/stdDist;
        
        %get p-value
        pValue(j,k) = 1-tcdf(testStat,numDegFree);
        
    end
end

%% ~~~ the end ~~~