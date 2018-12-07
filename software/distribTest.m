function [confValue]=distribTest(pop1,pop2)
% distribTest returns the percent confidence that two distributions are different using K-S test
%
% here we test whether two populations, pop1 and pop2, are different by
% bootstrapping. first pop1 is used to find a calibrated p-value: we sample
% n1 values (with replacement) twice, and compare the mean-subtracted
% sample distributions using a K-S test.  by doing this many times, we get
% a distribution of p-values. the 5th percentile represents the
% bootstrapped p-value, thresh.  then we perform the K-S test many times
% with pop1 and pop2.  the proportion of resulting p-values smaller than
% thresh represent the confidence that the two populations are in fact
% different.
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


nReps=1000;

n1=length(pop1);
n2=length(pop2);

p1=zeros(nReps,1);       
p2=zeros(nReps,1);
for i=1:nReps
    % sample both populations with replacement, but sample the first twice
    s1a=randsample(pop1,n1,true);
    s1b=randsample(pop1,n1,true);
    
    s2a=randsample(pop2,n2,true);
    
    [h p1(i)]=kstest2(s1a-mean(s1a),s1b-mean(s1b)); % calibration
    [h p2(i)]=kstest2(s1a-mean(s1a),s2a-mean(s2a));
end

thresh=prctile(p1,5);
confValue=100*(sum(p2<thresh)/length(p2));