function C = convnfft(A, B, varargin)
% CONVNFFT  FFT-BASED N-dimensional convolution.
%   C = CONVNFFT(A, B) performs the N-dimensional convolution of
%   matrices A and B. If nak = size(A,k) and nbk = size(B,k), then
%   size(C,k) = max([nak+nbk-1,nak,nbk]);
% 
%   C = CONVNFFT(A, B, SHAPE) controls the size of the answer C:
%       'full'   - (default) returns the full N-D convolution
%       'same'   - returns the central part of the convolution that
%                  is the same size as A.
%       'valid'  - returns only the part of the result that can be
%                  computed without assuming zero-padded arrays.
%                  size(C,k) = max([nak-max(0,nbk-1)],0).
%
%   C = CONVNFFT(..., SHAPE, DIMS) with DIMS is vector of dimensions where
%       the convolution will be carried out. By default DIMS is
%       [1:max(ndims(A),ndims(B))] (all dimensions). A and B must have the
%       same lengths on other dimensions.
%   C = CONVNFFT(..., SHAPE, DIMS, GPU)
%       GPU is boolean flag, see next
%
%   C = CONVNFFT(..., SHAPE, DIMS, ...)
%
%   Optional Arguments:
%     
%         UsePowerOfTwo: true/false
%         
%             rounds the dimension to the nearest power of 2 by zero-padding 
%             while doing the FFT. It is faster but requires more memory.
%             Default-value: true
%
% Class support for inputs A,B:
% float: double, single
%
% METHOD: CONVNFFT uses Fourier transform (FT) convolution theorem, i.e.
%         FT of the convolution is equal to the product of the FTs of the
%         input functions.
%         In 1-D, the complexity is O((na+nb)*log(na+nb)), where na/nb are
%         respectively the lengths of A and B.
%
% Usage recommendation:
%         In 1D, this function is faster than CONV for nA, nB > 1000.
%         In 2D, this function is faster than CONV2 for nA, nB > 20.
%         In 3D, this function is faster than CONVN for nA, nB > 5.
% 
% See also conv, conv2, convn.
% 
%   Author: Deepak Roy Chittajallu
%  
%       
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

    p = inputParser;

    p.addRequired( 'A', @isnumeric );
    p.addRequired( 'B', @isnumeric );
    p.parse(A, B);

    nd = max(ndims(A),ndims(B));

    p.addOptional( 'shape', 'full', @(x) (ischar(x) && ismember(x, {'full', 'same', 'valid'}) ) );
    p.addOptional( 'dims', 1:nd, @(x) (isnumeric(x) && ~any(x < 1 || x > nd)) );
    p.addParamValue( 'UseGPU', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'UsePowerOfTwo', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse(A, B, varargin{:});

    shape = p.Results.shape;
    dims = p.Results.dims;
    flagUseGPU = p.Results.UseGPU;
    flagPower2 = p.Results.UsePowerOfTwo;

    if flagUseGPU
        flagUseGPU = false; % for not this not supported
    end
        
    dims = reshape(unique(dims), 1, []); % row (needed for for-loop index)

    % IFUN function will be used later to truncate the result
    % M and N are respectively the length of A and B in some dimension
    switch lower(shape)
        case 'full',
            ifun = @(m,n) 1:m+n-1;
        case 'same',
            ifun = @(m,n) ceil((n-1)/2)+(1:m);
        case 'valid',
            ifun = @(m,n) n:m;
        otherwise
            error('convnfft: unknown shape %s', shape);
    end

    ABreal = isreal(A) && isreal(B);

    % make dimension a power of 2 if needed
    if flagPower2
        % faster FFT if the dimension is power of 2
        lfftfun = @(l) 2^nextpow2(l);
    else
        % slower, but smaller temporary arrays
        lfftfun = @(l) l;
    end

    % Compute the FFTs
    if flagUseGPU

        % pad arrays with zeros
        subsA(1:ndims(A)) = {':'};
        subsB(1:ndims(B)) = {':'};
        for dim=dims
            m = size(A,dim);
            n = size(B,dim);
            l = lfftfun(m+n-1);               
            
            if l < m
                subsA(1:ndims(A)) = {':'};
                subsA{dim} = 1:l;
                A = A(subsA{:});
            elseif l > m
                subsA(1:ndims(A)) = {':'};
                subsA{dim} = m+1:l;
                A(subsA{:}) = 0;
            end
            
            if l < n
                subsB(1:ndims(B)) = {':'};
                subsB{dim} = 1:l;
                B = B(subsB{:});
            elseif l > n
                subsB(1:ndims(B)) = {':'};
                subsB{dim} = m+1:l;
                B(subsB{:}) = 0;
            end
        end
        
        % make gpy arrays
        A = gpuArray(A);
        B = gpuArray(B);
        
        % do the fft
        subs(1:ndims(A)) = {':'};
        for dim=dims
            % We need to swap dimensions because GPU FFT works along the first dimension
            if dim~=1 % do the work when only required
                swap = 1:nd;
                swap([1 dim]) = swap([dim 1]);
                A = permute(A, swap);
                B = permute(B, swap);
            end
            A = fft(A);
            B = fft(B);
            subs{dim} = ifun(m,n);
        end
        
    else
        
        subs(1:ndims(A)) = {':'};
        for dim=dims

            m = size(A,dim);
            n = size(B,dim);
            
            l = lfftfun(m+n-1);
            
            A = fft(A,l,dim);
            B = fft(B,l,dim);
            subs{dim} = ifun(m,n);

        end            
        
    end    

    % multiply the ffts of A and B element-wise 
    C = A .* B;

    % translate C back from frequency domain
    for dim=dims
        C = ifft(C,[],dim);
    end

    % Truncate the results
    if ABreal
        % Make sure the result is real
        C = real(C(subs{:}));
    else
        C = C(subs{:});
    end
            
    if flagUseGPU
        C = gather(C);
    end
    
end

