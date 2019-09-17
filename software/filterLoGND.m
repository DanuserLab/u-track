function [ imLoGResponse, varargout ] = filterLoGND( imInput, sigma, varargin )
% filterLoGND: An ND implementation of the Laplacian of Gaussian (LoG) Filter
% 
%     [ imLoGResponse ] = filterLoGND(im, sigma, varargin );
% 
%     Required Input Arguments:
% 
%                            im: Input ND Image
%         
%                         sigma: standard deviation of the gaussian
%                                should be a scalar value
% 
%                                Caution: If the pixel spacing is not 
%                                specified then unit isotropic spacing is 
%                                assumed as a result of which the units of 
%                                sigma will be assumed to be pixels. So if 
%                                you want to provide sigma in physical
%                                space units then you should also specify
%                                the pixel spacing.
% 
%     Optional Input Arguments:
% 
%                       spacing: pixel spacing of the input image.
%                                can be a scalar value or an array of size
%                                equal to the number of image dimensions.
%                                Default: 1 (isotropic spacing)                      
%                 
%               borderCondition: specifies the way in which the image is 
%                                padded at the borders. This argument is
%                                supplied as input to the functuin 'padarrayXT'. 
%                                Default: 'symmetric'
%                                Options: 'symmetric', 'replicate', 'circular', 
%                                         'antisymmetric', or a constant value
% 
%      UseNormalizedDerivatives: true/false
%                                specifies whether the gaussian derivatives
%                                should be scale-normalized.
%                                Default: false
%                                    
%         UseNormalizedGaussian: true/false
%                                specifies whether the guassian kernel should
%                                be normalized or not
%                                Default: true
%                                   
%                        UseGPU: true/false
%                                A flag that specifies whether or not to
%                                use the GPU for convolution. % 
%                                   True - Uses GPU if installed
%                                   False - otherwise (default-value)
%                                       
%     Output Arguments:
% 
%                 imLoGResponse: LoG filtered image
%         
%                     logKernel: returns the LoG kernel used
%                                (optional output argument)   
% 
%     Example-1: Application of LoG to a 1D Step Edge
%         
%         step = @(x) ( double( x >= 0 ) ); % 1-D step edge
%         g = @(x,sigma) ( (sigma * sqrt(2*pi))^-1 * exp( -x.^2 / (2*sigma^2) ) ); % gaussian
%         truelogresponse = @(x,sigma) ( (-x/sigma^2) .* g(x,sigma) ); % true LoG response
%         
%         sampleSpacing = 0.01;
%         x = -20:sampleSpacing:20;
%         sigma = 2;
%         
%         figure; hold all;
%         
%         % plot step edge
%         plot( x, step(x), '-' );
%         
%         % plot true scale-normalized (multiply by sigma^2) LoG response
%         plot( x, sigma^2 * truelogresponse(x,sigma), '-' ); 
%         
%         % plot response of discrete LoG response 
%         plot( x, filterLoGND( step(x), sigma, ...
%                               'spacing', sampleSpacing, ...
%                               'UseNormalizedDerivatives', true ), '-' ); 
% 
%         title( 'Comparison of our filter with true LoG response for a 1D Sted Edge', 'FontWeight', 'bold' );
%         legend( '1D Step Edge', sprintf( 'True LoG Response (\\sigma = %d)', sigma ), 'Our  LoG Response' );          
%         
%     Example-2: Application of LoG to a flat 1D Ridge
%         
%         step = @(x) ( double( x >= 0 ) ); % 1-D step edge
%         ridge = @(x,w) ( step(x + w/2) - step(x - w/2) ); % 1-D ridge of width w
%         g = @(x,sigma) ( (sigma * sqrt(2*pi))^-1 * exp( -x.^2 / (2*sigma^2) ) ); % gaussian
%         dg = @(x,sigma) ( (-x/sigma^2) .* g(x,sigma) ); % first derivative of gaussian
%         d2g = @(x,sigma) ( ((x.^2 - sigma^2)/sigma^4) .* g(x,sigma) ); % second-derivative (laplacian) of gaussian
%         truelogresponse = @(x,w,sigma) ( dg(x + w/2,sigma) - dg(x - w/2, sigma) ); % true LoG response for ridge of width w
%         
%         w = 10; % width of ridge
%         sampleSpacing = w/100;
%         x = -4*w:sampleSpacing:4*w; 
%         datadims = 1;
% 
%         ridgesignal = ridge(x,w);
%         
%         % check the difference between the response of this filter and the true LoG response
%         figure; hold all;        
%         sigma = w/(2*sqrt(datadims)); % this gives maximal reponse at ridge-center
%         
%             % plot ridge
%             plot( x, ridgesignal, '-' );
% 
%             % plot scale-normalized (multiply by sigma^2) log kernel/filter
%             plot( x, -sigma^2 * d2g(x,sigma), '-' );
% 
%             % plot true scale-normalized (multiply by sigma^2) LoG response
%             plot( x, -sigma^2 * truelogresponse(x,w,sigma), '-' ); 
% 
%             % plot response of discrete LoG response 
%             plot( x, -filterLoGND( ridgesignal, sigma, ...
%                                    'spacing', sampleSpacing, ...
%                                    'UseNormalizedDerivatives', true ), '-' ); 
%         
%             title( 'Comparison between the response of our filter and true LoG response for a 1D Ridge', 'FontWeight', 'bold' );
%             legend( sprintf( '1D Ridge (width = %d)', w ), ... 
%                     sprintf( 'Normalized LoG Kernel (\\sigma = width/%.2f)', w/sigma), ...
%                     'True LoG Response', ...
%                     'Our  LoG Response' );          
%             
%         % check the difference between normalized and unnormalized gaussian derivatives
%         figure;
%         
%             % plot responses for unnormalize gaussian derivatives
%             subplot(2,1,1); hold all;
%             strLegend = {};
%             
%                 % plot ridge - scale to reduce height of ridge so the tiny
%                 % unnormalized LoG responses become visible
%                 plot( x, (1/500) * ridgesignal, '-' ); 
% 
%                 strLegend{end+1} = sprintf( '0.002 \\times 1D Ridge (width = %d)', w );
%                 
%                 % plot response of discrete LoG response 
%                 sigmaValues = (w/(2*sqrt(datadims))) * 2.^[-2:1];
%                 for i = 1:numel( sigmaValues )
%                     
%                     plot( x, -filterLoGND( ridgesignal, sigmaValues(i), ...
%                                            'spacing', sampleSpacing, ...
%                                            'UseNormalizedDerivatives', false ), '-' ); 
%                                        
%                     strLegend{end+1} = sprintf( 'LoG Response for \\sigma = w/%.2f', w/sigmaValues(i) );
%                 end
% 
%             title( 'LoG Responses of a 1D Ridge for Unnormalized Gaussian Derivatives', 'FontWeight', 'bold' );
%             legend( strLegend );          
% 
%             % plot responses for scale-normalized gaussian derivatives
%             subplot(2,1,2); hold all;
%             strLegend = {};
%             
%                 % plot ridge
%                 plot( x, ridgesignal, '-' );
% 
%                 strLegend{end+1} = sprintf( '1D Ridge (width = %d)', w );
%                 
%                 % plot response of discrete LoG response 
%                 sigmaValues = (w/(2*sqrt(datadims))) * 2.^[-2:1];
%                 for i = 1:numel( sigmaValues )
%                     
%                     plot( x, -filterLoGND( ridgesignal, sigmaValues(i), ...
%                                            'spacing', sampleSpacing, ...
%                                            'UseNormalizedDerivatives', true ), '-' ); 
%                                        
%                     strLegend{end+1} = sprintf( 'LoG Response for \\sigma = w/%.2f', w/sigmaValues(i) );
%                 end
% 
%             title( 'LoG Responses of a 1D Ridge for Scale-normalized Gaussian Derivatives', 'FontWeight', 'bold' );
%             legend( strLegend );          
%             
%         
%     Author: Deepak Roy Chittajallu
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
    p.CaseSensitive( false );
    p.addRequired( 'imInput', @(x) ( isnumeric(x) ) );    
    p.addRequired( 'sigma', @(x) ( isscalar(x) ) ); 
    p.parse( imInput, sigma );    
    
    p.addParamValue( 'spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(imInput)])) );
    p.addParamValue( 'borderCondition', 'symmetric', @(x) ( isscalar(x) || (ischar(x) && ismember( lower(x), { 'symmetric', 'replicate', 'antisymmetric' } ) ) ) );
    p.addParamValue( 'UseNormalizedGaussian', true, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'UseNormalizedDerivatives', false, @(x) (isscalar(x) && islogical(x)) );
    p.addParamValue( 'UseGPU', false, @(x) (isscalar(x) && islogical(x)) );
    p.parse( imInput, sigma, varargin{:} );

    spacing = p.Results.spacing;
    borderCondition = p.Results.borderCondition;
    flagNormalizeDerivatives = p.Results.UseNormalizedDerivatives;
    flagNormalizeGaussian = p.Results.UseNormalizedGaussian;
    flagUseGPU = p.Results.UseGPU;
    
    dims = ndims(imInput);   
    if numel(imInput)==max(size(imInput))
        dims = 1;
    end        

   sigma = sigma * ones( 1, dims );

    % adjust sigma according to pixel spacing
    sigma = sigma ./ spacing;
    
    % Compute LoG kernel
    w = ceil(4 * sigma);
    xrange = cell(1,dims);
    
    for i = 1:dims
       xrange{i} = -w(i):w(i); 
    end
    
    x = cell(1,dims);
    
    if dims > 1
        [x{:}] = ndgrid( xrange{:} );    
    else
        x{1} = xrange{1};
    end
    
    G = ones( size(x{1}) );
    C = zeros( size(x{1}) );
    for i = 1:dims
        
        x2 = x{i}.^2;
        
        G = G .* ( exp(-x2 / (2*(sigma(i))^2)) );
        
        if flagNormalizeDerivatives
            % use normalized derivatives -- essentially multiply by sigma(i)^2
            % the reason for doing this is to counter the fact that the 
            % amplitude of the gaussian smoothed derivatives decreases as
            % the sigma increases
            C = C + ((x2 - (sigma(i))^2) / (sigma(i))^2); 
        else
            C = C + ((x2 - (sigma(i))^2) / (sigma(i))^4);
        end
    end    

    if flagNormalizeGaussian
        G = G / sum( G(:) );
    end
    
    logKernel = C .* G;
    logKernel = logKernel - mean( logKernel(:) );
    
    % apply in fourier domain
    if dims > 1
        padsize = w;
    else
        padsize = [0, w];
    end
    imPadded = padarrayXT(imInput, padsize, borderCondition);
    imLoGResponse = convnfft( imPadded, logKernel, 'valid', 'UseGPU', flagUseGPU);
    
    % return logKernel if requested
    if nargout > 1
        varargout{1} = logKernel;
    end
    
end