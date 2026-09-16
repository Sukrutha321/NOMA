function [output,mapping] = bin2gray(x,modulation,order)
%BIN2GRAY Gray encode
%
%   Warning: BIN2GRAY will be removed in a future release. Use the
%   appropriate modulation object or function to remap constellation points 
%   instead. See <a href="matlab:helpview('comm', 'CA_B2G')">this release note</a> for more information.
%
%   Y = BIN2GRAY(X,MODULATION,M) generates a Gray encoded output with the
%   same dimensions as its input parameter X.  The input X can be a scalar,
%   vector, matrix, or 3-D array.  MODULATION is the modulation type, which
%   can be a string equal to 'qam', 'pam', 'fsk', 'dpsk', or 'psk'.  M is
%   the modulation order that must be an integer power of two.
%
%   [Y,MAP] = BIN2GRAY(X,MODULATION,M) generates a Gray encoded output, Y and
%   returns its Gray encoded constellation map, MAP.  The constellation map is a
%   vector of numbers to be assigned to the constellation symbols.
%
%   If you are converting binary coded data to Gray coded data and modulating
%   the result immediately afterwards, you should use the appropriate modulation
%   functions with the 'gray' option, instead of BIN2GRAY.
%
%   EXAMPLE: 
%     % To Gray encode a vector x with a 16-QAM Gray encoded constellation and
%     % return its map, use:
%     x=randi([0 15],1,100);
%     [y,map] = bin2gray(x,'qam',16);
%     
%     % Obtain the symbols for 16-QAM
%     data = 0:15;
%     symbols = qammod(data,16,'bin');
%
%     % Plot the constellation
%     scatterplot(symbols);
%     set(get(gca,'Children'),'Marker','d','MarkerFaceColor','auto');
%     hold on;
%     % Label the constellation points according to the Gray mapping
%     for jj=1:16
%       text(real(symbols(jj))-0.15,imag(symbols(jj))+0.15,...
%       dec2base(map(jj),2,4));
%     end
%     set(gca,'yTick',(-4:2:4),'xTick',(-4:2:4),...
%      'XLim',[-4 4],'YLim',...
%      [-4 4],'Box','on','YGrid','on', 'XGrid','on');
%
%   See also GRAY2BIN, PSKMOD, QAMMOD, PAMMOD, FSKMOD.

%   Copyright 1996-2021 The MathWorks, Inc.

%#codegen

coder.internal.warning('comm:shared:willBeRemovedRef', ...
        'BIN2GRAY', 'CA_B2G');

%% Begin validating inputs

% Typical error checking.
narginchk(3, 3)

%Validate numeric x data
if isempty(x)
    coder.internal.error('comm:bin2gray:InputEmpty');
end

% x must be a scalar, vector, matrix or 3D array
if length(size(x)) > 3
    coder.internal.error('comm:bin2gray:InputDimensions');
end

% x must be a finite non-negative integer
xVec = x(:);
if (any(xVec<0) || any(isinf(xVec)) || (~isreal(x)) || any(floor(xVec) ~= xVec))
    coder.internal.error('comm:bin2gray:InputError');
end

% Validate modulation type
if (~comm.internal.utilities.isCharOrStringScalar(modulation)) || (~strcmpi(modulation,'QAM')) && (~strcmpi(modulation,'PSK'))...
        && (~strcmpi(modulation,'FSK')) && (~strcmpi(modulation,'PAM')) && (~strcmpi(modulation,'DPSK'))
    coder.internal.error('comm:bin2gray:ModulationTypeError');
end

%Validate modulation order
if (order < 2) || (isinf(order) || ...
        (~isreal(order)) || (floor(log2(order)) ~= log2(order)))
    coder.internal.error('comm:bin2gray:ModulationOrderError');
end

% Check for overflows - when x is greater than the modulation order
if (max(xVec) >= order)
    coder.internal.error('comm:bin2gray:XError');
end

%% Start Gray code conversion 
[output, mapping] = comm.internal.utilities.bin2gray(x, modulation, order);
