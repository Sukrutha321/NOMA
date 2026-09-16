% This function generates Gray coded Rectangular M-QAM constellations (not
% normalized)
function [LutKey,kappa,Lambda] = rectGrayQAM_mod(M)
Lambda = 2^(ceil(0.5*log2(M)))-1; %width in x-axis
moh = 2^(floor(0.5*log2(M)))-1;   %width in y-axis
InPhase = -Lambda:2:Lambda;       %in-phase components
InPhase = ones(moh+1,Lambda+1)*diag(InPhase); 
QuadPhase = [moh:-2:-moh].';      %quad-phase components
LUT = InPhase + 1j*QuadPhase;     %Scatter matrix
symbolGray = bin2gray(0:1:M-1,'qam',M); %Gray coding
LutKey = LUT(:).';                %vector LUT
kappa = LutKey*LutKey'/M;         %normalization factor
% LUT = LUT./sqrt(kappa);
LutKey = LutKey(symbolGray+1);    %Gray coded QAM
% LutKey = LutKey./sqrt(kappa);     %Normalized QAM
% symbol = [0:1:M-1].';             %QAM Symbol
% data_LUT = de2bi(symbol,'left-msb'); %QAM Bits
% symbolLabel = num2str(symbol);    %Symbol label
% bitsLabel = dec2base(symbol, 2);  %Bits label
% scatterplot(LutKey)
% text(real(LutKey),imag(LutKey),bitsLabel,'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',6);

