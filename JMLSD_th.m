function [BERth,SqERth] = JMLSD_th(SNR_dB,a,M,T,channel,m_factor)
% clc 
% clear
% close all
% SNR_dB = 0:35;
% a = [0.1 0.9];
% M = [2 2];
% T = [1 2];
% channel = 'Nak';
% m_factor = [4 4];
SNR = 10.^(SNR_dB./10);
N = length(M);
maxT = max(T);
minT = min(T);
tau = maxT./T;
LUT = cell(N,maxT);
symLUT = cell(N,maxT);
if(strcmp(channel,'AWGN'))
    PL = 10.^([0 0]./10);
else
    PL = 10.^([0 6]./10);
end



for i = 1:N
    [LutKey{i},kappa(i),Lambda(i)] = rectGrayQAM_mod(M(i));
    LutKey{i} = LutKey{i}./sqrt(kappa(i));
    symKey{i} = 0:M(i)-1;
    tempSymMesh = cell(1,tau(i));
    tempLut_o = cell(1,tau(i));
    [tempLut_o{:}] = ndgrid(LutKey{i});
    [tempSymMesh{:}] = ndgrid(symKey{i});
    if(tau(i)==maxT)
        for j = 1:length(tempLut_o)
            LutKeyT{i,j} = tempLut_o{j}(:);
            symKeyT{i,j} = tempSymMesh{j}(:);
        end
    elseif(tau(i)==minT)
        for j = 1:maxT
            LutKeyT{i,j} = LutKey{i}.';
            symKeyT{i,j} = symKey{i}.';
        end
    else
        for j = 1:maxT
            LutKeyT{i,j} = tempLut_o{ceil(j./length(tempLut_o))}(:);
            symKeyT{i,j} = tempSymMesh{ceil(j./length(tempLut_o))}(:);
        end
    end
    
end

LUT_C = cell(1,maxT);
data_LUT = cell(1,maxT);
for i = 1:maxT
    temp = cell(N,1);
    [temp{:}] = LutKeyT{:,i};
    tempSym = cell(N,1);
    [tempSym{:}] = symKeyT{:,i};
    [LUT{:,i}] = ndgrid(temp{:});
    [symLUT{:,i}] = ndgrid(tempSym{:});
    LUT_C{i} = 0;
    for j = 1:N
        LUT{j,i} = LUT{j,i}(:);
        symLUT{j,i} = symLUT{j,i}(:);
        data_LUT{j,i} = de2bi(symLUT{j,i},'left-msb');
        LUT_C{i} = LUT_C{i} + sqrt(a(j)).*LUT{j,i};
    end
end


LUT_C = cell2mat(LUT_C);
% LUT_C = cell2sym(LUT_C)
% tab1 = table(LUT_C(:,1), LUT_C(:,2), LUT{1,1}, LUT{1,2}, LUT{2,1}, LUT{2,2}, symLUT{1,1}, symLUT{1,2}, symLUT{2,1}, symLUT{2,2});
% tab1 = renamevars(tab1,1:width(tab1),["x[1]","x[2]","x_1[1]","x_1[2]","x_2[1]","x_2[2]","d_1[1]","d_[2]","d_2[1]","d_2[2]"])
%%
for i = 1:size(LUT_C,1)
    for j = 1:size(LUT_C,1)
%         Ed(i,j) = sum(abs(LUT_C(i,:) - LUT_C(j,:)),2);
        Ed(i,j) = norm(LUT_C(i,:) - LUT_C(j,:),2);
        Ed_real(i,j) = norm(real(LUT_C(i,:) - LUT_C(j,:)),2);

    end
end

for i = 1:N
    dataLUT{i} = [data_LUT{i,:}];
    for m = 1:size(dataLUT{i},2)
        tempDataLUT = dataLUT{i}(:,m);
        for l = 1:size(LUT_C,1)
            for k = 1:size(LUT_C,1)
                bitFlip_perBit{i,m}(l,k) = logical(tempDataLUT(l)~= tempDataLUT(k));
            end
        end
        Sd(i,:) = sqrt(PL(i)./(2*SNR));
        d_new{i,m} = Ed(bitFlip_perBit{i,m});
        d_real{i,m} = Ed_real(bitFlip_perBit{i,m});
        gamma_new{i,m} = (d_new{i,m}/2).^2./Sd(i,:).^2;
%         gamma_new{i,m} = (d_new{i,m}/2).^2./d_real{i,m}./Sd(i,:).^2;
        if(strcmp(channel,'Ray'))
            PHI_new{i,m} = 0.5.*(1 - sqrt(gamma_new{i,m}./(2 + gamma_new{i,m})));
            PEPth_new{i}(m,:) = 1./size(LUT_C,1).*sum(PHI_new{i,m},1);
        elseif(strcmp(channel,'Nak'))
            mu_new{i,m} = sqrt(gamma_new{i,m}./(2.*m_factor(i) + gamma_new{i,m}));
            k_nak{i} = [0:1:m_factor(i)-1].';
            k_nak{i} = reshape(k_nak{i},[1,1,length(k_nak{i})]);
            PHI_new{i,m} = mu_new{i,m}.*sum(round(exp(gammaln(2.*k_nak{i}+1)-gammaln(k_nak{i}+1)-gammaln(k_nak{i}+1))).*((1 - mu_new{i,m}.^2)./4).^k_nak{i},3);
            %         sum_term = 0;
            %         for k_nak = 0:1:m_factor(i)-1
            %             sum_term = sum_term + round(exp(gammaln(2.*k_nak+1)-gammaln(k_nak+1)-gammaln(k_nak+1))).*((1 - mu{i}.^2)./4).^k_nak;
            %         end
            %         PHI{i} = mu{i}.*sum_term;
            PEPth_new{i}(m,:) = 1./size(LUT_C,1).*sum(0.5.*(1 - PHI_new{i,m}),1);
        elseif(strcmp(channel,'AWGN'))
            PEPth_new{i}(m,:) = 1./size(LUT_C,1)*sum(qfunc(sqrt(gamma_new{i,m})),1);
        end

    end
end
for i = 1:N
    BERth(i,:) = mean(PEPth_new{i},[1]);
    dataLUT{i} = [data_LUT{i,:}];
    for l = 1:size(LUT_C,1)
        for k = 1:size(LUT_C,1)
            bitFlip{i}(l,k) = logical(sum(dataLUT{i}(l,:) ~= dataLUT{i}(k,:))>0);
        end
    end
    Sd(i,:) = sqrt(PL(i)./(2*SNR));
    d{i} = Ed(bitFlip{i});
    gamma{i} = (d{i}/2).^2./Sd(i,:).^2;
    if(strcmp(channel,'Ray'))
        PHI{i} = 0.5.*(1 - sqrt(gamma{i}./(2 + gamma{i})));
        PEPth(i,:) = 1./size(LUT_C,1).*sum(PHI{i},1);
    elseif(strcmp(channel,'Nak'))
        mu{i} = sqrt(gamma{i}./(2.*m_factor(i) + gamma{i}));
        k_nak{i} = [0:1:m_factor(i)-1].';
        k_nak{i} = reshape(k_nak{i},[1,1,length(k_nak{i})]);
        PHI{i} = mu{i}.*sum(round(exp(gammaln(2.*k_nak{i}+1)-gammaln(k_nak{i}+1)-gammaln(k_nak{i}+1))).*((1 - mu{i}.^2)./4).^k_nak{i},3);
%         sum_term = 0;
%         for k_nak = 0:1:m_factor(i)-1
%             sum_term = sum_term + round(exp(gammaln(2.*k_nak+1)-gammaln(k_nak+1)-gammaln(k_nak+1))).*((1 - mu{i}.^2)./4).^k_nak;
%         end
%         PHI{i} = mu{i}.*sum_term;
        PEPth(i,:) = 1./size(LUT_C,1).*sum(0.5.*(1 - PHI{i}),1);
    elseif(strcmp(channel,'AWGN'))
        PEPth(i,:) = 1./size(LUT_C,1)*sum(qfunc(sqrt(gamma{i})),1);
    end
    
end
SqERth = PEPth;