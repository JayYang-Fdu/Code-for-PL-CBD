warning off
clear
clc
rng(23)
%% parameters setting
frameCfg = nrFrameCfg_load();
idxMCS = 25; %MCS 0-27

numRx = 64;
numTx = 64;
numTxPort = 1;

ModType = 256;
Qm = log2(ModType);
if Qm==2
    qam = 'QPSK';
elseif Qm==4
    qam = '16QAM';
elseif Qm==6
    qam = '64QAM';
elseif Qm==8
    qam = '256QAM';
end

MonteCalo = 1000;


Qm = frameCfg.mcsTable2(idxMCS+1,1);
Rate = frameCfg.mcsTable2(idxMCS+1,2);

numRB = 32;
numSubCarPerRB = 12;
numDataSymPerSlot = 13;
listTB = numRB;
G = numDataSymPerSlot*listTB*numSubCarPerRB*Qm*numTxPort;  % Generate frame length

v_approach=[
   1,... %SVD-MMSE
   0,... %UCD-VBLAST
   0,... %Soft-UCD-DP
   0,... %CBD
   0,... %PL-CBD
   0,... %SD
    ];
MethodName = strvcat( ...
    'SVD-MMSE', ...
    'UCD-VBLAST', ...
    'UCD-DP', ...
    'CBD',...
    'PL-CBD',...
    'SD'...
    );
blockArray = [2,1];
snrdB = 26:34;12:18;6:12;
BER = zeros(length(snrdB),length(v_approach));
BERBlocksoft = zeros(length(snrdB),length(blockArray));

%%
for mm = 1:MonteCalo
    if mod(mm,100)==0
        fprintf([ '\n', 'MonteCalo = %d ', datestr(now), '\n'], mm);
        % mm
    end
    % Generate channel
    H = (randn(numRx,numTx)+1j*randn(numRx,numTx))/sqrt(2); 
    [U,S,V] = svd(H);
    [Q,B,P] = bidiag(H);
    [TBS,frameCfg] = determine_TBsize(frameCfg,listTB,Rate,Qm,numTxPort,0);
    if v_approach(4)==1
        dpFlag = 1;
        [symbol,origBit] = nrTbSym_gen(frameCfg,TBS,Qm,G,numTxPort,dpFlag);
    else
        dpFlag = 0;
        [symbol,origBit] = nrTbSym_gen(frameCfg,TBS,Qm,G,numTxPort,dpFlag);
    end
    
    xk = reshape(symbol, numTx,[]);
    Gaussian = (randn(size(xk)) + 1j*randn(size(xk)))/sqrt(2);
    for ii = 1:length(snrdB)
        SigmaN = db2pow(-snrdB(ii))*numTx;
        noise = Gaussian*sqrt(SigmaN);
        %%
        if (v_approach(1)==1)
            Method = 'SVD-MMSE';
            Heff_SVD = H;
            yk = Heff_SVD*xk + noise;
            Wmmse=inv(Heff_SVD'*Heff_SVD+eye(numRx)*SigmaN)*Heff_SVD';
            mse = real(diag(inv(Heff_SVD'*Heff_SVD/SigmaN+eye(numRx))));
            snrOut = (1./mse-1);

            ErrorRate = Equalizer(frameCfg,idxMCS,TBS,Wmmse,Heff_SVD,SigmaN, ...
                Method,origBit,Qm,yk,snrOut);
            BER(ii,1) = ErrorRate + BER(ii,1);
        end
        %%
        if (v_approach(2)==1)
            Method = 'UCD-VBLAST';
            [P_ucd,W_ucd,snrOut] = ucd(U,S,V,numTx,SigmaN,0);
            Heff_UCD = H*P_ucd;
            yk = Heff_UCD*xk + noise;
            ErrorRate = Equalizer(frameCfg,idxMCS,TBS,W_ucd,Heff_UCD,SigmaN, ...
                Method,origBit,Qm,yk,snrOut);
            BER(ii,2) = ErrorRate + BER(ii,2);
        end
        %%
        if (v_approach(3)==1)
            Method = 'UCD-DP';
            %% GMD
            % [W_ucd_DP,R,P] = gmd(U,S,V);
            % [dpCodData,Heff_DP,M,d] = dp_transmit(H,P,W_ucd_DP,xk,Qm,numTx);
            % snrOut = R(1,1)/SigmaN;
            % yk = H*P*dpCodData+noise;

            %% UCD

            [P_ucd_DP,W_ucd_DP,snrOut] = ucd(U,S,V,numTx,SigmaN,0);
            [dpCodData,Heff_DP,M,d] = dp_transmit(H,P_ucd_DP,W_ucd_DP,xk,Qm,numTx);
            yk = H*P_ucd_DP*dpCodData + noise;
            ErrorRate = Equalizer(frameCfg,idxMCS,TBS,W_ucd_DP,Heff_DP,SigmaN, ...
                Method,origBit,Qm,yk,snrOut);
            BER(ii,4) = ErrorRate + BER(ii,4);
        end
        %%
        if (v_approach(4)==1)
            Method = 'CBD';
            Heff_CBD = H*P;
            yk = Heff_CBD*xk + noise;
            ErrorRate = Equalizer(frameCfg,idxMCS,TBS,Q,B,SigmaN, ...
                Method,origBit,Qm,yk,0);
            BER(ii,5) = ErrorRate + BER(ii,5);
        end
        %%
        if (v_approach(5)==1)
            Method = 'PL-CBD';
            ErrorRatePlCBD = zeros(blockArray);
            for bb = 1:length(blockArray)
                [Ql,Bl,Pl] =PL_CBD(H, blockArray(bb));
                yk = H*Pl*xk + noise;
                ErrorRatePlCBD(bb) = Equalizer(frameCfg,idxMCS,TBS,Ql,Bl,SigmaN, ...
                    Method,origBit,Qm,yk,0);
                BERBlocksoft(ii,bb) = ErrorRatePlCBD(bb)+BERBlocksoft(ii,bb);
            end
        end
        
        %%
        if (v_approach(6)==1)
            Method = 'SD';
            Heff_CBD = H;
            yk = Heff_CBD*xk + noise;
            ErrorRate = Equalizer(frameCfg,idxMCS,TBS,H,0,SigmaN, ...
                Method,origBit,Qm,yk,0);
            BER(ii,7) = ErrorRate + BER(ii,7);
        end
    end
end
BER = BER/MonteCalo;
BERBlocksoft = BERBlocksoft/MonteCalo;
%% plot
vLegend = [];
fig1 = figure;
color = ['r','m','c','b','y','k','w','g'];
line=['-','--',':','-.'];
marker =['o','+','*','s','d','^','v','x','.'];
set(fig1, 'WindowStyle', 'docked');
PlotType = strvcat('-ro', '-b+', '-cx', '-kv','-gs', ...
    '--ro', '--b+', '--cx', '--kv','--gs', ...
    '-.ro','-.b+', '-.cx', '-.kv','-.gs');
for aa=1:length(v_approach)
    if sum(BER(:,aa))<=0
        continue
    end
    vLegend = strvcat(vLegend,  MethodName(aa,:));
    semilogy(snrdB, BER(:,aa), PlotType(aa,:),'LineWidth',1.5,'MarkerSize',6);
    hold on;
end
xlabel('SNR [dB]');
ylabel('BER');
title(strcat(num2str(numTx),'×',num2str(numRx),', ',qam,', MCS=',num2str(idxMCS)))
grid on;
legend(vLegend,'Location','southwest')
ylim([1e-6 1])
