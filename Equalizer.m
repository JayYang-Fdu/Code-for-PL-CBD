function [BER] = Equalizer(frameCfg,MCS,TBS,Q,R,SigmaN,Method,origBit,Qm,yk,snrOut)
% INPUT:
%  frameCfg: frame config
%  TBS: TB Size
%  Q: left unitary matrix
%  R: effective channel
%  P: precoding matrix
%  SigmaN: varaience of noise
%  MethodName: choose a algorithm to simulate
%  origBit: uncoded message
%  Qm: bits in a symbol
%  yReceive: the data received
%  snrOut: snr of perlayer

% OUTPUT:
%  BER:  the BER of equalizer
Rate = frameCfg.mcsTable2(MCS+1,2);
numTxPort = 1;
if Qm==2
    QAM = 'QPSK';
    ModType = 4;
    d = 2/sqrt(2);
elseif Qm==4
    QAM = '16QAM';
    ModType = 16;
    d = 2/sqrt(10);
elseif Qm==6
    QAM = '64QAM';
    ModType = 64;
    d = 2/sqrt(42);
elseif Qm==8
    QAM = '256QAM';
    ModType = 256;
    d = 2/sqrt(170);
end

if TBS > 3824
    CRC = '24A';
    L = 24;
else
    CRC = '16';
    L = 16;
end

[Nr,N] = size(yk);  % N transmit times
graphType = frameCfg.graphType;
rvId = frameCfg.rvId;

if Method == "UCD-VBLAST"
    
    zk = zeros(Nr,N);
    for k =Nr:-1:1
        if k<Nr
            yk = yk-R(:,k+1)*decData_ucd(k+1,:);
        end
        zk(k,:) = Q(:,k)'*yk;
        CodedBitsHat_UCD = nrSymbolDemodulate(zk(k,:).',QAM,1)*snrOut;
        CodedBitsHat_UCD(CodedBitsHat_UCD>30)=30;
        CodedBitsHat_UCD(CodedBitsHat_UCD<-30)=-30;
        decBitLLR (k,:) = CodedBitsHat_UCD.';
        CodedBitsHat_UCD(CodedBitsHat_UCD>=0)=0;
        CodedBitsHat_UCD(CodedBitsHat_UCD<0)=1;

        ModData_ucd = nrSymbolModulate(CodedBitsHat_UCD,QAM);
        decData_ucd(k,:) = ModData_ucd.';
    end
    % decBitsLLR = reshape(decBitLLR.',[],1);
    for n = 1:size(decBitLLR,2)/Qm
        decBitsLLR((n-1)*Qm*Nr+(1:Qm*Nr)) = reshape(decBitLLR(:,(n-1)*Qm+(1:Qm)).',1,[]);
    end
    decodeBitsTb = LDPCdecodeForsimNr(frameCfg,decBitsLLR,TBS,L,Rate,rvId,QAM,numTxPort,graphType,CRC);
    BER = mean(origBit.'~=decodeBitsTb);
elseif Method == "SVD-MMSE"
    snrOut = repmat(snrOut.',[Qm,1]);
    snrOut = snrOut(:);
    yHat_MMSE = Q*yk./diag(Q*R);
    symbolHat_MMSE = reshape(yHat_MMSE,Nr*N,1);
    LLR_MMSE = nrSymbolDemodulate(symbolHat_MMSE,QAM,1);
    CodedBitsHat_MMSE = reshape(LLR_MMSE, [Nr*Qm,N]).*snrOut;
    % CodedBitsHat_MMSE = reshape(LLR_MMSE, [Nr*Qm,N]);
    decBitsLLR = reshape(CodedBitsHat_MMSE,1,[]);
    decodeBitsTb = LDPCdecodeForsimNr(frameCfg,decBitsLLR,TBS,L,Rate,rvId,QAM,numTxPort,graphType,CRC);
    BER = mean(origBit.'~=decodeBitsTb);
elseif Method == "UCD-DP"
    M = 2^(Qm/2);
    yHat_DP = Q'*yk;
    y1 = diag(diag(R))\yHat_DP;
    dpRxData = real(y1)-floor((real(y1)+M*d/2)/(M*d))*M*d + ...
           1j*(imag(y1)-floor((imag(y1)+M*d/2)/(M*d))*M*d);
    symbolHat_UCDDP = reshape(dpRxData,Nr*N,1);
    decBitsLLR = softDemod_dp_v2(Qm,symbolHat_UCDDP).*snrOut;
    decodeBitsTb = LDPCdecodeForsimNr(frameCfg,decBitsLLR,TBS,L,Rate,rvId,QAM,numTxPort,graphType,CRC);
    BER = mean(origBit.'~=decodeBitsTb);
elseif Method == "CBD"
    yHat_CBD = Q'*yk;
    LeCkP = zeros(Nr*N*Qm, 1);
    LCkY = MIMO_optimal_detection(yHat_CBD,LeCkP, SigmaN,R, ModType, Nr); % MAP
    LCkY(LCkY>30)=30;
    LCkY(LCkY<-30)=-30;
    decBitsLLR = -LCkY.';
    decodeBitsTb = LDPCdecodeForsimNr(frameCfg,decBitsLLR,TBS,L,Rate,rvId,QAM,numTxPort,graphType,CRC);
    BER = mean(origBit.'~=decodeBitsTb);
elseif Method == "PL-CBD"
    yHat_CBD = Q'*yk;
    LeCkP = zeros(Nr*N*Qm, 1);
    LCkY = MIMO_optimal_detection(yHat_CBD,LeCkP, SigmaN,R, ModType, Nr); % MAP
    LCkY(LCkY>30)=30;
    LCkY(LCkY<-30)=-30;
    decBitsLLR = -LCkY.';
    decodeBitsTb = LDPCdecodeForsimNr(frameCfg,decBitsLLR,TBS,L,Rate,rvId,QAM,numTxPort,graphType,CRC);
    BER = mean(origBit.'~=decodeBitsTb);
elseif Method == "SD"
    yHat_SD = yk;
    LCkY = nrSphereDecoder(Q,yHat_SD,QAM,'soft')/SigmaN; % SD
    
    LCkY(LCkY>30)=30;
    LCkY(LCkY<-30)=-30;
    decBitsLLR = LCkY.';
    decodeBitsTb = LDPCdecodeForsimNr(frameCfg,decBitsLLR,TBS,L,Rate,rvId,QAM,numTxPort,graphType,CRC);
    BER = mean(origBit.'~=decodeBitsTb);
    
end

% descrambleOut = reshape(decBitsLLR,[],1);
% [C,~,~,~,~] = NR_cbSize(frameCfg,TBS+L);
% deRateMatchOut = nrRateRecoverLDPC(descrambleOut,TBS,Rate/1024,rvId,QAM,numTxPort,C);
% % LDPC decode
% decodeBitsCb = nrLDPCDecode(deRateMatchOut,graphType,20,'DecisionType','soft');
% % LDPC code block desegment
% [decodeBitsTbCrc,~] = nrCodeBlockDesegmentLDPC(decodeBitsCb,graphType,TBS+L);
% % TB CRC decode
% [decodeBitsTb,~] = nrCRCDecode(decodeBitsTbCrc,CRC);
% % calculate the number of error bits
% BER = mean(origBit.'~=decodeBitsTb);



end