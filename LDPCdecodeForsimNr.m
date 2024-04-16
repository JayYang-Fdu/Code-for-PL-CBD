function [decodeBitsTb] = LDPCdecodeForsimNr(frameCfg,decBitsLLR,TBS,L,Rate,rvId,QAM,numTxPort,graphType,CRC)

descrambleOut = reshape(decBitsLLR,[],1);
[C,~,~,~,~] = NR_cbSize(frameCfg,TBS+L);
deRateMatchOut = nrRateRecoverLDPC(descrambleOut,TBS,Rate/1024,rvId,QAM,numTxPort,C);
% LDPC decode
decodeBitsCb = nrLDPCDecode(deRateMatchOut,graphType,20);
% LDPC code block desegment
[decodeBitsTbCrc,~] = nrCodeBlockDesegmentLDPC(decodeBitsCb,graphType,TBS+L);
% TB CRC decode
[decodeBitsTb,~] = nrCRCDecode(decodeBitsTbCrc,CRC);
% calculate the number of error bits

end