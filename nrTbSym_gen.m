function [symMod,origBit] = nrTbSym_gen(frameCfg,TBS,Qm,G,Nl,dpFlag)
%  Authors: Fengjie Li（18210720070@fudan.edu.cn）;
%  copyright - CSRL@fudan 2020/06/24
%  inputs  - frameCfg: frame configuration
%            TBS:transport block size
%            Qm:modulation order 
%            numRbPerTb:number of resource block per transport block
%            G:the total number of coded bits available for transmission of the transport block
%            Nl:the number of transmission layers that the transport block is mapped onto
%            tbIdx:transport block(TB) index           
%  outputs - symMod：symbol after modulation
%            origBit:original bits per transport block

% generate original bits
origBit =  randi([0 1],1,TBS);

% add TB CRC
graphType = frameCfg.graphType;
rvId = frameCfg.rvId;
if TBS > 3824
    CRC = '24A';
else
    CRC = '16';
end
tbCRCout = nrCRCEncode(origBit.',CRC);

% add CB CRC
cbCRCout = nrCodeBlockSegmentLDPC(tbCRCout,graphType);

% LDPC encode
codeBits = nrLDPCEncode(cbCRCout,graphType);
if Qm==2
    qam = 'QPSK';
elseif Qm==4   
    qam = '16QAM';
elseif Qm==6
    qam = '64QAM';
elseif Qm==8
    qam = '256QAM';
end

% LDPC rate matching and interleaving
ldpcRmOut = nrRateMatchLDPC(codeBits,G,rvId,qam,Nl);

% scrambling
% scrambleOut = NR_scrambling(frameCfg,nodeCfg,ldpcRmOut.',tbIdx,numRbPerTb,Qm);
scrambleOut = ldpcRmOut.';

% modulation
if dpFlag==1
    symMod = NR_modulation(scrambleOut.',Qm);
elseif dpFlag==0
    symMod = nrSymbolModulate(scrambleOut.',qam);
    symMod = symMod.';
end
end

