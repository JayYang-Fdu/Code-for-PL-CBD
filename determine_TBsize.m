function [TBS,frameCfg] = determine_TBsize(frameCfg,nPRB,R,Qm,numLayer,numExRE)
%  Authors: Fengjie Li（18210720070@fudan.edu.cn）;
%  copyright - CSRL@fudan 2020/06/24
%  inputs  - frameCfg:frame configuration
%            nPRB:number of physical resource block
%            R:coding rate(×1024)
%            Qm:modulation order
%            numLayer: the number of layer
%            numExRE:number of extra RE except DMRS
%  outputs - TBS:transport block size
%            frameCfg:frame configuration
%  This function calculates the transport block size 
%  reference: subclause 5.1.3.2 (Transport block size determination) of TS38.214
%% common parameters
numSubCarPerRB = frameCfg.numSubCarPerRB;
tbsTable = frameCfg.tbsTable;
%% determine the TBS as specified below
numDataSymPerSlot = frameCfg.numDataSymPerSlot;
temp = numSubCarPerRB*numDataSymPerSlot - numExRE;
numRE = min(156,temp)*nPRB;
Ninfo = numRE*R*Qm*numLayer/1024;
if Ninfo <= 3824
    n = max(3,floor(log2(Ninfo))-6);
    NinfoTemp = max(24,2.^n*floor(Ninfo/2.^n));
    idxTB = find(tbsTable - NinfoTemp<=0);
    TBS = tbsTable(idxTB(end));
elseif Ninfo > 3824
    n = floor(log2(Ninfo-24))-5;
    NinfoTemp = max(3840,2^n*round((Ninfo-24)/2^n));
    if (R/1024) <= 1/4
        C = ceil((NinfoTemp+24)/3816);
        TBS = 8*C*ceil((NinfoTemp+24)/8/C)-24;
    else
        if NinfoTemp > 8424
            C = ceil((NinfoTemp+24)/8424);
            TBS = 8*C*ceil((NinfoTemp+24)/8/C)-24;
        else
            TBS = 8*ceil((NinfoTemp+24)/8)-24;
        end
    end
end
% set LDPC graph type(according to clause 7.2.2 LDPC base graph selection in TS38.212)
if TBS<=292
    frameCfg.graphType = 2;
elseif (TBS<=3824)&&(R/1024<=2/3)
    frameCfg.graphType = 2;
elseif (R/1024)<=1/4
    frameCfg.graphType = 2;
else
    frameCfg.graphType = 1;
end
end

