function [frameCfg] = nrFrameCfg_load()
%  Authors: Fengjie Li（18210720070@fudan.edu.cn）;
%  copyright - CSRL@fudan 2020/07/03
%% NR frame configuration
frameCfg = struct(...
    'tbCRCType',                '24a',          ...
    'cbCRCType',                '24b',          ...
    'numSubCarPerRB',           12,             ...
    'Nf',                       4096,           ...
    'numSlotPerSubfrm' ,        2,              ...
    'dmrsSymPos0',              2,              ...
    'numDmrsSymPerSlot',        1,              ...
    'ptrsON',                   0,              ...
    'ptrsTimeDensity',          1,              ...   time density for PT-RS(clause 5.1.6.3 in TS38.214)
    'ptrsFreqDensity',          2,              ...   frequency density for PT-RS(clause 5.1.6.3 in TS38.214)
    'ptrsSymTimePos',           [0:1,3:13],     ...
    'dataSymTimePos',           [0:1,3:13],     ...
    'numPtrsSymPerSlot',        13,             ...
    'kREoffset00',              [0,2,1,3],      ...   (clause Table 7.4.1.2.2-1 in TS38.214)
    'numOFDMsymPerSlot',        14,             ...
    'numDataSymPerSlot',        13,             ...
    'lenOCCt',                  2,              ...
    'numOCCt',                  2,              ...
    'OCCt',                     [],             ...
    'rvId',                     0,              ...
    'q',                        0,              ...  % q = 0 for single-codeword transmission
    'normalCpLen0' ,            352,            ...
    'normalCpLen1' ,            288,            ...
    'fsHz',                     100*1e6  ,      ...  % 100*1e6
    'fcHz',                     3.5e9,          ...  % 3.5GHz
    'stbMaxRb',                 32,             ...
    'graphType',                [],             ...
    'TBsizeTable',              [],             ...
    'table2MCS',                []);

frameCfg.mcsTable2 = [2,120,  0.2344;2,193,0.3770;2,308,0.6016;2,449,0.8770;2,602,1.1758;...
                      4,378,  1.4766;4,434,1.6953;4,490,1.9141;4,553,2.1602;4,616,2.4063;4,658,2.5703;...
                      6,466,  2.7305;6,517,3.0293;6,567,3.3223;6,616,3.6094;6,666,3.9023;6,719,4.2129;6,772,  4.5234;6,822,4.8164;6,873,5.1152;...
                      8,682.5,5.3320;8,711,5.5547;8,754,5.8906;8,797,6.2266;8,841,6.5703;8,885,6.9141;8,916.5,7.1602;8,948,7.4063];
frameCfg.tbsTable = [24,32,40,48,56,64,72,80,88,96,104,112,120,128,136,144,152,160,168,176,184,192,208,224,240,256,272,288,304,320,...
                     336,352,368,384,408,432,456,480,504,528,552,576,608,640,672,704,736,768,808,848,888,928,984,1032,1064,1128,1160,1192,1224,1256,...
                     1288,1320,1352,1416,1480,1544,1608,1672,1736,1800,1864,1928,2024,2088,2152,2216,2280,2408,2472,2536,2600,2664,2728,2792,2856,2976,3104,3240,3368,3496,...
                     3624,3752,3824];
frameCfg.ldpcLiftingSizeTable = cell(1,8);
frameCfg.ldpcLiftingSizeTable{1} = [2,4,8,16,32,64,128,256];
frameCfg.ldpcLiftingSizeTable{2} = [3,6,12,24,48,96,192,384];
frameCfg.ldpcLiftingSizeTable{3} = [5,10,20,40,80,160,320];
frameCfg.ldpcLiftingSizeTable{4} = [7,14,28,56,112,224];
frameCfg.ldpcLiftingSizeTable{5} = [9,18,36,72,144,288];
frameCfg.ldpcLiftingSizeTable{6} = [11,22,44,88,176,352];
frameCfg.ldpcLiftingSizeTable{7} = [13,26,52,104,208];
frameCfg.ldpcLiftingSizeTable{8} = [15,30,60,120,240];

frameCfg.OCCt = [+1,+1;
                +1,-1];
end
