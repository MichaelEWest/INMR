%% Create HYP folders
clear
clc

dirP1A = 'D:\ML_Local_Data\Images\HyPOINT_VDPMasks_unorganized\Phase_1A\';
dirP1B = 'D:\ML_Local_Data\Images\HyPOINT_VDPMasks_unorganized\Phase_1B\';
dirP1As = 'D:\ML_Local_Data\Images\HyPOINT_VDPMasks_unorganized\Phase_1A_Scored_DRandJP\';
dirOLS = 'D:\ML_Local_Data\Images\HyPOINT_VDPMasks_unorganized\One_Last_Set\';
dirADDL = 'D:\ML_Local_Data\Images\HyPOINT_VDPMasks_unorganized\Additional_Images\';
dirRDM = 'D:\ML_Local_Data\Images\HyPOINT_VDPMasks_unorganized\Reader_Defect_Masks\';

dir_HYP = 'D:\ML_Local_Data\Images\HyPOINT\';

P1Alist = dir(dirP1A);
P1Blist = dir(dirP1B);
P1Aslist = dir(dirP1As);
OLSlist = dir(dirOLS);
ADDLlist = dir(dirADDL);

%Phase 1A
for nn =3:1:numel(P1Alist)
    PT = P1Alist(nn).name;
    if or(contains(PT,'A_Vent_1'),contains(PT,'A_Vent_2'))
        PT_trim = PT(1:18);
    end
    if or(contains(PT,'A_Vent.'),contains(PT,'A_Vent_m'))
        PT_trim = PT(1:16);
    end
    if contains(PT,'B_Vent')
        PT_trim = PT(1:16);
    end
    if exist(strcat(dir_HYP,PT_trim),'dir') ==0
        mkdir(strcat(dir_HYP,PT_trim))
    end
end

%Phase 1B
for nn =3:1:numel(P1Blist)
    PT = P1Blist(nn).name;
    if or(contains(PT,'A_Vent_1'),contains(PT,'A_Vent_2'))
        PT_trim = PT(1:18);
    end
    if or(contains(PT,'A_Vent.'),contains(PT,'A_Vent_m'))
        PT_trim = PT(1:16);
    end
    if contains(PT,'B_Vent')
        PT_trim = PT(1:16);
    end
    if exist(strcat(dir_HYP,PT_trim),'dir') ==0
        mkdir(strcat(dir_HYP,PT_trim))
    end
end

%Phase 1A (Scored)
for nn =3:1:numel(P1Aslist)
    PT = P1Aslist(nn).name;
    if or(contains(PT,'A_Vent_1'),contains(PT,'A_Vent_2'))
        PT_trim = PT(1:18);
    end
    if or(contains(PT,'A_Vent.'),contains(PT,'A_Vent_m'))
        PT_trim = PT(1:16);
    end
    if contains(PT,'B_Vent')
        PT_trim = PT(1:16);
    end
    if exist(strcat(dir_HYP,PT_trim),'dir') ==0
        mkdir(strcat(dir_HYP,PT_trim))
    end
end

%"One last set"
for nn =3:1:numel(OLSlist)
    PT = OLSlist(nn).name;
    if or(contains(PT,'A_Vent_1'),contains(PT,'A_Vent_2'))
        PT_trim = PT(1:18);
    end
    if or(contains(PT,'A_Vent.'),contains(PT,'A_Vent_m'))
        PT_trim = PT(1:16);
    end
    if contains(PT,'B_Vent')
        PT_trim = PT(1:16);
    end
    if exist(strcat(dir_HYP,PT_trim),'dir') ==0
        mkdir(strcat(dir_HYP,PT_trim))
    end
end

%"Additional Images"
for nn =3:1:numel(ADDLlist)

    PT = ADDLlist(nn).name;
    if or(contains(PT,'A_Vent_1'),contains(PT,'A_Vent_2'))
        PT_trim = PT(1:18);
    end
    if or(contains(PT,'A_Vent.'),contains(PT,'A_Vent_m'))
        PT_trim = PT(1:16);
    end
    if contains(PT,'B_Vent')
        PT_trim = PT(1:16);
    end
    if exist(strcat(dir_HYP,PT_trim),'dir') ==0
        mkdir(strcat(dir_HYP,PT_trim))
    end
end
%%

Hyplist =dir(dir_HYP);
RDMlist = dir(dirRDM);
%Find and move vent files

kk=0;
for nn = 3:1:numel(Hyplist)
    PT = Hyplist(nn).name;
    if strcmp(PT,'HYP0040111B_Vent') ==1

    end
    %Grab vent and mask
    for tt =3:1:numel(P1Alist)
        if contains(P1Alist(tt).name,Hyplist(nn).name)
            copyfile(strcat(dirP1A,P1Alist(tt).name),strcat(dir_HYP,PT,'/',P1Alist(tt).name))
        end
    end
    for tt =3:1:numel(P1Blist)
        if contains(P1Blist(tt).name,Hyplist(nn).name)
            %Grab vent and mask
            copyfile(strcat(dirP1B,P1Blist(tt).name),strcat(dir_HYP,PT,'/',P1Blist(tt).name))
        end
    end
    for tt =3:1:numel(P1Aslist)
        if contains(P1Aslist(tt).name,Hyplist(nn).name)
            %Grab vent and mask
            copyfile(strcat(dirP1As,P1Aslist(tt).name),strcat(dir_HYP,PT,'/',P1Aslist(tt).name))
        end
    end
    for tt =3:1:numel(OLSlist)
        if contains(OLSlist(tt).name,Hyplist(nn).name)
            %Grab vent and mask
            copyfile(strcat(dirOLS,OLSlist(tt).name),strcat(dir_HYP,PT,'/',OLSlist(tt).name))
        end
    end
    for tt =3:1:numel(ADDLlist)
        if contains(ADDLlist(tt).name,Hyplist(nn).name)
            %Grab vent and mask
            copyfile(strcat(dirADDL,ADDLlist(tt).name),strcat(dir_HYP,PT,'/',ADDLlist(tt).name))
        end
    end

    %Grab reader defects
    for tt = 3:1:numel(RDMlist)
        if contains(RDMlist(tt).name,Hyplist(nn).name)
            copyfile(strcat(dirRDM,RDMlist(tt).name),strcat(dir_HYP,PT,'/',RDMlist(tt).name))
        end
    end
    if numel(dir(strcat(dir_HYP,Hyplist(nn).name))) ~=6
        kk = kk+1;
        MissingDataList{kk} = Hyplist(nn).name;
    end
end




