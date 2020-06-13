%% Analyze the frequency of invasion outcomes
clear

infile = 'CRC_ResMEG_Nif20_ABE_CR1e+08_Ares1000_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed9861.mat';
load(infile)
fpcnt = 5;
CommCntC(fpcnt) = sum(NE0DC(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

fpcnt = 1;
CommCntC(fpcnt) = sum(NE0DC(2,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

fpcnt = 9;
CommCntC(fpcnt) = sum(NE0DC(3,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

infile = 'CRC_ResMEG_Nif20_ABE_CR1e+08_Ares1000_fp20_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed9861.mat';
load(infile)
fpcnt = 2;
CommCntC(fpcnt) = sum(NE0DC(2,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

fpcnt = 8;
CommCntC(fpcnt) = sum(NE0DC(3,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

infile = 'CRC_ResMEG_Nif20_ABE_CR1e+08_Ares1000_fp30_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed9861.mat';
load(infile)
fpcnt = 3;
CommCntC(fpcnt) = sum(NE0DC(2,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

fpcnt = 7;
CommCntC(fpcnt) = sum(NE0DC(3,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

infile = 'CRC_ResMEG_Nif20_ABE_CR1e+08_Ares1000_fp40_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed9861.mat';
load(infile)
fpcnt = 4;
CommCntC(fpcnt) = sum(NE0DC(2,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

fpcnt = 6;
CommCntC(fpcnt) = sum(NE0DC(3,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
nif = 8;
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC(:,fpcnt) = 1/CommCntC(fpcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];

%% Plots
fprng = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

mycolor=[0.3 0.4 0.5; 0.2 0.6 1; 0.5 0.1 0.7; 1 0.6 0.2];
figure
hh = area(fprng,CommStatCommADC');
set(hh(1), 'FaceColor',mycolor(1,:))
set(hh(2), 'FaceColor',mycolor(2,:))
set(hh(3), 'FaceColor',mycolor(3,:))
set(hh(4), 'FaceColor',mycolor(4,:))
xlim([0.1 0.9])
ylim([0 1])
xlabel('Fraction of facilitation among residents')
ylabel('Frequency')
