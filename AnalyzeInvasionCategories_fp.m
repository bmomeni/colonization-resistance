%% Analyze the frequency of invasion
clear

infile = 'CR_CommIntxn_Nif20_Dp_ExMT4C_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat';
load(infile)
ricnt = 5;
InvD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

ricnt = 1;
InvD = (InvEffBD1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(2,:)>1);
InvCommD = (InvEffBD1(1:Nif-1,NE0D(2,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

ricnt = 9;
InvD = (InvEffED1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(3,:)>1);
InvCommD = (InvEffED1(1:Nif-1,NE0D(3,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_CommIntxn_Nif20_Dp_ExMT4C_ABE_fp20_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat';
load(infile)
ricnt = 2;
InvD = (InvEffBD1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(2,:)>1);
InvCommD = (InvEffBD1(1:Nif-1,NE0D(2,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

ricnt = 8;
InvD = (InvEffED1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(3,:)>1);
InvCommD = (InvEffED1(1:Nif-1,NE0D(3,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_CommIntxn_Nif20_Dp_ExMT4C_ABE_fp30_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat';
load(infile)
ricnt = 3;
InvD = (InvEffBD1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(2,:)>1);
InvCommD = (InvEffBD1(1:Nif-1,NE0D(2,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

ricnt = 7;
InvD = (InvEffED1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(3,:)>1);
InvCommD = (InvEffED1(1:Nif-1,NE0D(3,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_CommIntxn_Nif20_Dp_ExMT4C_ABE_fp40_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat';
load(infile)
ricnt = 4;
InvD = (InvEffBD1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(2,:)>1);
InvCommD = (InvEffBD1(1:Nif-1,NE0D(2,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

ricnt = 6;
InvD = (InvEffED1(1:Nif-1,:)>=0.95);
numInvD = sum(InvD,1);
WeakInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==0);
DependentInv = (numInvD>0).*(numInvD<Nif-1).*(InvD(1,:)==1);
UnableInv = (numInvD==0);
StrongInv = (numInvD==Nif-1);
InvStatD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(3,:)>1);
InvCommD = (InvEffED1(1:Nif-1,NE0D(3,:)>1)>=0.95);
numInvCommD = sum(InvCommD,1);
WeakCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==0);
DependentCommInv = (numInvCommD>0).*(numInvCommD<Nif-1).*(InvCommD(1,:)==1);
UnableCommInv = (numInvCommD==0);
StrongCommInv = (numInvCommD==Nif-1);
InvStatCommD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

%% Plots
fprng = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
% figure
% bar(fprng,InvStatD',0.3,'stack')
% set(gca,'XTick',[0.1:0.2:1])
% xlim([0 1])
% xlabel('Fraction of facilitation in resident comm.')
% ylabel('Frequency')

figure
bar(fprng,InvStatCommD',0.3,'stack')
set(gca,'XTick',[0.1:0.2:1])
xlim([0 1])
xlabel('Fraction of facilitation in resident comm.')
ylabel('Frequency')

figure
plot(fprng,InvStatCommD(1,:),'k')
hold on
plot(fprng,sum(InvStatCommD(1:2,:)),'k')
plot(fprng,sum(InvStatCommD(1:3,:)),'k')
plot(fprng,sum(InvStatCommD(1:4,:)),'k')
xlim([0.1 0.9])
ylim([0 1])
set(gca,'XTick',0.1:0.2:0.9)
xlabel('Fraction of facilitation in resident comm.')
ylabel('Frequency')