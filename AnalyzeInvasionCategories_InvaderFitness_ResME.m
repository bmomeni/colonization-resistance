%% Analyze the frequency of invasion
clear

infile = 'CR_Nif20_Dp_ExMT4C_ABE_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI0_r0m10_r0Inv0_bt10_at50_Nc20_Nm10_qp30_qc30_qpi0_qci0_Nr20_Ns1000_rndseed6130.mat';
load(infile)
ricnt = 1;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv5_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 2;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv10_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 3;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv12_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 4;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 5;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];


infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv20_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 6;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv25_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 7;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv30_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 8;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

infile = 'CR_ResMEG_Nif20_ABE_CR1e+08_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv35_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_nGen200_Ns1000_rndseed7905.mat';
load(infile)
ricnt = 9;
InvAD = (InvEffAD1(1:Nif-1,:)>=0.95);
numInvAD = sum(InvAD,1);
WeakInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==0);
DependentInv = (numInvAD>0).*(numInvAD<Nif-1).*(InvAD(1,:)==1);
UnableInv = (numInvAD==0);
StrongInv = (numInvAD==Nif-1);
InvStatAD(:,ricnt) = 1/nSample*[sum(StrongInv); sum(DependentInv); sum(WeakInv); sum(UnableInv)];

CommCnt(ricnt) = sum(NE0D(1,:)>1);
InvCommAD = (InvEffAD1(1:Nif-1,NE0D(1,:)>1)>=0.95);
numInvCommAD = sum(InvCommAD,1);
WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
UnableCommInv = (numInvCommAD==0);
StrongCommInv = (numInvCommAD==Nif-1);
InvStatCommAD(:,ricnt) = 1/CommCnt(ricnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];

%% Plots
r0Irng = 1/0.1*[0 0.05 0.1 0.12 0.15 0.20 0.25 0.3 0.35];
% figure
% bar(r0Irng,InvStatAD','stack')
% xlabel('Normalized invader growth rate')
% ylabel('Frequency')

% categories: 'Strong', 'Dependent', 'Weak', 'Failed'
figure
bar(r0Irng,InvStatCommAD','stack')
xlabel('Normalized invader growth rate')
ylabel('Frequency')

figure
plot(r0Irng,InvStatCommAD(1,:),'k')
hold on
plot(r0Irng,sum(InvStatCommAD(1:2,:)),'k')
plot(r0Irng,sum(InvStatCommAD(1:3,:)),'k')
plot(r0Irng,sum(InvStatCommAD(1:4,:)),'k')
xlim([0 3.5])
ylim([0 1])
xlabel('Normalized invader growth rate')
ylabel('Frequency')