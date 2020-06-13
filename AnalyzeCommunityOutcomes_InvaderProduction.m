%% Analyze the frequency of invasion
clear

infile = {
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt0_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt1_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt2_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt5_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt20_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt50_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    'CRC_Nif20_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt100_at50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns1000_rndseed6130.mat'
    };

for btcnt = 1:8
    load(infile{btcnt})
    CommCntC(btcnt) = sum(NE0DC(1,:)>1);
    InvCommAD = (InvEffAD1(1:Nif-1,NE0DC(1,:)>1)>=0.95);
    NoInvCommAD = (InvEffADC(1:Nif-1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
    numInvCommAD = sum(InvCommAD,1);
    WeakCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==0);
    DependentCommInv = (numInvCommAD>0).*(numInvCommAD<Nif-1).*(InvCommAD(1,:)==1);
    UnableCommInv = (numInvCommAD==0);
    StrongCommInv = (numInvCommAD==Nif-1);
    InvStatCommAD(:,btcnt) = 1/CommCntC(btcnt)*[sum(StrongCommInv); sum(DependentCommInv); sum(WeakCommInv); sum(UnableCommInv)];
    nif = 8; %Invader at 0.3%
    AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    CommStatCommADC(:,btcnt) = 1/CommCntC(btcnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];
end


%% Plots
btrng = 1/50*[5 10 20 50 100 250 500];
% categories: 'Strong', 'Dependent', 'Weak', 'Failed'
% figure
% bar(qp,InvStatCommAD','stack')
% xlabel('Chance of mediator production by the invader')
% ylabel('Frequency')

% categories: 'displace', 'augment', 'perturb', 'resist'
mycolor=[0.3 0.4 0.5; 0.2 0.6 1; 0.5 0.1 0.7; 1 0.6 0.2];
figure
hh = area(btrng,CommStatCommADC(:,2:8)');
set(hh(1), 'FaceColor',mycolor(1,:))
set(hh(2), 'FaceColor',mycolor(2,:))
set(hh(3), 'FaceColor',mycolor(3,:))
set(hh(4), 'FaceColor',mycolor(4,:))
% hold on
% bar(0.02,CommStatCommADC(:,1),'stack')
set(gca,'XScale','log')
xlim([1e-1 10])
ylim([0 1])
xlabel('Normalized consumption rate')
ylabel('Frequency')
