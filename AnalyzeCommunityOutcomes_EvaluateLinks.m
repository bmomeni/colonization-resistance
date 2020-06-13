%% Analyze the frequency of invasion
clear

infile = {
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_btI10_atI50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_btI0_atI50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_btI10_atI0_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI50_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI20_r0m10_r0Inv15_bt10_at50_btI0_atI0_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI0_r0m10_r0Inv15_bt10_at50_btI10_atI50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI0_r0m10_r0Inv15_bt10_at50_btI0_atI50_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI0_r0m10_r0Inv15_bt10_at50_btI10_atI0_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    'CRC_if003_Dp_ExMTC_ABE_fp10_fpI10_CSD1e+04_DilTh1e+07_ExtTh0.1_Ksat10000_ri20_riI0_r0m10_r0Inv15_bt10_at50_btI0_atI0_Nc20_Nm10_qp30_qc30_qpi30_qci30_Nr20_Ns10000_rndseed2803.mat'
    };

for casecnt = 1:length(infile)
    load(infile{casecnt})
    CommCntC(casecnt) = sum(NE0DC(1,:)>1);
    InvCommAD = (InvEffAD1(1,NE0DC(1,:)>1)>=0.95);
    NoInvCommAD = (InvEffADC(1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
    nif = 1; % only invader at 0.3% simulated
    AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
    CommStatCommADC(:,casecnt) = 1/CommCntC(casecnt)*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];
end

%% Plots
LinkCases = 1:8;

% categories: 'Strong', 'Dependent', 'Weak', 'Failed'
% figure
% bar(Nifrng, InvStatCommADC','stack')
% xlabel('Fraction of positive influences on invader')
% ylabel('Frequency')
% xlim([0 0.6])
% legend('Strong', 'Dependent', 'Weak', 'Failed')

% categories: 'displace', 'augment', 'perturb', 'resist'
% wdth = logspace(-4,-1,Nif-1);
% figure
% for nf = 1:Nif-1
%     mycolor=[0.3 0.3 0.3; 0.7 0.1 0.7; 0.2 0.6 1; 1 0 0];
%     hh = bar([0 Nifrng(nf)], [CommStatCommADC(:,nf) CommStatCommADC(:,nf)]','stack');
%     set(hh(1), 'FaceColor',mycolor(1,:))
%     set(hh(2), 'FaceColor',mycolor(2,:))
%     set(hh(3), 'FaceColor',mycolor(3,:))
%     set(hh(4), 'FaceColor',mycolor(4,:))
%     hold on
% end
% set(gca,'XScale','log')
% xlabel('Propagule size (relative to resident community)')
% ylabel('Frequency')
% xlim([0 1])
% legend('Displace','Augment','Perturb','Resist')
% 
% figure
% plot(Nifrng,CommStatCommADC(1,:),'k')
% hold on
% plot(Nifrng,sum(CommStatCommADC(1:2,:)),'k')
% plot(Nifrng,sum(CommStatCommADC(1:3,:)),'k')
% plot(Nifrng,sum(CommStatCommADC(1:4,:)),'k')
% set(gca,'XScale','log')
% xlim([1e-3 0.6])
% ylim([0 1])
% xlabel('Propagule size (relative to resident community)')
% ylabel('Frequency')

mycolor=[0.3 0.4 0.5; 0.2 0.6 1; 0.5 0.1 0.7; 1 0.6 0.2];
figure
hh = bar(LinkCases,CommStatCommADC','stack');
set(hh(1), 'FaceColor',mycolor(1,:))
set(hh(2), 'FaceColor',mycolor(2,:))
set(hh(3), 'FaceColor',mycolor(3,:))
set(hh(4), 'FaceColor',mycolor(4,:))
xlim([0 9])
ylim([0 1])
set(gca,'XTick',1:8,'XTickLabel',{'all links','no b','no a','no ab','no r','no br','no ar','no abr'})
ylabel('Frequency')

% figure
% area(Nifrng,CommStatCommADC','FaceColor','flat');
% set(gca,'XScale','log')
% xlim([1e-3 0.6])
% ylim([0 1])
% xlabel('Propagule size (relative to resident community)')
% ylabel('Frequency')
