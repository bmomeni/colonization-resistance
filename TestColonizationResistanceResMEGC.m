%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% ExMT2: corrected the error in ExMT, now rIntMat only includes links in R

clear
rndseed0 = 2418;
rng(rndseed0,'twister');

nSample = 1; % # of samples being screened
rndseed = round(100000*rand(1,nSample));

Nif = 1; % number of invader fractions tested
InvFracRng = 0.003; %logspace(-4,-0.004364,Nif); % fraction of invader
nGen = 200;
nInitialCell = 1e4; % total initial cells
dilTh = 1e7; % coculture dilution threshold
GenPerRound = log(dilTh/nInitialCell)/log(2);
nRound = round(nGen/GenPerRound); % number of rounds of propagation

nCellType = 20; % # of cell types in the initial pool
nMediator = 10; % # of mediators
kSatLevel = 1e4; % interaction strength saturation level of each population
extTh = 0.1; % population extinction threshold
r0m = 0.1; % mean basal growth rate, 1/hr
r0d = 0.04; % basal growth rate deviation, 1/hr
r0Inv = 0.25; % basal growth rate of invader, 1/hr
ri0 = 0.2; % maximum interaction strength, 1/hr
ri0I = ri0; % maximum interaction strength of invader, 1/hr
fp = 0.1; % fraction of interactions that are positive
fpI = 0.1; % fraction of interactions that are positive
tau0 = 0; % in hours
tauf = 80; % in hours
dtau = 3.6; % in hours, cell growth update and uptake timescale
at = 0.5; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
bt = 0.1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
qp = 0.3; % probability of production link per population
qc = 0.3; % probability of influence link per population
qpi = 0.3; % probability of production link per population
qci = 0.3; % probability of influence link per population

CR = 1e7; % Resource concentration, fmole/ml
Ares0 = 1e3;
Ares = (0.5+0.5*rand(1,nCellType))*Ares0; % Resource consumption rate, fmole per cell
Kres0 = 1e8;
Kres = (0.5+0.5*rand(1,nCellType))*Kres0; % Resource consumption affinity, fmole/ml

r0T = zeros(nCellType,nSample); %matrix of zeros: 4 rows because 4 species and 1000 columns because 1000 samples being screened
SiT = zeros(nMediator,nCellType,nSample);
AT = zeros(nMediator,nCellType,nSample);
BT = zeros(nMediator,nCellType,nSample);
CmpADT = zeros(nCellType,Nif,nSample);
CmpADCT = zeros(nCellType,Nif,nSample);
CmpADIT = zeros(nCellType+1,Nif,nSample);
CmpADCIT = zeros(nCellType+1,Nif,nSample);
rintAT = zeros(nMediator,nCellType,nSample);
V0DT = zeros(3,nCellType,nSample);
VDT = zeros(3,nCellType,nSample);

NE0D = zeros(3,nSample);
NE0DC = zeros(3,nSample);

CmpFADIF = zeros(nCellType+1,nSample);
CmpFBDIF = zeros(nCellType+1,nSample);
CmpFEDIF = zeros(nCellType+1,nSample);

InvEffAD1 = zeros(Nif,nSample);
InvEffBD1 = zeros(Nif,nSample);
InvEffED1 = zeros(Nif,nSample);

InvEffADC = zeros(Nif,nSample);
InvEffBDC = zeros(Nif,nSample);
InvEffEDC = zeros(Nif,nSample);

for ns = 1:nSample
    disp(ns)
    tic
    
    rng(rndseed(ns),'twister');
    r0 = (r0m-r0d/2) + r0d*rand(nCellType,1); % population reproduction rates, per hour
    kSatMat = kSatLevel * (0.5 + rand(nMediator, nCellType)); % population levels for influence saturation
    
    %% Parameters
    % Network configuration
    % NetworkConfig_Balanced(Nc,Nm,q): link between Nc and Nm present with
    % a probability q
    R = NetworkConfig_Binomial(nCellType,nMediator,qc);
    P = NetworkConfig_Binomial(nCellType,nMediator,qp);
    
    % interaction matrix
    alpha = at * (0.5+rand(nCellType,nMediator)); % consumption rates
    beta = bt * (0.5+rand(nCellType,nMediator)); % mediator release rates
    A = (R.*alpha)';
    B = (P.*beta)';
    
    rIntMatA = R .* DistInteractionStrengthMT_PA(nCellType, nMediator, ri0); % matrix of interaction coefficients, 50/50
    rIntMatB = R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, fp); % matrix of interaction coefficients, more negative
    rIntMatE = R .* DistInteractionStrengthMT_PB(nCellType, nMediator, ri0, 1-fp); % matrix of interaction coefficients, more positive
    
    Cmp = 1 / nCellType * ones(1,nCellType); % cell distribution; population ratios
    
    %% Simulating dynamics, Dp, depletable
    [NeADC, CmpADC, NeAD, CmpAD] = WellmixedInteraction_ResME_Plot(nGen,r0,Cmp,rIntMatA,nInitialCell,kSatMat,A,B,Ares,CR,Kres,extTh,dilTh,tauf,dtau,11);
    
    V0AD = zeros(1,nCellType);
    V0AD(NeAD) = 1;
    
    V0ADC = zeros(1,nCellType);
    V0ADC(NeADC) = 1;
    
    NE0D(:,ns) = sum(V0AD,2);
    NE0DC(:,ns) = sum(V0ADC,2);
    
    CmpADT(NeAD,ns) = CmpAD;
    
    CmpADCT(NeADC,ns) = CmpADC;

    Cmp0AD = zeros(1,nCellType);
    Cmp0AD(NeAD) = CmpAD;
    
    Cmp0ADC = zeros(1,nCellType);
    Cmp0ADC(NeADC) = CmpADC;
    
    r0T(:,ns) = r0;
    SiT(:,:,ns) = kSatMat;
    AT(:,:,ns) = A;
    BT(:,:,ns) = B;
    rintAT(:,:,ns) = rIntMatA';
    
    % Properties of the invader
    r0I = [r0; r0Inv];
    RcI = NetworkConfig_Binomial(1,nMediator,qci);
    RpI = NetworkConfig_Binomial(1,nMediator,qpi);
    kSatMatI = [kSatMat, kSatLevel*(0.5+rand(nMediator,1))];
    rIntMatAI = [rIntMatA; RcI.*DistInteractionStrengthMT_PB(1,nMediator,ri0I,fpI)];
    %     rIntMatBI = [rIntMatB; RcI.*DistInteractionStrengthMT_PB(1,nMediator,ri0I,fpI)];
    %     rIntMatEI = [rIntMatE; RcI.*DistInteractionStrengthMT_PB(1,nMediator,ri0I,fpI)];
    AI = [A'; at*RcI.*(0.5+rand(1,nMediator))]';
    BI = [B'; bt*RpI.*(0.5+rand(1,nMediator))]';
    AresI = [Ares, (0.5+0.5*rand(1))*Ares0]; % Resource cnsumption rate, fmole per cell
    KresI = [Kres, (0.5+0.5*rand(1))*Kres0]; % Resource cnsumption affinity, fmole/ml
    
    for nif = 1:Nif
        InvFrac = InvFracRng(nif);
        
        Cmp0ADI = [(1-InvFrac)*Cmp0ADC, InvFrac];
        
        %% Simulating dynamics after introduction of invader, Dp, depletable
        [NeADCIF, CmpADCIF, NeADIF, CmpADIF] = WellmixedInteraction_ResME_Plot(nGen,r0I,Cmp0ADI,rIntMatAI,nInitialCell,kSatMatI,AI,BI,AresI,CR,KresI,extTh,dilTh,tauf,dtau,12);
        
        CmpADIT(NeADIF,nif,ns) = CmpADIF;
        InvEffAD1(nif,ns) = shiftdim(CmpADIT(nCellType+1,nif,ns),1)/InvFrac;
        
        CmpADCIT(NeADCIF,nif,ns) = CmpADCIF;
        InvEffADC(nif,ns) = shiftdim(CmpADCIT(nCellType+1,nif,ns),1)/InvFrac;
        
    end
    
    toc
end

% rith = 0.001;
% rInt = rIntMatAI';
% Nm = nMediator;
% Nc = nCellType+1;
% rintEn = zeros(Nm,Nc);
% AtEn = zeros(Nm,Nc);
% BtEn = zeros(Nm,Nc);
% Ncrng = 1:Nc;
% NeD = [NeADCIF 21];
% rintEn(:,NeD) = rInt(:,NeD);
% rintEn = rintEn.*((sum(BI(:,NeD),2)>0)*ones(1,Nc));
% rintEn(abs(rintEn)<rith) = 0; % remove weak interactions
% AtEn(:,NeD) = AI(:,NeD);
% AtEn = AtEn.*((sum(BI(:,NeD),2)>0)*ones(1,Nc));
% BtEn(:,NeD) = BI(:,NeD);
% BtEn = BtEn.*((sum(AI(:,NeD),2)>0)*ones(1,Nc));
% 
% rth = 0.001;
% PlotInteractionNetworkExMT4_StrengthBasalFitness(r0I,rIntMatAI,AI,BI,rth,1)
% PlotInteractionNetworkExMT4_StrengthBasalFitness(r0I,rintEn',AtEn,BtEn,rth,2)


CommCntC = sum(NE0DC(1,:)>1);
InvCommAD = (InvEffADC(1,NE0DC(1,:)>1)>=0.95);
NoInvCommAD = (InvEffADC(1,NE0DC(1,:)>1)<1); % finding cases that invader frequency increased
nif = 1; %Invader at 0.3%
DisplaceCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
AugmentCommInv = InvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
PerturbCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)>shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
ResistCommInv = NoInvCommAD(nif,:).*(sum(CmpADCT(:,NE0DC(1,:)>1)>0,1)<=shiftdim(sum(CmpADCIT(1:nCellType,nif,NE0DC(1,:)>1)>0,1),1));
CommStatCommADC = 1/CommCntC*[sum(DisplaceCommInv); sum(AugmentCommInv); sum(PerturbCommInv); sum(ResistCommInv)];
disp(CommStatCommADC)