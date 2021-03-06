function [Ne, Cmp, Ne0, Cmp0] = WellmixedInteraction_DpMM_ExMTCC(nRound,r0Vector,cellRatioArray,intMat,nInitialCell,kSatVector,A,B,kSatLevel,ExtTh,DilTh,tauf,dtau)

%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
% 5: populations growth rates < 0.9 of max rate in the past 20 generations thrown out from coexistence
% rndseed = 1389;
% rand('twister',rndseed)

% nCellType = 15; % # of cell types
% nMediator = 6; % # of mediators
% nRound = 15; % number of rounds of propagation
% r0 = 0.08+0.04*rand(Nc,1); % population reproduction rates, per hour
% nInitialCell = 1e4; % total initial cells
% kSatLevel = 1e7; % interaction strength saturation level of each population
% ExtTh = 0.1; % population extinction threshold
% DilTh = 1e10; % coculture dilution threshold
tau0 = 0;
% tauf = 250; % in hours
% dtau = 0.01; % in hours, cell growth update and uptake timescale
% at = 0.1; % avg. consumption values (fmole per cell); alpha_ij: population i, resource j
% bt = 1; % avg. production rates (fmole per cell per hour); beta_ij: population i, resource j
% mp = 3; % average number of production links per population
% mc = 2; % average number of consumption links per population

% intMat : % matrix of interaction coefficients

[nCellType, nMediator] = size(intMat);

%% Parameters
% R = zeros(Nc,Nm);
% rndc = rand(Nc,Nm);
% R(rndc <= mc/Nm) = 1;
% P = zeros(Nc,Nm);
% rndp = rand(Nc,Nm);
% P(rndp <= mp/Nm) = 1;

% interaction matrix
% alpha = at*(0.5+rand(Nc,Nm)); % consumption rates
% beta = bt*(0.5+rand(Nc,Nm)); % mediator release rates
% A = (R.*alpha)';
% B = (P.*beta)';

%% Initial state 
% cellRatioArray = 1 / nCellType * ones(1,nCellType) % cell distrbution
cMedVector = zeros(nMediator,1); % concentrations of interaction mediators

%% Cell-growth time-course
tauScaleArray = tau0:dtau:tauf;
nTauScale = length(tauScaleArray);

tc = zeros(1,nRound*nTauScale);
Xc = zeros(nCellType,nRound*nTauScale);
Cc = zeros(nMediator,nRound*nTauScale);

cct = 0;
nCellVector = nInitialCell * cellRatioArray'; % initial number of each cell type

for iRound = 1 : nRound
    cMedVector = nInitialCell / sum(nCellVector) * cMedVector;
    nCellVector = nInitialCell * cellRatioArray'; % initial number of each cell type
    
    tau0 = 0; % in hours
    tau = tau0;
    
    nCellOnEachScale = zeros(nCellType,nTauScale);
    cMedOnEachScale = zeros(nMediator,nTauScale);
    rzm = zeros(nCellType,nTauScale);
    count = 0;
    while (tau<=tauf-dtau) && (sum(nCellOnEachScale(:,max(count,1)))<DilTh)
        
        count = count+1;
        tau = tauScaleArray(count);
            
        rIntPerCellVector = r0Vector + ((intMat < 0) .* intMat) * (cMedVector ./ kSatVector) + ((intMat >= 0).* intMat) * (cMedVector ./ (cMedVector + kSatVector));
                    
        nCellVector = nCellVector + dtau * (rIntPerCellVector .* nCellVector);
        nCellVector(nCellVector < ExtTh) = 0;
        
        Ce = cMedVector * ones(1, nCellType); % nM * nC matrix, each column is cMed
        AMM = A.*Ce./(Ce+kSatLevel);
        cMedVector = cMedVector + dtau*(B*nCellVector - AMM*nCellVector);
        cMedVector(cMedVector<0) = 0;
        
        nCellOnEachScale(:,count) = nCellVector;
        cMedOnEachScale(:,count) = cMedVector;
        rzm(:,count) = rIntPerCellVector;
        
    end
    cellRatioArray = 1/sum(nCellOnEachScale(:,count))*nCellOnEachScale(:,count)';
    
    if cct ==0
        tc(cct+1:cct+count) = tauScaleArray(1:count);
    else
        tc(cct+1:cct+count) = tc(cct) + tauScaleArray(1:count);
    end
    Xc(:,cct+1:cct+count) = nCellOnEachScale(:,1:count);
    Cc(:,cct+1:cct+count) = cMedOnEachScale(:,1:count);
    cct = cct+count;

end
indx = 1:nCellType;
Ne0 = indx(nCellOnEachScale(:,count)>0);
Cmp0 = cellRatioArray(Ne0);
% get Cmp as percentage each cell type contributes to the total community
if sum(Cmp0) > 0
    Cmp_sum = zeros(1,size(Cmp0,2));
    Cmp_sum(1,:) = sum(Cmp0);

    Cmp0 = Cmp0./Cmp_sum;
end

nGen = log(DilTh/nInitialCell)/log(2);
r = (nCellOnEachScale(:,count)./nCellOnEachScale(:,1)).^(20/nGen);
stp = (r > abs(0.9*max(r)));
Ne = indx(stp);
Cmp = cellRatioArray(Ne);
% get Cmp as percentage each cell type contributes to the total community
if sum(Cmp) > 0
    Cmp_sum = zeros(1,size(Cmp,2));
    Cmp_sum(1,:) = sum(Cmp);

    Cmp = Cmp./Cmp_sum;
end

return;
