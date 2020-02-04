function [Ne, Cmp, Ne0, Cmp0] = WellmixedInteraction_ResME(nGen,r0Vector,cellRatioArray,intMat,nInitialCell,kSatMat,A,B,Ares,CR,Kres,ExtTh,DilTh,tauf,dtau)
%% ResME: Resource limited (resource supplied in the medium)
% ResA: vector describing the consumption of resource by each species (per cell)
% CR: amount of resource supplied in the media
% Kres: resource dependent logistic growth K

%% ExMT5: kSatMat (different per species) replaces kSatVactor
%% ExMT5: removed kSatLevel, using kSatMat instead, to match the formulation in the manuscript (v4.4)

%% Well-mixed model for growth of interacting species
% UIC: Uniform initial condition
% Ex: explicitly including the mediators
% MT: multi-target mediators
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

cct = 0;
Res = CR;
nCellVector = nInitialCell * cellRatioArray'; % initial number of each cell type
AccGen = 0;
iRound = 0;
while AccGen < nGen
    iRound = iRound + 1;
    Res = nInitialCell/sum(nCellVector)*Res + (1 - nInitialCell/sum(nCellVector))*CR;
    cMedVector = nInitialCell/sum(nCellVector) * cMedVector;
    nCellVector = nInitialCell * cellRatioArray'; % initial number of each cell type
    
    tau0 = 0; % in hours
    tau = tau0;
    
    nCellOnEachScale = zeros(nCellType,nTauScale);
    cMedOnEachScale = zeros(nMediator,nTauScale);
    ResOnEachScale = zeros(1,nTauScale);
    rzm = zeros(nCellType,nTauScale);
    count = 0;
    while (tau<=tauf-dtau) && (sum(nCellOnEachScale(:,max(count,1)))<DilTh)
        
        count = count+1;
        tau = tauScaleArray(count);
        cMedMat = cMedVector*ones(1,nCellType);
     
        rIntPerCellVector = (r0Vector + (((intMat < 0) .* intMat)./ kSatMat') * cMedVector + (((intMat >= 0).* intMat)./(cMedMat + kSatMat)') * cMedVector).*(Res./(Res+Kres'));
        
        dRes = dtau * Ares * (rIntPerCellVector .* nCellVector);
        nCellVector = nCellVector + dtau * (rIntPerCellVector .* nCellVector);
        nCellVector(nCellVector < ExtTh) = 0;
        Res = Res - dRes;
        Res = (Res>0)*Res;
        
        Ce = cMedVector * ones(1, nCellType); % nM * nC matrix, each column is cMed
        AMM = A.*Ce./(Ce+kSatMat);
        cMedVector = cMedVector + dtau*(B-AMM)*((Res./(Res+Kres')).*nCellVector);
        cMedVector(cMedVector<0) = 0;
        
        nCellOnEachScale(:,count) = nCellVector;
        cMedOnEachScale(:,count) = cMedVector;
        rzm(:,count) = rIntPerCellVector;
        ResOnEachScale(count) = Res;
        
    end
    cellRatioArray = 1/sum(nCellOnEachScale(:,count))*nCellOnEachScale(:,count)';
    AccGen = AccGen + log(sum(nCellOnEachScale(:,count))/nInitialCell)/log(2);
    
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

nGen = log(sum(nCellOnEachScale(:,count))/nInitialCell)/log(2);
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
