%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
%%    Description 

%     This code provides the fuison results of PANchromatic (PAN) and MultiSpectral (MS) images
%     based on Particle Swarm Optimization (PSO) algorithm proposed in the following reference: 

%     [1] A. Azarang and H. Ghassemian, "An adaptive multispectral image fusion
%     using particle swarm optimization," in Proc. Iranian Conf. Elec. Eng.
%     (ICEE), May 2017, pp. 1708-1712.

%     [2] S. Rahmani, M. Strait, D. Merkurjev, M. Moeller, and T. Wittman,
%     "An adaptive IHS Pan-sharpening method," IEEE Geosci. Remote Sens.
%     Lett., vol. 7, no. 4, pp. 746-750, Oct. 2010.

%     [3] Y. Leung, J. Liu, and J. Zhang, "An improved adaptive intensity-huesaturation
%     method for the fusion of remote sensing images," IEEE Geosci.
%     Remote Sens. Lett., vol. 11, no. 5, pp. 985-989, May 2014.
%%    The steps invovled in the proposed method is summerized as: 

%     1) Pre-processing the datasets,
%     2) Estimating the primitive detail map, 
%     3) Extracting the primitive PAN and MS edge detectors,
%     4) Optimizing the weights of edge detectors using the PSO algorithm
%     and ERGAS loss function to be applied on primitive detail map.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Loading the images and pre-processing steps
addpath QuickBird_Data  %% Dataset path
load PAN;               %% loading the MS image
load  MS;               %% Loading the PAN image

%% Make the PAN and MS data ready for the processing

MSWV_db  = double(MS);
PANWV_db = double(PAN);
MS_ORG   = double(MS);

%% Resizing, Upsampling the MS data to the size of PAN

MSWV_US  = imresize(MSWV_db,  1/4, 'bicubic');
MSWV_US  = imresize(MSWV_US,  4,   'bicubic');
MSWV_DG  = MSWV_US;
PANWV_DS = imresize(PANWV_db, 1/4, 'bicubic');
PANWV_US = imresize(PANWV_DS, 4,   'bicubic');

%% Seperating the spectral bands

R   = MSWV_US(:,:,1); 
G   = MSWV_US(:,:,2); 
B   = MSWV_US(:,:,3); 
NIR = MSWV_US(:,:,4); 

%% Data Normialization

for i=1:size(MSWV_US,3)
    bandCoeffs(i)      = max(max(MSWV_US(:,:,i)));
    MSWV_US(:,:,i)     = MSWV_US(:,:,i)/bandCoeffs(i);
end

P = PANWV_DS;
panCoeff = max(max(P));
P = P/panCoeff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem Definition
CostFunction=@(x) ERGAS_Index(x);          % Cost Function

nVar = 4;                                  % Number of Decision Variables

VarSize = [1 nVar];                        % Size of Decision Variables Matrix

VarMin = 0;                                % Lower Bound of Variables
VarMax = 1;                                % Upper Bound of Variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSO Parameters

MaxIt = 25;         % Maximum Number of Iterations

nPop = 10;          % Population Size (Swarm Size)

% PSO Parameters
% w=0.9;            % Inertia Weight
% wdamp=0.9958;     % Inertia Weight Damping Ratio
% c1=2.0;           % Personal Learning Coefficient
% c2=2.0;           % Global Learning Coefficient

%     If you would like to use Constriction Coefficients for PSO,
%     uncomment the following block and comment the above set of parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constriction Coefficients
phi1  = 2.05;
phi2  = 2.05;
phi   = phi1+phi2;
chi   = 2/(phi-2+sqrt(phi^2-4*phi));
w     = chi;          % Inertia Weight
wdamp = 1;            % Inertia Weight Damping Ratio
c1    = chi*phi1;     % Personal Learning Coefficient
c2    = chi*phi2;     % Global Learning Coefficient

% Velocity Limits

VelMax = 0.1*(VarMax-VarMin);
VelMin = -VelMax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization

empty_particle.Position = [];
empty_particle.Cost     = [];
empty_particle.Velocity = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

particle = repmat(empty_particle,nPop,1);

GlobalBest.Cost = 0;

for i = 1:nPop
    
    % Initialize Position
    particle(i).Position = unifrnd(VarMin,VarMax,VarSize);
    
    % Initialize Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Evaluation
    particle(i).Cost = CostFunction(particle(i).Position);
    
    % Update Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost > GlobalBest.Cost
        
    GlobalBest = particle(i).Best;
        
    end
    
end

BestCost = zeros(MaxIt,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PSO Main Loop

for it = 1:MaxIt
%     figure
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        %% Plot Nodes 
%         plot(particle(i).Position,'*r','MarkerSize',8)
%         axis tight
%         hold on
%         grid on
%         Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;
            
         % Update Global Best
         if particle(i).Best.Cost < GlobalBest.Cost
                
         GlobalBest = particle(i).Best;
                
         end
            
         end
        
    end
    
    BestCost(it) = GlobalBest.Cost;
    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
    w = w*wdamp;
    
end

BestSol = GlobalBest;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pansharpening Framework

% Edge detecotrs 

lamda = 10^-9;
eps   = 10^-10;

% PAN weight in edge detector

F_P   = expEdge(P, lamda,eps);

% MS weight in edge detector
Red_W = max(max(R));      % Red component
R     = R/Red_W;
F_R   = expEdge(R, lamda, eps);


Green_W = max(max(G));    % Green component
G       = G/Green_W;
F_G     = expEdge(G, lamda,eps);


Blue_W  = max(max(B));    % Blue component
B       = B/Blue_W;
F_B     = expEdge(B, lamda,eps);


NIR_W=max(max(NIR));      % NIR component
NIR=NIR/NIR_W;
F_NIR = expEdge(NIR, lamda,eps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracting primitive detail map

W = impGradDes(MSWV_US,P);   % Optimal weights for spectral bands

I     = W(1)*MSWV_US(:,:,1) + W(2)*MSWV_US(:,:,2) + W(3)*MSWV_US(:,:,3) + W(4)*MSWV_US(:,:,4); 
P     = (P-mean(P(:)))*std(I(:))/std(P(:)) + mean(I(:));  % Histogram matching

W_R   = 4*R./(R + B + G + NIR);
W_B   = 4*B./(R + B + G + NIR);
W_G   = 4*G./(R + B + G + NIR);
W_NIR = 4*NIR./(R + B + G + NIR);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fusion result

F_PSO(:,:,1) = MSWV_US(:,:,1) + W_R.*(BestSol.Position(1)*F_P + (1-BestSol.Position(1))*F_R).*(P-I);
F_PSO(:,:,2) = MSWV_US(:,:,2) + W_G.*(BestSol.Position(2)*F_P + (1-BestSol.Position(2))*F_G).*(P-I);
F_PSO(:,:,3) = MSWV_US(:,:,3) + W_B.*(BestSol.Position(3)*F_P + (1-BestSol.Position(3))*F_B).*(P-I);
F_PSO(:,:,4) = MSWV_US(:,:,4) + W_NIR.*(BestSol.Position(4)*F_P + (1-BestSol.Position(4))*F_NIR).*(P-I);


for i = 1:size(MS_ORG, 3)
    F_PSO(:,:,i)       = F_PSO(:,:,i)*bandCoeffs(i);
end
figure, imshow(uint8(MSWV_DG(:,:,1:3)),'border','tight');
figure, imshow(uint8(PANWV_DS),'border','tight');
figure, imshow(uint8(F_PSO(:,:,1:3)),'border','tight');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective Assesment of Fusion Results

addpath Objective_Evaluation

Methods = {'PSO'};
ERGAS   = ERGAS(MS_ORG,F_PSO, 4);
SAM     = SAM(MS_ORG,F_PSO);
RASE    = RASE(MS_ORG,F_PSO);
RMSE    = RMSE(MS_ORG,F_PSO);
UIQI    = uqi(MS_ORG,F_PSO);
CC      = CC(MS_ORG,F_PSO);
T       = table(ERGAS, SAM, RASE, RMSE, UIQI, CC, 'RowNames', Methods)

% End of Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%