function ERGAS_ind = ERGAS_Index(x)

%%    Description 

%     This code provides the fuison results of PANchromatic (PAN) and MultiSpectral (MS) images
%     based on Particle Swarm Optimization (PSO) algorithm proposed in the following reference: 

%     [1] A. Azarang and H. Ghassemian, "An adaptive multispectral image fusion
%     using particle swarm optimization," in Proc. Iranian Conf. Elec. Eng.
%     (ICEE), May 2017, pp. 1708-1712.

%     [2] S. Rahmani, M. Strait, D. Merkurjev, M. Moeller, and T. Wittman,
%     "An adaptive IHS Pan-sharpening method," IEEE Geosci. Remote Sens.
%     Lett., vol. 7, no. 4, pp. 746-750, Oct. 2010.

%% Load images and pre-processing steps
addpath QuickBird_Data
load PAN;           %% loading the MS image
load  MS;           %% Loading the PAN image

%% Make the PAN and MS data ready for the processing

MSWV_db  = double(MS);
PANWV_db = double(PAN);
MS_ORG   = double(MS);

%% Resizing, Upsampling the MS data to the size of PAN

MSWV_US  = imresize(MSWV_db,  1/4, 'bicubic');
MSWV_US  = imresize(MSWV_US,  4,   'bicubic');
MS_OBJ = MSWV_US;
PANWV_DS = imresize(PANWV_db, 1/4, 'bicubic');
PANWV_US = imresize(PANWV_DS, 4, 'bicubic');

%% Seperating the spectral bands

R   = MSWV_US(:,:,1); 
G   = MSWV_US(:,:,2); 
B   = MSWV_US(:,:,3); 
NIR = MSWV_US(:,:,4); 

%% Data Normialization

for i=1:size(MSWV_US,3)
    bandCoeffs(i)      =  max(max(MSWV_US(:,:,i)));
    MSWV_US(:,:,i)     =  MSWV_US(:,:,i)/bandCoeffs(i);
end

P = PANWV_DS;
panCoeff = max(max(P));
P = P/panCoeff;

lamda = 10^-9;
eps   = 10^-10;

% PAN weight in edge detector

F_P   = expEdge(P, lamda,eps);

% MS weight in edge detector
Red_W = max(max(R));      %% Red component
R     = R/Red_W;
F_R   = expEdge(R, lamda, eps);

Green_W = max(max(G));    %% Green component
G       = G/Green_W;
F_G     = expEdge(G, lamda,eps);

Blue_W  = max(max(B));    %% Blue component
B       = B/Blue_W;
F_B     = expEdge(B, lamda,eps);

NIR_W=max(max(NIR));      %% NIR component
NIR=NIR/NIR_W;
F_NIR = expEdge(NIR, lamda,eps);


W = impGradDes(MSWV_US,P);

I = W(1)*MSWV_US(:,:,1) + W(2)*MSWV_US(:,:,2) + W(3)*MSWV_US(:,:,3) + W(4)*MSWV_US(:,:,4);

P = (P-mean(P(:)))*std(I(:))/std(P(:)) + mean(I(:));

F_PSO(:,:,1) = MSWV_US(:,:,1) + (x(1)*F_P+(1-x(1))*F_R).*(P-I);
F_PSO(:,:,2) = MSWV_US(:,:,2) + (x(2)*F_P+(1-x(2))*F_G).*(P-I);
F_PSO(:,:,3) = MSWV_US(:,:,3) + (x(3)*F_P+(1-x(3))*F_B).*(P-I);
F_PSO(:,:,4) = MSWV_US(:,:,4) + (x(4)*F_P+(1-x(4))*F_NIR).*(P-I);


    for i=1:size(MS_ORG, 3)
        
    F_PSO(:,:,i)=F_PSO(:,:,i)*bandCoeffs(i);
    
    end
    
    ERGAS_ind = ERGAS(MS_OBJ, F_PSO, 4);
    
end
