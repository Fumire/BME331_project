clc
clear all
close all

%% ----- Parameter setting ----- %%

Mat = 256; % Matrix size (256 x 256)

TR = [50 100 200 300 500 700 1000 2000 3000 5000]; % Repetetion time (TR) of RARE-VTR sequence (T1), [ms]

TE_msme = [8 16 24 32 40 48 56 64 72 80 88 96 104 112 120 128 136 144 152 160 168 176 ...
    184 192 200 208 216 224 232 240 248 256 264 272 280 288 296 304 312 320 328 ...
    336 344 352 360 368 376 384]; % Echo time (TE) of MSME sequence (T2), [ms]

TE_mge = [2.7 6.1 9.5 12.9 16.3 19.7 23.1 26.5 29.9 33.3 36.7 40.1 43.5 46.9 50.3]; % TE of MGE sequence (T2*), [ms]

load Fid_rat_6week_rarevtr.mat % Load FID data (RARE-VTR)

%% ----- Fill k-space using the given FID data ----- %%


%%%%%%%%%% Fill out this part %%%%%%%%%%


%% ----- Apply 2D Inverse Fourier transform  ----- %%


%%%%%%%%%% Fill out this part %%%%%%%%%%


%% ----- T1 fitting (RARE-VTR)   ----- %%


%%%%%%%%%% Fill out "¡Ú" part %%%%%%%%%%


MR_signals = ¡Ú; % MR signal array to fitting (RARE-VTR)

opt = fitoptions('Method','NonlinearLeastSquares');
opt.StartPoint=[0 1800]; 
opt.Lower=[0 0]; 
opt.Upper=[inf inf];

f=fittype('S0*(¡Ú¡Ú)','independent','¡Ú','coefficients',{'S0','¡Ú'},'options',opt); 
% ¡Ú¡Ú = equation part 
[myfit,goodness] = fit(TR',MR_signals',f);


%% ----- T2 fitting (MSME)   ----- %%


%%%%%%%%%% Fill out "¡Ú" part %%%%%%%%%%


MR_signals = ¡Ú; % MR signal array to fitting (MSME)

opt = fitoptions('Method','NonlinearLeastSquares');
opt.StartPoint=[0 50 0]; 
opt.Lower=[0 0 0]; 
opt.Upper=[inf inf inf];

f=fittype('S0*(¡Ú¡Ú)+C','independent','¡Ú','coefficients',{'S0','¡Ú','C'},'options',opt); % C = noise. 
% ¡Ú¡Ú = equation part 
[myfit,goodness] = fit(TE_msme',MR_signals',f);


%% ----- T2* fitting (MGE)   ----- %%

%%%%%%%%%% Fill out "¡Ú" part %%%%%%%%%%

MR_signals = ¡Ú; % MR signal array to fitting (MGE)

opt = fitoptions('Method','NonlinearLeastSquares');
opt.StartPoint=[0 30]; 
opt.Lower=[0 0]; 
opt.Upper=[inf inf];

f=fittype('S0*(¡Ú¡Ú)','independent','¡Ú','coefficients',{'S0','¡Ú'},'options',opt); 
% ¡Ú¡Ú = equation part 
[myfit,goodness] = fit(TE_mge',MR_signals',f);

