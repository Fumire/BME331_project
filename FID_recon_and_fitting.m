clc
clear
close all

%% ----- Parameter setting ----- %%

Mat = 256; % Matrix size (256 x 256)

TR = [50 100 200 300 500 700 1000 2000 3000 5000]; % Repetetion time (TR) of RARE-VTR sequence (T1), [ms]
size_TR = size(TR);

TE_msme = [8 16 24 32 40 48 56 64 72 80 88 96 104 112 120 128 136 144 152 160 168 176 ...
    184 192 200 208 216 224 232 240 248 256 264 272 280 288 296 304 312 320 328 ...
    336 344 352 360 368 376 384]; % Echo time (TE) of MSME sequence (T2), [ms]
size_msme = size(TE_msme);

TE_mge = [2.7 6.1 9.5 12.9 16.3 19.7 23.1 26.5 29.9 33.3 36.7 40.1 43.5 46.9 50.3]; % TE of MGE sequence (T2*), [ms]
size_mge = size(TE_mge);

%% ----- Fill k-space using the given FID data ----- %%

vtr_6week_data = load("FID_data/Fid_rat_6week_rarevtr").Fid_rat_6week_rarevtr;
tmp = size(vtr_6week_data);
vtr_6week_data = complex(vtr_6week_data(1:2:tmp(1), :), vtr_6week_data(2:2:tmp(1), :));
vtr_6week_data = reshape(vtr_6week_data, Mat * Mat, size_TR(2));

msme_6week_data = load("FID_data/Fid_rat_6week_msme").Fid_rat_6week_msme;
tmp = size(msme_6week_data);
msme_6week_data = complex(msme_6week_data(1:2:tmp(1), :), msme_6week_data(2:2:tmp(1), :));
msme_6week_data = reshape(msme_6week_data, Mat, Mat * size_msme(2));

%% ----- Apply 2D Inverse Fourier transform  ----- %%

vtr_6week_image_data = zeros(size_TR(2), Mat, Mat);
parfor i = 1:size_TR(2)
    tmp = transpose(reshape(vtr_6week_data(:, i), Mat, Mat));
    tmp2 = zeros(Mat, Mat);
    
    for j = 1:Mat
        [row, r] = quorem(sym(j), 2);
        
        if r == 0
            row = row * -1
        end
        row = row + Mat / 2 + 1
        
        tmp2(row, :) = tmp(j, :);
    end
    
    tmp2 = ifftshift(ifftn(tmp2));
    tmp2 = squeeze(sqrt(sum(abs(tmp2).^2, 4)));
    
    image = imagesc(tmp2);
    title(strcat("TR = ", int2str(TR(i)), " ms"));
    axis image;
    colormap gray;
    
    saveas(image, strcat("images/vtr_", int2str(i), ".png"));
    
    vtr_6week_image_data(i, :, :) = tmp2;
end

msme_6week_image_data = zeros(size_msme(2), Mat, Mat);
for i = 1:size_msme(2)
    msme_6week_image_data(i, :, :) = msme_6week_data(:, i:size_msme(2):(size_msme(2) * Mat));
end

parfor i = 1:size_msme(2)
    tmp = ifftshift(ifftn(msme_6week_image_data(i, :, :)));
    tmp = squeeze(sqrt(sum(abs(tmp).^2, 4)));
    tmp = transpose(tmp);
    
    image = imagesc(tmp);
    title(strcat("TE = ", int2str(TE_msme(i)), " ms"));
    axis image;
    colormap gray;
    
    saveas(image, strcat("images/msme_", int2str(i), ".png"));
    msme_6week_image_data(i, :, :) = tmp;
end

return;

%% ----- T1 fitting (RARE-VTR)   ----- %%


%%%%%%%%%% Fill out "this" part %%%%%%%%%%


MR_signals = ''; % MR signal array to fitting (RARE-VTR)

opt = fitoptions('Method','NonlinearLeastSquares');
opt.StartPoint=[0 1800]; 
opt.Lower=[0 0]; 
opt.Upper=[inf inf];

f=fittype('S0*','independent','asdf','coefficients',{'S0','asdf'},'options',opt); 
% �ڡ� = equation part 
[myfit,goodness] = fit(TR',MR_signals',f);


%% ----- T2 fitting (MSME)   ----- %%


%%%%%%%%%% Fill out "��" part %%%%%%%%%%


MR_signals = ''; % MR signal array to fitting (MSME)

opt = fitoptions('Method','NonlinearLeastSquares');
opt.StartPoint=[0 50 0]; 
opt.Lower=[0 0 0]; 
opt.Upper=[inf inf inf];

f=fittype('S0*(�ڡ�)+C','independent','��','coefficients',{'S0','��','C'},'options',opt); % C = noise. 
% �ڡ� = equation part 
[myfit,goodness] = fit(TE_msme',MR_signals',f);


%% ----- T2* fitting (MGE)   ----- %%

%%%%%%%%%% Fill out "��" part %%%%%%%%%%

MR_signals = ''; % MR signal array to fitting (MGE)

opt = fitoptions('Method','NonlinearLeastSquares');
opt.StartPoint=[0 30]; 
opt.Lower=[0 0]; 
opt.Upper=[inf inf];

f=fittype('S0*(�ڡ�)','independent','��','coefficients',{'S0','��'},'options',opt); 
% �ڡ� = equation part 
[myfit,goodness] = fit(TE_mge',MR_signals',f);

