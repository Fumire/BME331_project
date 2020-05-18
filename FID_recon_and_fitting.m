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

mge_6week_data = load("FID_data/Fid_rat_4month_mge").Fid_rat_4month_mge;
tmp = size(mge_6week_data);
mge_6week_data = complex(mge_6week_data(1:2:tmp(1), :), mge_6week_data(2:2:tmp(1), :));
mge_6week_data = reshape(mge_6week_data, Mat, Mat * size_mge(2));

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

mge_6week_image_data = zeros(size_mge(2), Mat, Mat);
for i = 1:size_mge(2)
    mge_6week_image_data(i, :, :) = mge_6week_data(:, i:size_mge(2):(size_mge(2) * Mat));
end

parfor i = 1:size_mge(2)
    tmp = ifftshift(ifftn(mge_6week_image_data(i, :, :)));
    tmp = squeeze(sqrt(sum(abs(tmp).^2, 4)));
    tmp = transpose(tmp);
    
    image = imagesc(tmp);
    title(strcat("TE = ", num2str(TE_mge(i)), " ms"));
    axis image;
    colormap gray;
    
    saveas(image, strcat("images/mge_", int2str(i), ".png"));
    mge_6week_image_data(i, :, :) = tmp;
end

%% ----- T1 fitting (RARE-VTR)   ----- %%

MR_signals = zeros(1, size_TR(2));
parfor i = 1:size_TR(2)
    MR_signals(1, i) = mean(vtr_6week_image_data(i, :, :), 'all');
end

opt = fitoptions('Method', 'NonlinearLeastSquares');
opt.StartPoint = [0 1800];
opt.Lower = [0 0];
opt.Upper = [inf inf];

f = fittype('M0 * (1-exp(-t / T1))', 'independent', {'t'}, 'dependent', {'Mz'}, 'coefficients', {'M0', 'T1'}, 'options', opt);
[myfit, goodness] = fit(TR', MR_signals', f);

plot(myfit, TR, MR_signals);
title(strcat("T1 = ", num2str(myfit.T1), " ms; R^2 = ", num2str(goodness.rsquare)));

saveas(gcf, "images/T1_fit.png");
close all;

%% ----- T2 fitting (MSME)   ----- %%

MR_signals = zeros(1, size_msme(2));
parfor i = 1:size_msme(2)
    MR_signals(1, i) = mean(msme_6week_image_data(i, :, :), 'all');
end

opt = fitoptions('Method', 'NonlinearLeastSquares');
opt.StartPoint = [0 50];
opt.Lower = [0 0];
opt.Upper = [inf inf];

f = fittype('M0 * exp(-t/T2)', 'independent', 't', 'dependent', 'Mxy', 'coefficients', {'M0', 'T2'}, 'options', opt);
[myfit, goodness] = fit(TE_msme', MR_signals', f);

plot(myfit, TE_msme, MR_signals);
title(strcat("T2 = ", num2str(myfit.T2), " ms; R^2 = ", num2str(goodness.rsquare)));

saveas(gcf, "images/T2_fit.png");
close all;

%% ----- T2* fitting (MGE)   ----- %%

MR_signals = zeros(1, size_mge(2));
parfor i = 1:size_mge(2)
    MR_signals(1, i) = mean(mge_6week_image_data(i, :, :), 'all');
end

opt = fitoptions('Method', 'NonlinearLeastSquares');
opt.StartPoint = [0 30];
opt.Lower = [0 0];
opt.Upper = [inf inf];

f = fittype('M0 * exp(-t / T2)', 'independent', 't', 'dependent', 'Mxy', 'coefficients', {'M0', 'T2'}, 'options', opt);
[myfit, goodness] = fit(TE_mge', MR_signals', f);

plot(myfit, TE_mge, MR_signals);
title(strcat("T2 = ", num2str(myfit.T2), " ms; R^2 = ", num2str(goodness.rsquare)));

saveas(gcf, "images/T2star_fit.png");
close all

%% ----- T1 Mapping ----- %%

vtr_map_data = zeros(Mat, Mat);
tmp = size_TR(2);

parfor i = 1:Mat
    tmp2 = zeros(1, Mat);
    for j = 1:Mat
        MR_signals = zeros(1, tmp);
        for x = 1:tmp
            MR_signals(1, x) = vtr_6week_image_data(x, i, j);
        end

        opt = fitoptions('Method', 'NonlinearLeastSquares');
        opt.StartPoint = [0 1800];
        opt.Lower = [0 0];
        opt.Upper = [inf inf];

        f = fittype('M0 * (1-exp(-t / T1))', 'independent', {'t'}, 'dependent', {'Mz'}, 'coefficients', {'M0', 'T1'}, 'options', opt);
        [myfit, goodness] = fit(TR', MR_signals', f);
        tmp2(1, j) = min([myfit.T1 2000]);
    end
    vtr_map_data(i, :) = tmp2;
end

image = imagesc(vtr_map_data);
title("T1 map");
axis image;
colorbar;

saveas(image, "images/T1_map.png");

%% ----- T2 Mapping ----- %%

msme_map_data = zeros(Mat, Mat);
tmp = size_msme(2);

parfor i = 1:Mat
    tmp2 = zeros(1, Mat);
    for j = 1:Mat
        MR_signals = zeros(1, tmp);
        for x = 1:tmp
            MR_signals(1, x) = msme_6week_image_data(x, i, j);
        end
        
        opt = fitoptions('Method', 'NonlinearLeastSquares');
        opt.StartPoint = [0 50];
        opt.Lower = [0 0];
        opt.Upper = [inf inf];

        f = fittype('M0 * exp(-t/T2)', 'independent', {'t'}, 'dependent', {'Mxy'}, 'coefficients', {'M0', 'T2'}, 'options', opt);
        [myfit, goodness] = fit(TE_msme', MR_signals', f);
        tmp2(1, j) = min([50 myfit.T2]);
    end
    msme_map_data(i, :) = tmp2;
end

image = imagesc(msme_map_data);
title("T2 map")
axis image;
colorbar;

saveas(image, "images/T2_map.png");

%% ----- T2* Mapping ----- %%

mge_map_data = zeros(Mat, Mat);
tmp = size_mge(2);

parfor i = 1:Mat
    tmp2 = zeros(1, Mat);
    for j = 1:Mat
        MR_signals = zeros(1, tmp);
        for x = 1:tmp
            MR_signals(1, x) = mge_6week_image_data(x, i, j);
        end
        
        opt = fitoptions('Method', 'NonlinearLeastSquares');
        opt.StartPoint = [0 50];
        opt.Lower = [0 0];
        opt.Upper = [inf inf];

        f = fittype('M0 * exp(-t/T2)', 'independent', {'t'}, 'dependent', {'Mxy'}, 'coefficients', {'M0', 'T2'}, 'options', opt);
        [myfit, goodness] = fit(TE_mge', MR_signals', f);
        tmp2(1, j) = min([50 myfit.T2]);
    end
    mge_map_data(i, :) = tmp2;
end

image = imagesc(mge_map_data);
title("T2* map");
axis image;
colorbar;

saveas(image, "images/T2star_map.png");
