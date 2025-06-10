% MXB201 Project Part II Code.

%% Initialisation
clear
d = dir('supplied_data/faces/*.pgm');
N = length(d);
I = imread([d(1).folder, '/', d(1).name]);
[rows,cols] = size(I);
M = rows*cols;
A = zeros(M, N);  % big matrix, whose columns are the images

%% Read images as columns of the matrix
for j = 1:N
    I = imread([d(j).folder, '/', d(j).name]);
    A(:,j) = I(:);
end

%% Calculate and visualise mean face
Abar = mean(A, 2);
meanFace = reshape(Abar, [rows cols]);
figure
hold on
title('Mean face');
imshow(mat2gray(meanFace));

%% Calculate mean-centred SVD
Acentered = A - Abar;
[U, Sigma, V] = svd(Acentered, 'econ');

%% Visualise first 20 eigenfaces
figure
layout = tiledlayout('flow');
for j = 1:20
    nexttile
    imshow(mat2gray(reshape(U(:, j), [rows cols])))
    title(j)
end
title(layout, "Eignenfaces")

%% Find good nu value
sigma = diag(Sigma);
figure
tiledlayout('flow')
ax1 = nexttile;
semilogy(ax1, sigma / sigma(1))
ax1.YLim = [0.001, 1];
grid on
xlabel('\nu')
ylabel('\sigma_\nu / \sigma_1')
title('Relative contribution of left singular values')

%% Calculate coordinate vectors
nu = 100;
Up = U(:, 1:nu);
Sp = Sigma(1:nu, 1:nu);
Vtp = V(1:nu, 1:nu)';
Aapprox = Up*Sp*Vtp;

%% Demonstrate rudimentary moustache detector
% Moustache face is eigenface 13
moustaches = [239, 551, 254];
nomoustaches = [439, 363, 388]; % also 179
demofaces = [moustaches, nomoustaches];
figure
layout = tiledlayout('flow');
for j = demofaces
    v = Sigma(1:nu, 1:nu) * V(j, 1:nu)';
    if (abs(v(13)) >= 750)
        word = "has";
    else
        word = "doesn't have";
        % continue; % skip while testing
    end
    nexttile
    imshow(mat2gray(reshape(A(:, j), [rows cols])))
    title(sprintf('Face %d %s\n a moustache.\nScore of %.4f', j, word, v(13)))
end
title(layout, "Moustache detector")