% MXB201 Project Part I

%% Initialisation
clear
load partI

[X,Y,num_dirs] = size(S);
assert(isequal(size(g), [num_dirs 3]));

% These arrays will be be filled in during the main loop below
MD  = nan(X, Y);    % mean diffusion
FA  = nan(X, Y);    % fractional anistropy
PDD = nan(X, Y, 3); % principal diffusion direction

%% Compute the diffusion tensor for each pixel
for x = 1:X
    for y = 1:Y

        % If not within the mask, skip the pixel
        if ~mask(x, y), continue; end
 
        % Handling bad data
        % If signal value S0 is negative, skip the pixel
        if S0(x, y) < 0, continue; end
        
        % Solving least squares problem

        % Solving for matrix A
        gx = g(:,1);
        gy = g(:,2);
        gz = g(:,3);
        A = zeros(num_dirs, 6);
        A(:,1) = gx.^2;
        A(:,2) = gy.^2;
        A(:,3) = gz.^2;
        A(:,4) = 2 * gx .* gy;
        A(:,5) = 2 * gx .* gz;
        A(:,6) = 2 * gy .* gz;
       
        % Solving for vector B
        B = zeros(num_dirs, 1);
        for i = 1:length(B)
        B(i) = (log((S(x, y, i))./(S0(x, y))))/(-b);
        end

        % Filtering out complex results (handling bad data)
        index = find(imag(B) == 0);
        B = B(index);
        A = A(index, :);

        % Solving for matrix Dbar
        Dbar = A \ B;
        
        % Forming diffusion tensor
        D = [Dbar(1) Dbar(4) Dbar(5); Dbar(4) Dbar(2) Dbar(6); Dbar(5) Dbar(6) Dbar(3)];

        % Finding eigenvalues and eigenvectors
        [U, V] = eig(D);
        [P,ind] = sort(diag(V), "descend");
        VS = V(ind,ind);
        US = U(:,ind);
        lambda1 = VS(1, 1);
        lambda2 = VS(2, 2);
        lambda3 = VS(3, 3);
        v1 = US(:,1);
        v2 = US(:,2);
        v3 = US(:,3);

        % Calculating MD
        MD(x, y) = (lambda1+lambda2+lambda3)/3;
        
        % Calculating FA
        FA(x, y) = (sqrt(3)/sqrt(2)) * (((sqrt((lambda1 - MD(x, y))^2 + (lambda2 - MD(x, y))^2 + (lambda3 - MD(x, y))^2)))/(sqrt(lambda1^2 + lambda2^2 + lambda3^2)));

        % Calculating PDD
        xp = [1 0 0];
        yp = [0 1 0];
        zp = [0 0 1];

        % Finding the angle v1 makes with the x axis
        alpha = acos((dot(v1, xp))./(norm(v1).*norm(xp)));

        % Finding the angle v1 makes with the y axis
        beta = acos((dot(v1, yp))./(norm(v1).*norm(yp)));

        % Finding the angle v1 makes with the z axis
        gamma = acos((dot(v1, zp))./(norm(v1).*norm(zp)));

        PDD(x, y, 1) = FA(x, y) * cos(alpha);

        PDD(x, y, 2) = FA(x, y) * cos(beta);

        PDD(x, y, 3) = FA(x, y) * cos(gamma);

    end
end

%% Plot mean diffusivity, fractional anisotropy and principal diffusion direction maps

% Plotting mean diffusivity map
MD = MD*255;
figure, imshow(MD, 'InitialMagnification', 'fit')

% Plotting fractional anisotropy map
figure, imshow(FA, 'InitialMagnification', 'fit')

% Plotting principal diffusion direction map
PDD = abs(PDD);
figure, imshow(PDD, 'InitialMagnification', 'fit')
