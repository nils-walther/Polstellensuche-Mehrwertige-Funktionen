f = @(z) sqrt(z.^2 - 1);

[X, Y] = meshgrid(-3:0.05:3, -3:0.05:3);
Z = X + 1i*Y;

figure;
pcolor(X, Y, angle(f(Z)));
shading interp;
colormap hsv;
axis equal tight;
title('Branch Cuts');