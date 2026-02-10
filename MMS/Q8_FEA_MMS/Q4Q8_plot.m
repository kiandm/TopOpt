clear; clc;
nels      = [2 4 8 16 32 64 128 256];

%Q4
uErrQ4    = [1.19E+00 2.34E-01 6.24E-02 1.59E-02 4.01E-03 1.00E-03 2.51E-04 6.27E-05];
sigErrQ4  = [9.55E-01 4.35E-01 2.24E-01 1.13E-01 5.66E-02 2.83E-02 1.42E-02 7.09E-03];

uRateQ4   = mean((log(uErrQ4(2:end))-log(uErrQ4(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))
sigRateQ4 = mean((log(sigErrQ4(2:end))-log(sigErrQ4(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))

%Q8
uErrQ8 = [1.055997e+00 2.842328e-02 3.366641e-03 4.146436e-04 5.159117e-05 6.441009e-06 8.048750e-07 1.006015e-07];
sigErrQ8 = [9.502253e-01 1.068669e-01 2.347848e-02 5.768017e-03 1.437917e-03 3.592645e-04 8.980352e-05 2.245011e-05];

uRateQ8   = mean((log(uErrQ8(2:end))-log(uErrQ8(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))
sigRateQ8 = mean((log(sigErrQ8(2:end))-log(sigErrQ8(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))

figure; clf;
loglog(nels,uErrQ4,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
loglog(nels,uErrQ8,'k-s','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
legend('4 noded element', '8 noded element', 'Location','southwest');
title('Displacement error')
xlabel('1/h');
ylabel('error')

figure; clf;
loglog(nels,sigErrQ4,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
loglog(nels,sigErrQ8,'k-s','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
legend('4 noded element', '8 noded element', 'Location','southwest');
title('Stress error')
xlabel('1/h');
ylabel('error')

% uRateQ4 = -2.0303
% sigRateQ4 = -1.0105
% uRateQ8 = -3.2437
% sigRateQ8 = -2.1003

