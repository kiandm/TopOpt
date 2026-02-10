clear; clc;
nels = [2 4 8 16 32 64 128 256];
% 2x2 Gauss
%uErr = [1.146898e+00 3.432782e-02 5.234100e-03 6.791710e-04 8.561134e-05 1.072329e-05 1.341091e-06 1.676576e-07];
%sigErr = [9.292097e-01 5.957046e-02 6.090875e-03 7.059150e-04 8.644716e-05 1.074943e-05 1.341908e-06 1.676832e-07];

% 3x3 Gauss - not much difference
uErr = [1.055997e+00 2.842328e-02 3.366641e-03 4.146436e-04 5.159117e-05 6.441009e-06 8.048750e-07 1.006015e-07];
sigErr = [9.502253e-01 1.068669e-01 2.347848e-02 5.768017e-03 1.437917e-03 3.592645e-04 8.980352e-05 2.245011e-05];

figure(1); clf;
loglog(nels,uErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
%legend('4 noded element','8 noded element');
xlabel('1/h');
ylabel('error')

figure(2); clf;
loglog(nels,sigErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
%legend('4 noded element','8 noded element');
xlabel('1/h');
ylabel('error')

uRate   = mean((log(uErr(2:end))-log(uErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))
sigRate = mean((log(sigErr(2:end))-log(sigErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))