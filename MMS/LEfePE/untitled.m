nels = [2 4 8 16 32 64 128 256];
uErr = [1.19E+00 2.34E-01 6.24E-02 1.59E-02 4.01E-03 1.00E-03 2.51E-04 6.27E-05];
sigErr = [9.55E-01 4.35E-01 2.24E-01 1.13E-01 5.66E-02 2.83E-02 1.42E-02 7.09E-03];

figure(1); clf;
loglog(nels,uErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
legend('4 noded element');
title('log-log displacement error convergence')
xlabel('1/h');
ylabel('relative L2error')

figure(2); clf;
loglog(nels,sigErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
set(gca,'FontSize',16);
legend('4 noded element');
title('log-log stress error convergence')
xlabel('1/h');
ylabel('relative L2 error')

uRate   = (log(uErr(2:end))-log(uErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1)))
sigRate = (log(sigErr(2:end))-log(sigErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1)))
