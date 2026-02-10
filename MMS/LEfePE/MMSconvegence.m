clear;

nels = 2.^(1:9);

sigErr = [1.0205e+00   5.1948e-01   2.5428e-01   1.2775e-01   6.4203e-02...
          3.2212e-02   1.6139e-02   8.0789e-03   4.0420e-03];
        
uErr   = [1.0932e-01   6.8969e-02   2.3588e-02   6.6602e-03   1.7696e-03...
          4.5620e-04   1.1583e-04   2.9184e-05   7.3245e-06];        

      
uEr8   = [1.1283e-01   1.1021e-02   9.2444e-04   6.3813e-05   4.1935e-06...
          2.6871e-07   1.7008e-08   1.0698e-09   6.6828e-11];  

sigEr8 = [8.5341e+13   8.5360e-02   2.4293e-02   5.9836e-03   1.4811e-03...
          3.6892e-04   9.2122e-05   2.3022e-05   5.7550e-06];
      
disEr8 = [6.1424e-01   1.5822e-02   1.8791e-03   2.2875e-04   2.8025e-05...
          3.4651e-06   4.3071e-07   5.3685e-08   6.7012e-09];
      
figure(1); clf;
loglog(nels,sigErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
loglog(nels(2:end),sigEr8(2:end),'k-s','MarkerSize',8,'MarkerFaceColor','k'); 
set(gca,'FontSize',16);
legend('4 noded element','8 noded element');
xlabel('1/h');
ylabel('error')


figure(2); clf;
loglog(nels,uErr,'k-o','MarkerSize',8,'MarkerFaceColor','w'); hold on;
loglog(nels,uEr8,'k-s','MarkerSize',8,'MarkerFaceColor','k');
loglog(nels,disEr8,'r-s','MarkerSize',8,'MarkerFaceColor','r');
set(gca,'FontSize',16);
legend('4 noded: nodal','8 noded: nodal','8 noded: volume');
xlabel('1/h');
ylabel('error')


% convergence rates
uRate   = mean((log(uErr(2:end))-log(uErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))
sigRate = mean((log(sigErr(2:end))-log(sigErr(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))

uRate8   = mean((log(uEr8(2:end))-log(uEr8(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))
sigRate8 = mean((log(sigEr8(2:end))-log(sigEr8(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1))))
disRate8 = (log(disEr8(2:end))-log(disEr8(1:end-1)))./(log(nels(2:end))-log(nels(1:end-1)));