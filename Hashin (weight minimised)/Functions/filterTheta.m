function th_phys = filterTheta(th, H, Hs)
% pi-periodic fibre-orientation filter: average the DOUBLE angle, then halve.  Result is in (-pi/2, pi/2].
% (Fibre orientation theta and theta+pi are the same physical fibre; averaging cos(theta), sin(theta) is wrong across 0/pi.)
th_phys = 0.5*atan2((H*sin(2*th))./Hs, (H*cos(2*th))./Hs);
end