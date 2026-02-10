function xnew = OC_update(x, dc, volfrac)
    l1 = 0; l2 = 1e9; move = 0.2;
    xnew = x;
    while (l2 - l1) / (l2 + l1) > 1e-6
        lmid = 0.5*(l1 + l2);
        xcand = max(0.001, max(x - move, min(1.0, min(x + move, x .* sqrt(-dc ./ lmid)))));
        if mean(xcand) - volfrac > 0
            l1 = lmid;
        else
            l2 = lmid;
        end
        xnew = xcand;
    end
end
