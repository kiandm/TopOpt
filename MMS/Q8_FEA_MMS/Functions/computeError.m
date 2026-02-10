function errors = computeError(nodes, enodes, edof, U, U_exact, De, mms)
    % Initialize error accumulators
    l2_error_squared = 0;
    l2_exact_squared = 0;
    stress_l2_error_squared = 0;
    stress_exact_squared = 0;
    total_area = 0;  
    % Gauss quadrature - update to 3x3
    gp = [-sqrt(3/5) 0 sqrt(3/5)];
    w = [5/9 8/9 5/9];
    % Loop over elements
    for e = 1:size(enodes, 1)
        % Compute element error contributions
        elem_errors = computeElementError(e, nodes, enodes, edof, U, U_exact, De, mms, gp, w);       
        % Accumulate errors
        l2_error_squared = l2_error_squared + elem_errors.l2_error_squared;
        l2_exact_squared = l2_exact_squared + elem_errors.l2_exact_squared;
        stress_l2_error_squared = stress_l2_error_squared + elem_errors.stress_error_squared;
        stress_exact_squared = stress_exact_squared + elem_errors.stress_exact_squared;
        total_area = total_area + elem_errors.area;
    end   
    % Compute final errors
    errors.l2_error = sqrt(l2_error_squared);
    errors.l2_exact = sqrt(l2_exact_squared);
    errors.rel_l2_error = errors.l2_error / errors.l2_exact;
    errors.stress_error = sqrt(stress_l2_error_squared);
    errors.stress_exact = sqrt(stress_exact_squared);
    errors.rel_stress_error = errors.stress_error / errors.stress_exact;    
    errors.total_area = total_area;
end

function elem_errors = computeElementError(e, nodes, enodes, edof, U, ~, De, mms, gp, w)
    % Initialize element errors
    elem_errors.l2_error_squared = 0;
    elem_errors.l2_exact_squared = 0;
    elem_errors.stress_error_squared = 0;
    elem_errors.stress_exact_squared = 0;
    elem_errors.area = 0;   
    % Get element data
    dofs = edof(e,:);
    n = enodes(e,:);
    xe = nodes(n,2);
    ye = nodes(n,3);
    Ue = U(dofs);   
    % Gauss integration for error
    for ixi = 1:3
        xi = gp(ixi);
        wx = w(ixi);
        for ieta = 1:3
            eta = gp(ieta);
            wy = w(ieta);           
            % Shape functions and derivatives
            [N, dNdxi] = shapeFcn(xi, eta);           
            % Jacobian
            J = [dNdxi(1,:)*xe, dNdxi(1,:)*ye;
                 dNdxi(2,:)*xe, dNdxi(2,:)*ye];
            detJ = abs(det(J));
            elem_errors.area = elem_errors.area + detJ * wx * wy;            
            % Shape function derivatives in physical coordinates
            dNdx = J \ dNdxi;            
            % B matrix
            B = zeros(3,16);
            for i = 1:8
                B(:,2*i-1:2*i) = [ dNdx(1,i),            0;
                                   0,            dNdx(2,i);
                                   dNdx(2,i),  dNdx(1,i) ];
            end           
            % Physical coordinates at Gauss point
            x_gp = N * xe;
            y_gp = N * ye;           
            % FEA solution at Gauss point
            ux_fe_gp = N * Ue(1:2:16);
            uy_fe_gp = N * Ue(2:2:16);            
            % Exact solution at Gauss point
            ux_ex_gp = mms.ux_exact(x_gp, y_gp);
            uy_ex_gp = mms.uy_exact(x_gp, y_gp);            
            % FEA stress at Gauss point
            strain_fe = B * Ue;
            stress_fe = De * strain_fe;            
            % Exact stress at Gauss point
            stress_ex = [mms.sigma_xx(x_gp, y_gp);
                         mms.sigma_yy(x_gp, y_gp);
                         mms.sigma_xy(x_gp, y_gp)];            
            % Displacement errors
            error_x = ux_fe_gp - ux_ex_gp;
            error_y = uy_fe_gp - uy_ex_gp;            
            % Stress errors
            stress_error = stress_fe - stress_ex;
            stress_error_norm_sq = sum(stress_error.^2);
            stress_exact_norm_sq = sum(stress_ex.^2);            
            % Accumulate L2 displacement errors
            elem_errors.l2_error_squared = elem_errors.l2_error_squared + ...
                (error_x^2 + error_y^2) * detJ * wx * wy;            
            elem_errors.l2_exact_squared = elem_errors.l2_exact_squared + ...
                (ux_ex_gp^2 + uy_ex_gp^2) * detJ * wx * wy;            
            % Accumulate stress errors
            elem_errors.stress_error_squared = elem_errors.stress_error_squared + ...
                stress_error_norm_sq * detJ * wx * wy;       
            elem_errors.stress_exact_squared = elem_errors.stress_exact_squared + ...
                stress_exact_norm_sq * detJ * wx * wy;
        end
    end
end