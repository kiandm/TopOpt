function plotGeometry(nodes, enodes, BCs_xy, LC, force_in, disp_out)
% nodes: [id, x, y, dofx, dofy]
% enodes: element connectivity (Q4)
% BCs_xy: nodes with supports
% force_in: input force nodes
% disp_out: output displacement nodes

figure; hold on; axis equal; box on;
title('Geometry with BCs and Loads');
xlabel('X'); ylabel('Y');

%% --- Plot mesh (elements)
nel = size(enodes,1);
for e = 1:nel
    n = enodes(e,:);
    xy = nodes(n,2:3);
    patch(xy(:,1), xy(:,2), [0.9 0.9 0.9], ...
        'EdgeColor',[0.6 0.6 0.6], 'FaceAlpha',0.3);
end

%% --- Plot support BCs (blue squares)
if ~isempty(BCs_xy)
    plot(nodes(BCs_xy,2), nodes(BCs_xy,3), ...
        'bs','MarkerFaceColor','b','MarkerSize',8, ...
        'DisplayName','Roller Supports');
end
if ~isempty(LC)
    plot(nodes(LC,2), nodes(LC,3), ...
        'bs','MarkerFaceColor','b','MarkerSize',8, ...
        'DisplayName','Supports');
end

%% --- Plot input load nodes (red arrows)
if ~isempty(force_in)
    for k = 1:numel(force_in)
        x0 = nodes(force_in(k),2);
        y0 = nodes(force_in(k),3);

        % Example: downward arrow (scale 0.2)
        quiver(x0, y0, 1, 0, ...
            'Color','r','LineWidth',2, ...
            'MaxHeadSize',3, 'DisplayName','Input load');
    end
end

%% --- Plot output nodes (green arrows)
if ~isempty(disp_out)
    for k = 1:numel(disp_out)
        x0 = nodes(disp_out(k),2);
        y0 = nodes(disp_out(k),3);

        % Example: upward arrow for desired output direction
        quiver(x0, y0, -1, 0, ...
            'Color','g','LineWidth',2, ...
            'MaxHeadSize',3, 'DisplayName','Output direction');
    end
end

%% Legend
legend({'Elements','Supports','Input Load','Output Port'}, ...
        'Location','southoutside');

end