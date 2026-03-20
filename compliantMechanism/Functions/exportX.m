% rho_cut = 0.5;
% thickness = 3;
% 
% solid = x_phy_best_filtered > rho_cut;
% 
% V = [];    % vertices
% F = [];    % faces
% vo = 0;    % vertex offset
% 
% for e = 1:size(enodes,1)
% 
%     if ~solid(e)
%         continue
%     end
% 
%     % Node coordinates of the element
%     el_nodes = enodes(e,:);
%     xy = nodes(el_nodes, 2:3);   % (4 × 2)
% 
%     % Bottom surface (z = 0)
%     vb = [xy, zeros(4,1)];
% 
%     % Top surface (z = thickness)
%     vt = [xy, thickness*ones(4,1)];
% 
%     % Append vertices
%     V = [V; vb; vt];
% 
%     % Local indices
%     b = vo + (1:4);
%     t = vo + (5:8);
% 
%     % Triangulate bottom quad
%     Fb = [
%         b(1) b(2) b(3);
%         b(1) b(3) b(4)
%     ];
% 
%     % Triangulate top quad
%     Ft = [
%         t(1) t(3) t(2);
%         t(1) t(4) t(3)
%     ];
% 
%     % Side wall triangles
%     Fs = [
%         b(1) b(2) t(2); b(1) t(2) t(1);
%         b(2) b(3) t(3); b(2) t(3) t(2);
%         b(3) b(4) t(4); b(3) t(4) t(3);
%         b(4) b(1) t(1); b(4) t(1) t(4)
%     ];
% 
%     F = [F; Fb; Ft; Fs];
%     vo = vo + 8;
% 
% end
% 
% TR = triangulation(F, V);
% stlwrite('topopt.stl', TR.ConnectivityList, TR.Points);

%% =========================
%  INPUTS
%  nodes   = [id, x, y]
%  enodes  = [n1 n2 n3 n4]
%  x       = density per element
% ==========================
rho_cut = 0.5;
gaussian_sigma = 4;     % smoothing amount (increase for smoother boundaries)
img_res = 800;          % resolution of rasterized level set
thickness = 3;          % extrusion height

%% =========================
%  1. Compute centroids of all FE elements
% =========================
E = length(x_phy_best_filtered);
cx = zeros(E,1);
cy = zeros(E,1);

for e = 1:E
    vx = nodes(enodes(e,:), 2);
    vy = nodes(enodes(e,:), 3);
    cx(e) = mean(vx);
    cy(e) = mean(vy);
end

solid = x_phy_best_filtered > rho_cut;

%% =========================
%  2. Rasterize onto a uniform grid for level-set processing
% =========================
xi = linspace(min(cx), max(cx), img_res);
yi = linspace(min(cy), max(cy), img_res);
[XI, YI] = meshgrid(xi, yi);

% nearest-neighbour assignment of densities to grid

F_interp = scatteredInterpolant(cx, cy, double(solid), 'nearest', 'nearest');
Z = F_interp(XI, YI);
Z(isnan(Z)) = 0;

%% =========================
%  3. Gaussian smoothing of the binary mask 
% =========================
Zsmooth = imgaussfilt(double(Z), gaussian_sigma);

% threshold to recover a smooth boundary
bw = Zsmooth > 0.5;

%% =========================
%  4. Extract smoothed contour curve
% =========================
C = contourc(double(bw), [0.5 0.5]);
% parse contour output
idx = 2 : (C(2,1) + 1);
poly = C(:, idx);        % poly(1,:) = x coords, poly(2,:) = y coords

px = XI(1, round(poly(1,:)));  % map contour to real-world coordinates
py = YI(round(poly(2,:)), 1);

%% =========================
%  5. Build alphaShape from smooth boundary
% =========================
shp = alphaShape(px(:), py(:), 2.0);
shp = alphaShape(shp.Points, shp.Alpha*1.1);  % optional smoothing tighten

[tri, pts2] = alphaTriangulation(shp);

%% =========================
%  6. Extrude shape into 3D
% =========================
numPts = size(pts2,1);

pts_bottom = [pts2, zeros(numPts,1)];
pts_top    = [pts2, thickness*ones(numPts,1)];

V = [pts_bottom; pts_top];

tri_bottom = tri;
tri_top = tri + numPts;

% Build side wall triangles
side_faces = [];
for k = 1:size(tri,1)
    t = tri(k,:);
    a = t(1); b = t(2); c = t(3);
    A = a + numPts; B = b + numPts; C = c + numPts;

    side_faces = [side_faces;
        a b B;   a B A;
        b c C;   b C B;
        c a A;   c A C];
end

F = [tri_bottom; tri_top; side_faces];

%% =========================
%  7. Export final smooth STL
% =========================
TR = triangulation(F, V);
stlwrite('smooth_topopt.stl', TR.ConnectivityList, TR.Points);

disp("Smooth STL written to smooth_topopt.stl");
