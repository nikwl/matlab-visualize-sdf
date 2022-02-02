
% Generate our polygons
complete_gon = polyshape([0 0 1 1], [1 0 0 1]); % Square
% complete_gon = polyshape(circle(0.5, [0.5, 0.5])); % Circle
% primitive_gon = polyshape([0.0 0.0 0.5 0.5],[1.0 0.0 0.0 1.0]); % Just big enough primitive
primitive_gon = polyshape([-0.1 -0.1 0.5 0.5],[1.1 -0.1 -0.1 1.1]); % Slightly bigger primitive
% primitive_gon = polyshape([-5 -5 0.5 0.5],[1.1 -0.1 -0.1 1.1]); % Medium primitive
% primitive_gon = polyshape([-5 -5 0.5 0.5],[5 -5 -5 5]); % Really big primitive
% primitive_gon = polyshape([-5 -5 100 0.5 0.5 100],[5 -5 -5 0 1 5]); % Really big primitive that wraps around the complete object

% Compute the subraction 
fractured_gon = subtract(complete_gon, primitive_gon);

% Get sdf/occ values
[qpx, qpy, X, Y] = query_points(128);
C_sdf = get_sdf(complete_gon, qpx, qpy);
F_sdf = get_sdf(fractured_gon, qpx, qpy);
P_sdf = get_sdf(primitive_gon, qpx, qpy);

% Compute max(C, -T)
[max_v, max_arg] = max([C_sdf, -P_sdf], [], 2);

% Plot names, values, and isosurfaces
% names = ["sdf(C)", "sdf(B)", "sdf(T)", "max(sdf(C), -sdf(T))"];
% values = [C_sdf, B_sdf, T_sdf, max_v];
names = ["sdf(C)", "sdf(F)", "sdf(P)", "max(sdf(C), -sdf(P))", "error", "argmax C = (0) -P = (1)"];
values = [C_sdf, F_sdf, P_sdf, max_v,  abs(F_sdf - max_v), max_arg-1];

% Levels at which to extract the isosurface
levels = [0.0, 0.0, 0.0, 0.0, 0.0, 0.5];

% Make all plots use same color scale?
rescale = [true, true, true, true, false, false];

% Overlay the polygons on the plot in red
overlay_polygons = [false, false, false, true, true, true];

% Plot
figure
for i = 1:size(values,2)
    
    subplot(2,3,i);
    
    % Plot the values as a filled contour
    [~,h] = contourf(...
        reshape(qpx, [length(X), length(X)]),...
        reshape(qpy, [length(X), length(X)]),...
        reshape(values(:,i), [length(X), length(X)]), 50);
    set(h,'LineColor','none');
    title(names(i));
    hold on
    
    if (overlay_polygons(i))
        p1 = plot(primitive_gon);
        p1.EdgeColor = 'r';
        p1.FaceColor = 'none';
        p1.LineWidth = 2.0;

        p2 = plot(complete_gon);
        p2.EdgeColor = 'g';
        p2.FaceColor = 'none';
        p2.LineWidth = 2.0;
    end
    
    % Plot the isosurface extracted at the given label
    [C,h] = contour(...
        reshape(qpx, [length(X), length(X)]),...
        reshape(qpy, [length(X), length(X)]),...
        reshape(values(:,i), [length(X), length(X)]), [levels(i), levels(i)]);
    set(h,'LineColor','k');
    
    % Rescale the color values so that they all fall in the same range
    if (rescale(i))
        caxis([min(min(values)), max(max(values))]);
    end
    colorbar
    xlim([-0.2, 1.2]);
    ylim([-0.2, 1.2]);
     
    if (overlay_polygons(i))
        legend([p1, p2, h], 'Primitive', 'Complete', 'Isosurface');
    else
        legend([h], 'Isosurface');
    end
end

%% FUNCTIONS

function [coords] = circle(radius, center)
[x0,y0] = deal(center(1),center(2)); % Circle radius and center
t=linspace(0,360,1000).'; t(end)=[]; % circle angular samples
coords= [cosd(t), sind(t)]*radius+[x0,y0];
end

% Generate query points
function [qpx, qpy, X, Y] = query_points(density)
x_scan = -0.2:(1/density):1.2;
y_scan = -0.2:(1/density):1.2;
[X,Y] = meshgrid(x_scan,y_scan);
qpx = reshape(X, [length(X)^2, 1]);
qpy = reshape(Y, [length(X)^2, 1]);
end

% The sigmoid function
function [y] = sigmoid(x)
e = 1; % Exponent to raise sigmoid to
y = 1./(1 + exp(-1.*(x)).^e);
end

% The signed distance function
function [d_min] = get_sdf(gon, qpx, qpy)
d_min = p_poly_dist(qpx, qpy, gon.Vertices(:,1), gon.Vertices(:,2));
end

% The truncated signed distance function
function [d_min] = get_tsdf(gon, qpx, qpy, trunc)
d_min = p_poly_dist(qpx, qpy, gon.Vertices(:,1), gon.Vertices(:,2));
d_min = clamper(d_min, trunc);
end

% Clips the values between -/+ clipv
function [arr] = clamper(arr, clipv)
arr(arr>clipv) = clipv;
arr(arr<-clipv) = -clipv;
end

% The occupancy function
function [d_min] = get_occ(gon, qpx, qpy)
d_min = p_poly_dist(qpx, qpy, gon.Vertices(:,1), gon.Vertices(:,2));
d_min(d_min >= 0) = 0;
d_min(d_min < 0) = 1;
end
