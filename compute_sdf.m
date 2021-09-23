
% Generate our three polygons
complete_gon = polyshape([0 0 1 1], [1 0 0 1]);
broken_gon = polyshape([0 0 0.5 0.5],[1 0 0 1]);
restoration_gon = subtract(complete_gon, broken_gon);

% Generate query points
density = 128;
x_scan = -0.2:(1/density):1.2;
y_scan = -0.2:(1/density):1.2;
[X,Y] = meshgrid(x_scan,y_scan);
qpx = reshape(X, [length(X)^2, 1]);
qpy = reshape(Y, [length(X)^2, 1]);

% Get sdf/occ values
C_sdf = get_sdf(complete_gon, qpx, qpy);
B_sdf = get_sdf(broken_gon, qpx, qpy);
R_sdf = get_sdf(restoration_gon, qpx, qpy);

% C_occ = get_occ(complete_gon, qpx, qpy);
% R_occ = get_occ(restoration_gon, qpx, qpy);

% Raise sigmoid to this exponent
e = 1; 

% Plot names, values, and isosurfaces
names = ["sigmoid(sdf(C))", "sigmoid(sdf(R))", "-(sigmoid(sdf(C)) - sigmoid(sdf(R)))", "-(sigmoid(sdf(B)) - 1)"];
values = [sigmoid(C_sdf, e), sigmoid(R_sdf, e), -(sigmoid(C_sdf, e) - sigmoid(R_sdf, e)), -(sigmoid(B_sdf, e) - 1)];
levels = [0.5, 0.5, 0.5, 0.5];

% Make all plots use same color scale?
rescale = false;

% Overlay the polygons on the plot
overlay_polygons = false;

figure
for i = 1:size(values,2)
    subplot(2,2,i);
    
    % Plot the values as a filled contour
    [~,h] = contourf(...
        reshape(qpx, [length(X), length(X)]),...
        reshape(qpy, [length(X), length(X)]),...
        reshape(values(:,i), [length(X), length(X)]), 50);
    set(h,'LineColor','none');
    title(names(i));
    hold on
    
    % Plot the restoration polygon
    if (overlay_polygons)
        p = plot(restoration_gon);
        p.EdgeColor = 'r';
        p.FaceColor = 'none';

        % Plot the broken polygon
        p = plot(broken_gon);
        p.EdgeColor = 'r';
        p.FaceColor = 'none';
    end
    
    % Plot the isosurface extracted at the given label
    [C,h] = contour(...
        reshape(qpx, [length(X), length(X)]),...
        reshape(qpy, [length(X), length(X)]),...
        reshape(values(:,i), [length(X), length(X)]), [levels(i), levels(i)]);
    set(h,'LineColor','k');
    
    % Rescale the color values so that they all fall in the same range
    if (rescale)
        caxis([min(min(values)), max(max(values))]);
    end
    colorbar
end

%% 
figure
x = linspace(-5, 5, 1000);
e = 100;
y = sigmoid(x, e);
plot(x, y);
title("sigmoid");

%%

% The sigmoid function
function [y] = sigmoid(x, e)
y = 1./(1 + exp(-1.*(x)).^e);
end

% The signed distance function
function [d_min] = get_sdf(gon, qpx, qpy)
d_min = p_poly_dist(qpx, qpy, gon.Vertices(:,1), gon.Vertices(:,2));
end

% The signed distance function
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
