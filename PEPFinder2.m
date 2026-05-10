clear; clc; close all;

global zero eps omega;
zero = 1e-20;          
omega = 5;           
eps = 5e-2;           
top_k = 30;            

function [is_intersect, pt] = segments_intersect(A, B, C, D)
    global zero;

    u = B - A;
    v = D - C;
    w = A - C;

    a = dot(u, u);
    b = dot(u, v);
    c = dot(v, v);
    d = dot(u, w);
    e = dot(v, w);
    denom = a * c - b * b;

    if abs(denom) < zero
        is_intersect = false;
        pt = [];
        return;
    end

    s = (b * e - c * d) / denom;
    t = (a * e - b * d) / denom;

    if (s > zero && s <= 1 - zero) && (t > zero && t <= 1 - zero)
        pt = A + s * u;
        is_intersect = true;
    else
        is_intersect = false;
        pt = [];
    end
end

function intersections = self_intersection(points, omega)
    global zero eps;
    N = size(points, 1);
    intersections = [];
    idx = 1;

    for i = 1:N-1
        A = points(i, :);
        B = points(i+1, :);

        for j = i + omega:N-1
            C = points(j, :);
            D = points(j+1, :);

            [intersect, pt] = segments_intersect(A, B, C, D);
            if intersect
                dis_mat = pdist2([A;B], [C;D]);
                [min_dist, min_idx] = min(dis_mat(:));
                [p1, p2] = ind2sub(size(dis_mat), min_idx);
                
                if min_dist < eps
                    intersections(idx).segments = [i, i+1, j, j+1];
                    intersections(idx).distance = min_dist;
                    intersections(idx).centroid = pt;
                    aa=[A;B];
                    bb=[C;D];
                    intersections(idx).points = [aa(p1,:); bb(p2,:)];
                    idx = idx + 1;
                end
            end
        end
    end
end


function display_curve(curve, intersections, top_k)
    global omega;
    figure('Position', [100, 100, 900, 700]);
    hold on; grid on; box on;

    plot3(curve(:,1), curve(:,2), curve(:,3), 'o-', 'LineWidth', 1, 'DisplayName', 'All');
    
    if ~isempty(intersections)
        [~, idx] = sort([intersections.distance]);
        intersections = intersections(idx);

        show_num = min(top_k, length(intersections));
        
        for i = 1:show_num
            intersect = intersections(i);
            segments = intersect.segments;
            centroid = intersect.centroid;
            dist = intersect.distance;
            points = intersect.points;
            
            i1 = segments(1); i2 = segments(2);
            j1 = segments(3); j2 = segments(4);

            i11=i2-omega;
            if i11<1
                i11=1;
            end
            i12=i1+omega;
            if i12>length(curve)
                i12=length(curve);
            end
            j11=j2-omega;
            if j11<1
                j11=1;
            end
            j12=j1+omega;
            if j12>length(curve)
                j12=length(curve);
            end
            
            niche=curve([[i11:i12],[j11:j12]],:);
            min_dim=min(niche);
            max_dim=max(niche);

            curve([i1,i2],:)
            curve([j1,j2],:)

            fprintf('相交线段：%d-%d 与 %d-%d\n', i1, i2, j1, j2);
            fprintf('交点：');
            fprintf('%.4f', centroid);
            fprintf('\n距离：%.6f\n\n', dist);

            plot3(curve([i1,i2],1), curve([i1,i2],2), curve([i1,i2],3), ...
                'g-', 'LineWidth', 2, 'DisplayName', 'Segement A');
            plot3(curve([j1,j2],1), curve([j1,j2],2), curve([j1,j2],3), ...
                'm-', 'LineWidth', 2, 'DisplayName', 'Segement B');

            plot3(points(1,1), points(1,2), points(1,3), ...
                'gd', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Point A');
            plot3(points(2,1), points(2,2), points(2,3), ...
                'md', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'Point B');

            xlabel('Chemical potential 1'); ylabel('Chemical potential 2'); zlabel('Chemical potential 3');
            title('The detection of equilibrium chemical potential pairs');
            legend('Location', 'best');
            view(3);
            rotate3d on;

            ax2=axes('Position', [0.3, 0.25, 0.25, 0.2]);
            grid on;
            box on;
            hold on;

            plot3(ax2,niche(:,1), niche(:,2), niche(:,3), ...
                'o-', 'MarkerSize', 8, 'MarkerFaceColor', 'b',  'DisplayName', 'Crossing Area');

            plot3(ax2,curve([i1,i2],1), curve([i1,i2],2), curve([i1,i2],3), ...
                'gd', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Point A');
            
            plot3(ax2,curve([j1,j2],1), curve([j1,j2],2), curve([j1,j2],3), ...
                'md', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'Point B');

            xlabel('CP 1'); ylabel('CP 2'); zlabel('CP 3');

            ax = gca;
            disableDefaultInteractivity(ax);              
 
            set(gcf,'WindowButtonUpFcn',@beginDrag);  
            set(ax2,'ButtonDownFcn',@(s,~)beginDrag(s));
        end
    else
        fprintf('❌ No intersections**\n');
    end
end

curve=readmatrix("data/data4.csv");
if (size(curve,2)<3)
    curve(:,3)=0;
end

intersections = self_intersection(curve, omega);

display_curve(curve, intersections, 1);
if isempty(intersections)
    fprintf('No crossed segments\n');
end