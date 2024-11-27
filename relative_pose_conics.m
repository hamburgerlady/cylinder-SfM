function [v1, v2, rotC1, rotC2, results] = relative_pose_conics(linesOrig, P1_orig, scaleval, scaletype, flip_rectify, e3q3)
% RELATIVE POSE CONICS
% [v1, v2, rotC1, rotC2, results, locked_dist_results] = relative_pose_conics(linesOrig, P1_orig, scaleval, scaletype, e3q3, plots)
%
% linesOrig Rectified lines in struct with values:
%             l_111: [3×1 double]
%             l_112: [3×1 double]
%             l_211: [3×1 double]
%             l_212: [3×1 double]
%             l_121: [3×1 double]
%             l_122: [3×1 double]
%             l_221: [3×1 double]
%             l_222: [3×1 double]
%           where line l_ijk is from camera i, cylinder j and the k:th
%           line. Normalise such that l_ijk(3) == 1.
% P1_orig   Absolute pose of first camera.
% scaleval  Value of scale setting parameter.
% scaletype Type of scale setting parameter. 
%             'p2p_dist' distance between cameras.
%             'C1_r'     radius of first cylinder.
% e3q3      Use the solver e3q3 from Kukelova et al.
% plots     Show plots

% Rectify
[linesRec, rotC1, rotC2] = rectify_data(linesOrig, flip_rectify); % relative rotations

% Run solver
[solsRec, v, P2_R_rec, T2_R_rec] = relative_pose_conics_rectified_system(linesRec, e3q3);
% sols: P2_t (x,y,z); T2_t (x,z) y=0, r_2^2
% 32 solutions in total

% Return to original system. Check and discard solutions
[v1, v2, results] = get_results_orig_sys(P1_orig, scaleval, scaletype, solsRec,linesRec, rotC1, rotC2, v,  P2_R_rec, T2_R_rec);

end

% Subfunctions
function [v1, v2, results] = get_results_orig_sys(P1_orig, scaleval, scaletype, solsRec, linesRec, rotC1, rotC2, v,  P2_R_rec, T2_R_rec)
    P1_R = P1_orig(1:3,1:3);
P1_t = pflat(null(P1_orig));

% v1
v1 = P1_R'*rotC1'*[0 1 0]';

% v2
v2 = P1_R'*rotC1'*v(1:3);

if isempty(solsRec)
    results = [];
    return
else
    results = cell(size(solsRec,2),1);
    
    for i = 1:size(solsRec,2)
        sol = solsRec(:,i);
        
        scale = calculate_scale(scaleval, scaletype, sol, linesRec);

        results{i}.scale = scale;
        
        % P2
        P2_R = rotC2'*P2_R_rec*rotC1*P1_R;
        P2_t_rec = sol(1:3);
        P2_t = P1_R'*rotC1'*P2_t_rec*scale + P1_t(1:3);
        results{i}.P2   = [P2_R -P2_R*P2_t];
        
        % T1_t
        T1_t_rec = [0 0 1]';
        results{i}.T1_t     = P1_R'*rotC1'*T1_t_rec*scale + P1_t(1:3);
        
        % T2_t
        T2_t_rec = [sol(4) 0 sol(5)]';
        results{i}.T2_t     = P1_R'*rotC1'*T2_t_rec*scale + P1_t(1:3);
        
        % r2
        results{i}.r2 = sqrt(sol(6))*scale;
        
        % r1
        x = abs(1/linesRec.l_111(1));
        r1_rec = abs(x)/sqrt(x^2+1);
        results{i}.r1 = r1_rec*scale;
        
        % T2_R
        results{i}.T2_R = T2_R_rec*rotC1*P1_R; % True?
        
        % T1_R
        results{i}.T1_R = eye(3)*rotC1*P1_R;
        
    end
    
    
end
end

function [ll_rec, rot] = rectify_lines(guidelines, ll, flip)

[rot1, rot2] = find_rotation(guidelines);

% Original order
x1 = pflat(cross(guidelines(:,1), [0 1 0]'));
x2 = pflat(cross(guidelines(:,2), [0 1 0]'));

ll_cyl = rot1*[guidelines(:,1) guidelines(:,2)];

% After rectification
x1_hat = pflat(cross(ll_cyl(:,1), [0 1 0]'));
x2_hat = pflat(cross(ll_cyl(:,2), [0 1 0]'));

if sign(x1(1)-x2(1)) == sign(x1_hat(1) -x2_hat(1)) % check same order
    rot = rot1;
    if flip
        rot = rot2;
    end
    %disp('rot1')
else
    rot = rot2;
    if flip
        rot = rot1;
    end
    %disp('rot2')
end

ll_rec = rot*ll;


end

function [lines, rotC1, rotC2, rotC1_1] = rectify_data(lines, flip_rectify)

%Cam 1 rectify
[llcam1_pre, rotC1_1] = rectify_lines([lines.l_111 lines.l_112], ...
    [lines.l_111 lines.l_112 lines.l_121 lines.l_122],0);

% cam 1 center
rotC1_2 = find_rotation_center(llcam1_pre(:,1:2));
rotC1 = rotC1_2*rotC1_1;
llcam1 = pflat(rotC1_2*llcam1_pre);

% cam 2 rectify
[llcam2, rotC2] = rectify_lines([lines.l_211 lines.l_212], ...
    [lines.l_211 lines.l_212 lines.l_221 lines.l_222], flip_rectify);
llcam2 = pflat(llcam2);

% Save in lines variable
lines.l_111 = llcam1(:,1);
lines.l_112 = llcam1(:,2);
lines.l_121 = llcam1(:,3);
lines.l_122 = llcam1(:,4);

lines.l_211 = llcam2(:,1);
lines.l_212 = llcam2(:,2);
lines.l_221 = llcam2(:,3);
lines.l_222 = llcam2(:,4);

end

function rot = find_rotation_center(lines)
lines = pflat(lines);

% ax + 0y + 1 = 0
x1 = -1/lines(1,1);
x2 = -1/lines(1,2);

if x2<x1
    temp = x1; x1 = x2; x2 = temp;
end

theta = atan(x2) - atan(x1);
phi = -(atan(x1) + 1/2*theta);
rot = [cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];

end

function [Rvp1, Rvp2] = find_rotation(ll)
% Camera 1
[~,~,V] = svd(ll');
vp = V(:,end); % vanishing point
if vp(2)<0
    vp=-vp;
end
r2p = vp/norm(vp);

% Choosing a basis
z = [0 0 1]';
zp = (r2p'*z)*r2p;
zn = z-zp;
zn = zn/norm(zn);

Rvp1 = [cross(r2p,zn)   r2p zn]';
Rvp2 = [cross(-r2p,zn) -r2p zn]';

end

function scale = calculate_scale(scaleval, scaletype, sol, linesRec)

if strcmp(scaletype, 'p2p_dist')
    d_rec = norm(sol(1:3));
    scale = scaleval/d_rec;
elseif strcmp(scaletype, 'C1_r')
    x = abs(1/linesRec.l_111(1));
    r1_rec = abs(x)/sqrt(x^2+1);
    scale = scaleval/r1_rec;
else
    keyboard;
end

end
