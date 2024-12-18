function [sols, v, P2_R, C2_R] = relative_pose_conics_rectified_system(lines, e3q3)
% RELATIVE POSE CONICS RECTIFIED SYSTEM
% [sols, v, P2_R, C2_R] = relative_pose_conics_rectified_system(lines, e3q3, plots)
%
% lines     Rectified lines in struct with values:
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
% e3q3      Use the solver e3q3 from Kukelova et al.
% plots     Show plots

%% Calculate rotations
[P2_R] = calculateP2_R(lines);
[C2_R,v] = calculateC2_R(lines);
%[P2_R] = calculateP2_R_v2(lines, v);
%[P2_R] = mycalculateP2_R(lines);

%% Solver 1
data0(1:4) = rot2quat(P2_R);

% Lines cylinder1 (idx: camera, cylinder, line)
data0(5) = lines.l_111(1);
data0(6) = lines.l_211(1);
data0(7) = lines.l_212(1);

sols_P2_xz = solver_initialization_camera_translation(data0);

% Pick correct solution
P2_R_red = P2_R([1 3],[1 3]);

true_sol_xz = [];

for i = 1:4
    p = sols_P2_xz(:,i);
    % P = [R -R*t1];
    proj_c = [P2_R_red -P2_R_red*p]*[sols_P2_xz; ones(1,4)];
    proj_c(:,i) = [];
    if ~sum(proj_c(2,:) < 0)
        true_sol_xz = p;
    end
end

%true_sol_xz

sols_P2_xz = true_sol_xz;
%% Solver 2
    
    data1(1:4) = rot2quat(C2_R);
    data1(5:8) = rot2quat(P2_R); 
    
    % Lines cylinder2 (idx: camera, cylinder, line)
    data1(9:10) = lines.l_121(1:2); 
    data1(11:12) = lines.l_122(1:2);
    data1(13:14) = lines.l_221(1:2);
    data1(15:16) = lines.l_222(1:2);
    
    data1(17:18) = [sols_P2_xz(1) sols_P2_xz(2)];
    
    if size(data1,1) == 1
        data1 = data1';
    end
    
    [data_special, r_eq] = get_data_for_coeffs(data1);

    if e3q3 % Using the e3q3 solver from Kukelova et al.
        ds = data_special;
        
        ff = 50; % Normalisation konstant
        ee = [2 2 2 2 2 2 1 1 1 0];
        % The function takes a 3x10 matrix corresponding to the coefficients of the three quadratics.
        % The order of the monomials is
        % x^2,    xy,    xz,   y^2,    yz,    z^2,     x,      y,      z,    1.0;
        
        c1 = [     0,     0,     0, ds(1), ds(2),  ds(3),     0,      0,      0,      0]/min(ds(1:3));
        c2 = [ ds(4), ds(5), ds(6), ds(8), ds(9), ds(11), ds(7), ds(10), ds(12), ds(13)]/min(ds(4:13));
        c3 = [ds(14),ds(15),ds(16),ds(18),ds(19), ds(21),ds(17), ds(20), ds(22), ds(23)]/min(ds(14:23));
        
        
        c1 = c1.*(ff.^ee);
        c2 = c2.*(ff.^ee);
        c3 = c3.*(ff.^ee);
        
        c1 = c1/max(abs(c1));
        c2 = c2/max(abs(c2));
        c3 = c3/max(abs(c3));
        
        sols_e3q3 = re3q3_mex([c1; c2; c3]);
        sols2 = sols_e3q3*ff;
        
    else
        sols2 = solver_initialization_rest_of_variables_4vTo3v(data_special);
    end
    
    m = monomials(r_eq);
    c = coeffs(r_eq);
    
    r = -(c(1)*sols2(2,:).^2 + c(2)*sols2(2,:).*sols2(3,:) + c(3)*sols2(3,:).^2);
    
    n_sols = size(sols2,2);
    
    sols = [repmat(sols_P2_xz(1),1,n_sols); sols2(1,:); repmat(sols_P2_xz(2),1,n_sols); sols2(2:end,:); r];
end

function [P2_R] = calculateP2_R(lines)
% Calculate Camera rotation for P2 (and theta)
l = lines;

A1 = [l.l_121(1)*l.l_122(2) - l.l_121(2)*l.l_122(1) + l.l_121(2)*l.l_221(1) - l.l_122(2)*l.l_221(1), ...
    l.l_121(1)*l.l_122(2)*l.l_221(1) - l.l_121(2)*l.l_122(1)*l.l_221(1) - l.l_121(2) + l.l_122(2)];

y1 = -(- l.l_121(1)*l.l_221(2) + l.l_122(1)*l.l_221(2));

A2 = [l.l_121(1)*l.l_122(2) - l.l_121(2)*l.l_122(1) + l.l_121(2)*l.l_222(1) - l.l_122(2)*l.l_222(1), ...
    l.l_121(1)*l.l_122(2)*l.l_222(1) - l.l_121(2)*l.l_122(1)*l.l_222(1) - l.l_121(2) + l.l_122(2)];

y2 = -(- l.l_121(1)*l.l_222(2) + l.l_122(1)*l.l_222(2));

ab = [A1;A2]\[y1;y2];
ab = ab/norm(ab); % a och b

a = ab(1);
b = ab(2);

P2_R = [a 0 b; 0 1 0; -b 0 a];
end

function [P2_R] = calculateP2_R_v2(lines, v)
% normal of plane  perpendicular to vanishing point

% Multiply(Transpose(Pi_221), v2) = [(a*l_221_x - b)*v_x + l_221_y*v_y + (b*l_221_x + a)*v_z])
% =  (l_221_x*v_x + v_z)*a + (-v_x + l_221_x*v_z)*b + l_221_y*v_y
% Multiply(Transpose(Pi_222), v2) = [(a*l_222_x - b)*v_x + l_222_y*v_y + (b*l_222_x + a)*v_z])
% =  (l_222_x*v_x + v_z)*a + (-v_x + l_222_x*v_z)*b + l_222_y*v_y

l_221_x = lines.l_221(1);
l_221_y = lines.l_221(2);
l_222_x = lines.l_222(1);
l_222_y = lines.l_222(2);

v_x = v(1);
v_y = v(2);
v_z = v(3);

A = [(l_221_x*v_x + v_z)  (-v_x + l_221_x*v_z); ...
    (l_222_x*v_x + v_z)  (-v_x + l_222_x*v_z)];

y = -[l_221_y*v_y; l_222_y*v_y];

[ab] = A\y;
ab = ab/norm(ab); % a och b

a = ab(1);
b = ab(2);

P2_R = [a 0 b; 0 1 0; -b 0 a];
end

function [C2_R,v] = calculateC2_R(lines)
% Calculate C2_R, rotation of cylinder 2

% a = P2_R(1,1);
% b = P2_R(1,3);

% normals of planes tangenting T2
n_121 = lines.l_121;
n_122 = lines.l_122;
% n_221 = [a*lines.l_221(1) - b, lines.l_221(2), b*lines.l_221(1) + a]';
% n_222 = [a*lines.l_222(1) - b, lines.l_222(2), b*lines.l_222(1) + a]';

v = cross(n_121, n_122);
v = v/norm(v);
% pflat(cross(n_221, n_222)) == pflat(cross(n_121, n_122))

% Pick v and one normal as basis for rotation, e_y maps onto v

R_x = n_121/norm(n_121);
R_y = v;
R_z = cross(R_x, R_y);

C2_R  = [R_x, R_y, R_z]';

end

function [P2_R] = mycalculateP2_R(lines)
l1 = cross(lines.l_121,lines.l_122);
l2 = cross(lines.l_221,lines.l_222);
l1 = l1/norm(l1);
l2 = l2/norm(l2);
lx = l1(1);
%ly = l1(2);
lz = l1(3);
jx = l2(1);
%jy = l2(2);
jz = l2(3);
a = (jx*lx + jz*lz)/(lx^2 + lz^2);
b = -(jx*lz - jz*lx)/(lx^2 + lz^2);
n = sqrt(a^2+b^2);
a = a/n;
b = b/n;
R = [a 0 -b;0 1 0;b 0 a];
P2_R = R;
end

function [data_special, r_eq] = get_data_for_coeffs(data0)

nbr_unknowns = 3; % Variabler
nbr_generic_coeffs = 18; % length(data0)

if nargin < 1 || isempty(data0)
    % no input, generate a random integer instance
    data0 = randi(30, nbr_generic_coeffs, 1); % sänk från 30 om coeffs börjar bli för höga.
end

C2_R_quat = data0(1:4);
P2_R_quat = data0(5:8);

n = 9;

% Lines cylinder2 (camera, cylinder, line)
l_121 = [data0(n:n+1); 1]; n = n+2;
l_122 = [data0(n:n+1); 1]; n = n+2;
l_221 = [data0(n:n+1); 1]; n = n+2;
l_222 = [data0(n:n+1); 1];


xx = create_vars(nbr_unknowns+1);
P2_t = [data0(17) xx(1) data0(18)]';
% r1pow2 = l111(1)^2/(1+l111(1)^2)*p_z^2; % Likformighet, men får inte
% använda kvoter
T2_t = xx(2:3);
r2pow2 = xx(4);

% Cameras
P1 = [eye(3) zeros(3,1)];
R2 = quat2rot(P2_R_quat);
P2 = [R2 (-R2*P2_t)];

C2_R = quat2rot(C2_R_quat);
A = [C2_R' [T2_t(1); 0; T2_t(2)]; 0 0 0 1];
C_0inv = [-r2pow2  0  0       0; ...
    0       0  0       0; ...
    0       0 -r2pow2  0; ...
    0       0  0       1];

D2 = A*C_0inv*A';
% D2 = zp_reduce(D2, p);

% Equations (camera, cylinder, line)
temp_eq(1,1) = l_121'*P1*D2*P1'*l_121;
temp_eq(2,1) = l_122'*P1*D2*P1'*l_122;
temp_eq(3,1) = l_221'*P2*D2*P2'*l_221;
temp_eq(4,1) = l_222'*P2*D2*P2'*l_222;

% temp_eq = zp_reduce(temp_eq,p);
c1 = coeffs(temp_eq(1,1));
c2 = coeffs(temp_eq(2,1));
c3 = coeffs(temp_eq(3,1));
c4 = coeffs(temp_eq(4,1));

eqs(1,1) = temp_eq(2,1)*c1(4) - temp_eq(1,1)*c2(4);
eqs(2,1) = temp_eq(3,1)*c1(4) - temp_eq(1,1)*c3(10);
eqs(3,1) = temp_eq(4,1)*c1(4) - temp_eq(1,1)*c4(10);


cc1 = coeffs(eqs(1));
cc2 = coeffs(eqs(2));
cc3 = coeffs(eqs(3));

data_special(1:3) = cc1;
data_special(4:13) = cc2;
data_special(14:23) = cc3;

data_special = data_special';

r_eq = temp_eq(1,1)/c1(4);

end

