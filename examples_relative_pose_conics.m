function [min_errs] = examples_relative_pose_conics(type, noise, e3q3, plots, seed)
% EXAMPLES RELATIVE POSE CONICS
% min_errs = examples_relative_pose_conics(type, noise, e3q3, seed, dbug, plots)
%
% type =    1 random synthetic data
%           2 Rollercoaster data
%
% noise     Add noise to lines
% e3q3      Use e3q3 solver (Kukelova et al) for second problem.
% plots     show plots
% seed      seed for random number generator

% Handle inputs
if nargin > 4
    rng(seed)
    disp(['For reference: seed = ' num2str(seed)]);
end

if nargin < 4, plots = 0; end
if nargin < 3, e3q3  = 0; end
if nargin < 2, noise = 0; end
if nargin < 1, type  = 1; end

addpath(genpath('../cylinder-SfM'))

%% Setup and pre-calculations

if type == 1
    [trueValues, linesOrig] = setup_synthetic(plots); % Seems to work
    
    % add noise
    ll1n = addnoise([linesOrig.l_111 linesOrig.l_112 linesOrig.l_121 linesOrig.l_122], noise, plots, 11);
    ll2n = addnoise([linesOrig.l_211 linesOrig.l_212 linesOrig.l_221 linesOrig.l_222], noise, plots, 12);
    linesOrig.l_111 = pflat(ll1n(:,1)); % Planes from cylinder nr 1 in camera 1
    linesOrig.l_112 = pflat(ll1n(:,2));
    linesOrig.l_121 = pflat(ll1n(:,3)); % Planes from cylinder nr 2 in camera 1
    linesOrig.l_122 = pflat(ll1n(:,4));
    
    linesOrig.l_211 = pflat(ll2n(:,1)); % Planes from cylinder nr 1 in camera 2
    linesOrig.l_212 = pflat(ll2n(:,2));
    linesOrig.l_221 = pflat(ll2n(:,3)); % Planes from cylinder nr 2 in camera 2
    linesOrig.l_222 = pflat(ll2n(:,4));
    
    % Send possible scale parameter: distance between cameras
    t1 = pflat(null(trueValues.P1));
    t2 = pflat(null(trueValues.P2));
    d_true = norm(t1-t2);
    scaleval = d_true;
    scaletype = 'p2p_dist';
    
    % Send possible scale parameter: radius cylinder 1
    scaleval = trueValues.T1_r;
    scaletype = 'C1_r';
    
    flip_rectify = 0;
    [v1, v2, rotC1, rotC2, results] = relative_pose_conics(linesOrig, trueValues.P1,  scaleval, scaletype, flip_rectify, e3q3);
    [min_errs, idx] = check_answers(trueValues, v2, rotC1, rotC2, results, plots);
    
    if max(min_errs > 10^-3)
        flip_rectify = 1; % Rectification of cameras in different directions. (edge case)
        [v1, v2, rotC1, rotC2, results] = relative_pose_conics(linesOrig, trueValues.P1,  scaleval, scaletype, flip_rectify, e3q3);
        [min_errs2, idx2] = check_answers(trueValues, v2, rotC1, rotC2, results, plots);
        
        [~,pick] = min([sum(min_errs), sum(min_errs2)]);
        if pick == 2
            min_errs = min_errs2;
            disp('second option')
        else
            disp('first option')
        end
    end
    
else % type = 2
    no_cyl = 5;
    radii = zeros(no_cyl,no_cyl); % Save estimated radii
    for ii = 1:no_cyl % Cylinder 1
        for jj = 1:no_cyl % cylinder 2
            if ii ~= jj
                [trueValues, linesOrig] = setup_real(ii, jj, 1);
                
                % Send possible scale parameter
                t1 = pflat(null(trueValues.P1));
                t2 = pflat(null(trueValues.P2));
                d_true = norm(t1-t2);
                scaleval = d_true;
                scaletype = 'p2p_dist';
                
                % Run solvers
                [v1, v2, rotC1, rotC2, results] = relative_pose_conics(linesOrig, trueValues.P1,  scaleval, scaletype, 0,  e3q3);
                
                %evaluate
                for i = 1:length(results)
                    [err_P_translation(i)] = error_translation_camera(results{i}.P2, trueValues);
                    [err_rotation_P2(i)] = error_rotation_camera2(results{i}.P2, trueValues); 
                end
                
                [value, idx] = min(err_P_translation);
                min_errs(:,ii,jj) = [err_P_translation(idx) err_rotation_P2(idx)]';
                radii(ii,jj) = results{idx}.r2;
                
                if plots
                    P2_est   = results{idx}.P2;
                    T1_t_est = results{idx}.T1_t;
                    T2_t_est = results{idx}.T2_t;
                    r2_est   = results{idx}.r2;
                    r1_est   = results{idx}.r1;
                    T2_R_est = results{idx}.T2_R;
                    T1_R_est = results{idx}.T1_R;
                    
                    % plots 15
                    start = -10;
                    finish = 30;
                    dist1 = 0.03;
                    dist2 = 1;
                    
                    pts_1 = cylinder_to_pts(T1_R_est, T1_t_est, r1_est, dist1, start, finish);
                    pts_2 = cylinder_to_pts(T2_R_est, T2_t_est, r2_est, dist2, start, finish);
          
                    % Plot some cylinders
                    if sum([ii jj] == [1 5]) == 2 || sum([ii jj] == [2 3]) == 2
                        figure(ii*100+jj*10+1)
                        hold on;
                        impts1{1} = pflat(trueValues.K*trueValues.P1*pts_1); % Project cylindrs into camera
                        rita(impts1{1},'.'); % draw cylinders
                        title('Estimated cylinders in camera1')
                        
                        impts1{2} = pflat(trueValues.K*trueValues.P1*pts_2); % Project cylindrs into camera
                        rita(impts1{2},'.'); % draw cylinders
                        axis([0 1920 0 1080])
                        axis ij
                    end
                    
                    if sum([ii jj] == [1 5]) == 2 || sum([ii jj] == [2 3]) == 2
                        figure(ii*100+jj*10+2)
                        hold on;
                        impts1{1} = pflat(trueValues.K*P2_est*pts_1); % Project cylindrs into camera
                        rita(impts1{1},'.'); % draw cylinders
                        title('Estimated cylinders in estimated camera2')
                        
                        impts1{2} = pflat(trueValues.K*P2_est*pts_2); % Project cylindrs into camera
                        rita(impts1{2},'.'); % draw cylinders
                        axis([0 1920 0 1080])
                        axis ij
                        
                    end
                end
            end
        end
    end
    variance_radius = var(radii)
    radii
end
end

function [min_errs, idx] = check_answers(trueValues, v2, rotC1, rotC2, results, plots)
if isempty(results)
    min_errs = NaN
else
    % check values
    P1 = trueValues.P1;
    P2 = trueValues.P2;
    
    % Rectified system
    P1_rec_true = rotC1*P1;
    P2_rec_true = rotC2*P2;
    
    % check rectified cameras
    if plots
        figure(4);
        hold on
        color = [0.5 0 0];
        R1 = P1_rec_true(1:3,1:3);
        R2 = P2_rec_true(1:3,1:3);
        t1 = pflat(null(P1_rec_true));
        t2 = pflat(null(P2_rec_true));
        pose = rigid3d(R1,t1(1:3)');
        plotCamera('AbsolutePose',pose,'Opacity',0, 'AxesVisible', 1, 'Color', color,'Label', 'Rec P1');
        pose = rigid3d(R2,t2(1:3)');
        plotCamera('AbsolutePose',pose,'Opacity',0, 'AxesVisible', 1, 'Color', color, 'Label', 'Rec P2');
        hold off
    end
    
    for i = 1:length(results)
        err_radius(i) = error_radius(results{i}.r2, trueValues); % OK
        err_vp(i) = error_degrees_vp(v2, trueValues); %  OK
        
        err_C_translation(i) = error_translation_cylinders(results{i}.T2_t, v2, trueValues);
        err_P_translation(i) = error_translation_camera(results{i}.P2, trueValues);
        err_rotation_P2(i) = error_rotation_camera2(results{i}.P2, trueValues); % Not a good measurement
    end
    
    if plots
        figure(5); clf;
        title('errors')
        hold on
        plot(log(err_radius), 'DisplayName', 'radius')
        plot(log(err_vp), 'DisplayName', 'vanishing p.')
        plot(log(err_P_translation), 'DisplayName', 'tr cam2')
        plot(log(err_rotation_P2), 'DisplayName', 'rot cam2')
        plot(log(err_C_translation), 'DisplayName', 'tr C')
        legend
        
        
    end
    [~, idx] = min(err_P_translation);
    min_errs = [err_radius(idx) err_vp(idx) err_P_translation(idx) err_rotation_P2(idx) err_C_translation(idx)]';
end
end

function [trueValues, lines] = setup_real(c1, c2, plots)
load('roller_coaster_two_views_im10_23.mat')

% Calibrate lines
ll1n = K'*jj; % ll1
ll2n = K'*ll; % ll2

lines.l_111 = pflat(ll1n(:,c1*2-1)); % Planes from cylinder nr 1 in camera 1
lines.l_112 = pflat(ll1n(:,c1*2));
lines.l_211 = pflat(ll2n(:,c1*2-1)); % Planes from cylinder nr 1 in camera 2
lines.l_212 = pflat(ll2n(:,c1*2));

lines.l_121 = pflat(ll1n(:,c2*2-1)); % Planes from cylinder nr 2 in camera 1
lines.l_122 = pflat(ll1n(:,c2*2));
lines.l_221 = pflat(ll2n(:,c2*2-1)); % Planes from cylinder nr 2 in camera 2
lines.l_222 = pflat(ll2n(:,c2*2));

trueValues.P1 = P1; % P1
trueValues.P2 = P0; % P2
trueValues.K = K;

if plots
    % Plot examples
    if sum([c1 c2] == [1 5]) == 2 || sum([c1 c2] == [2 3]) == 2
        figure(c1*100+c2*10+1)
        hold on
        imagesc(im1)
        rital(jj)
        axis equal
        figure(c1*100+c2*10+2)
        hold on
        imagesc(im0)
        rital(ll)
        axis equal
    end

end
end

function [trueValues, lines] = setup_synthetic(plots)
% Simulate cylinders and cameras
if plots
    figure(4); clf;
end 

% Simulate cylinders
cylinder1 = load_cylinder(plots, 1); % rotation in xz-plane, free translation
cylinder2 = load_cylinder(plots, 2); % rotation in xz-plane, free translation

% Simulate cameras.
% P = [R (-R*t)];

% Camera 1
R1 = moderatelyRandomRotationMatrix();
t1 = 3*randn(3,1);
cameras{1} = [R1 (-R1*t1)];

% Camera 2
R2 = moderatelyRandomRotationMatrix();
t2 = 3*randn(3,1);
cameras{2} = [R2 (-R2*t2)];

% Calculate lines

lines1 = lines_from_conic(cylinder1.iC1,cylinder1.iC2, cameras); % cylinder 1
lines2 = lines_from_conic(cylinder2.iC1,cylinder2.iC2, cameras); % cylinder 2

% Save lines
% Nomenclature: camera, cylinder, line_nr
lines.l_111 = pflat(lines1(:,1,1)); % Lines from cylinder nr 1 in camera 1
lines.l_112 = pflat(lines1(:,2,1));
lines.l_211 = pflat(lines1(:,1,2)); % Lines from cylinder nr 1 in camera 2
lines.l_212 = pflat(lines1(:,2,2));

lines.l_121 = pflat(lines2(:,1,1)); % Lines from cylinder nr 2 in camera 1
lines.l_122 = pflat(lines2(:,2,1));
lines.l_221 = pflat(lines2(:,1,2)); % Lines from cylinder nr 2 in camera 2
lines.l_222 = pflat(lines2(:,2,2));

% Save true values
trueValues.p_z = cylinder1.c(3);
trueValues.P1_R = R1;
trueValues.P1_t = t1;
trueValues.P2_t = t2;
trueValues.P2_R = R2;
trueValues.P1 = [R1 -R1*t1];
trueValues.P2 = [R2 -R2*t2];
trueValues.T1_t = cylinder1.translation;
trueValues.T1_R = cylinder1.R;
trueValues.T1_r = cylinder1.r;
trueValues.T2_t = cylinder2.translation;
trueValues.T2_R = cylinder2.R;
trueValues.T2_r = cylinder2.r;
trueValues.v1 = cylinder1.v;
trueValues.v2 = cylinder2.v;

if plots
    pts_1 = cylinder1.pts; % points on surface for plotting
    pts_2 = cylinder2.pts;
    
    trueValues.pts_1 = pts_1;
    trueValues.pts_2 = pts_2;
end

%% Plot cameras, conics and lines.
if plots
    
    % 3D plot
    figure(4); hold on;
    %Cameras
    pose = rigid3d(R1,t1');
    plotCamera('AbsolutePose',pose,'Opacity',0, 'AxesVisible', 1, 'Label', 'True P1');
    pose = rigid3d(R2,t2');
    plotCamera('AbsolutePose',pose,'Opacity',0, 'AxesVisible', 1, 'Label', 'True P2');
    hold off; axis equal; xlabel('x'); ylabel('y'); zlabel('z');
    
    % Project cylinders and lines into cam1
    figure(11); clf;
    rita(pflat(cameras{1}*pts_1),'.');% Project cylinder1 points into camera1
    hold on;
    rita(pflat(cameras{1}*pts_2),'.');% Project cylinder2 points into camera1
    title('Project original cylinders into camera1 (setup)')
    % draw lines
    axis manual;
    rital([lines.l_111 lines.l_112 lines.l_121 lines.l_122]);% lines.l_131 lines.l_132]); %camera, cylinder, line
    
    legend('cylinder1', 'cylinder2', ...
        'line1, c1', 'line2, c1', 'line1, c2', 'line2, c2')
    
    % Project cylinders and lines into cam2
    figure(12); clf;
    rita(pflat(cameras{2}*pts_1),'.');% Project cylinder1 points into camera2
    title('project original cylinders into camera2 (setup)')
    hold on;
    rita(pflat(cameras{2}*pts_2),'.');% Project cylinder2 points into camera2
    % Draw lines
    axis manual;
    rital([lines.l_211 lines.l_212 lines.l_221 lines.l_222]); % lines.l_231 lines.l_232]);
    
    legend('cylinder1', 'cylinder2', ...
        'line1, c1', 'line2, c1', 'line1, c2', 'line2, c2')
end

end

function R = moderatelyRandomRotationMatrix()
% Rotation around y
theta = rand*0.3;
R1 = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)];
% Rotation around z
theta = rand*0.3;
R2 = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];
% Rotation around x
theta = rand*0.3;
R3 = [1 0 0; 0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];

R = R1*R2*R3;
end

function err_vp = error_degrees_vp(v2, trueValues)
v_true = trueValues.v2; v_true = v_true(1:3);
v_est = v2;

v_true = v_true/norm(v_true);
v_est  = v_est(1:3)/norm(v_est(1:3));

if sign(v_est(1)) ~= sign(v_true(1))
    v_est = -v_est;
end
sc_prod = v_true'*v_est;
err_vp = abs(acosd(sc_prod));
end

function err_radius = error_radius(r2, trueValues)
err_radius = abs(trueValues.T2_r - r2);
end

function err_C_translation = error_translation_cylinders(T2_t, v2, trueValues)
tr_true = trueValues.T2_t;
tr_est = T2_t;
v = v2;

d_vec = tr_est-tr_true;
% project d_vec onto v:
d_proj = (d_vec'*v)/norm(v)^2*v;
err_vec = d_vec-d_proj;

err_C_translation = norm(err_vec);
end

function [err_P_translation] = error_translation_camera(P2, trueValues)

P2_t_true = pflat(null(trueValues.P2)); P2_t_true = P2_t_true(1:3);
P2_t_est = pflat(null(P2)); P2_t_est = P2_t_est(1:3);

err_P_translation = norm(P2_t_true - P2_t_est);
end

function [err_rotation_P2 ] = error_rotation_camera2(P2, trueValues)
R_true = trueValues.P2(1:3,1:3);
R_est  = P2(1:3,1:3);

err_rotation_P2 = abs(acosd((trace(R_est*R_true')-1)/2));
end

function pts = cylinder_to_pts(T_R, T_t, T_r, density, start, finish)

A = [T_R' T_t; 0 0 0 1];

% Plot cylinders
[th,y]=meshgrid(0:0.1:(2*pi),start:density:finish);

% Cylinder
X = T_r*cos(th);
Y = y;
Z = T_r*sin(th);
pts = A*[X(:) Y(:) Z(:) ones(size(Z(:)))]'; % Cylinder. jfr x = A*x_t


end

function cylinders = load_cylinder(plots, cyl_nr)

% Simulate 3D cylinders
r = 0.3*rand+0.5; % radius
dist = 10+rand*2; % distance to origo

c1 = (rand-1)*5; % x coordinate
c3 = sqrt(dist^2-c1^2); % z coordinate
c =[c1 0 c3]; %  x,y,z coordinates for center line intersection with y-plane

% Find cylinder matrix C and its dual iC without rotations
ir = 1/(r^2);
C1= [diag([ir 0 ir]) (-[c1*ir;0;c3*ir]); ...
    (-[c1*ir;0;c3*ir]') (c1^2*ir + c3^2*ir -1)];

% Fix so we can get the dual evethough C1 is not invertible
C2=diag([0 1 0 0]); % The plane z=0
iC1tmp = zeros(4,4);
iC1tmp([1 3 4],[1 3 4]) = inv(C1([1 3 4],[1 3 4]));
iC1 = iC1tmp;
iC2 = C2;

% Apply rotations
R = randrot([0;0;2.5]+ 0.1*randn(3,1));
T = [R zeros(3,1);zeros(1,3) 1];
%C1 = T'*C1*T; % The C-matrix containing the restrictions of cylinder-points.
%OBS! inv(T)=T' jfr T=inv(A) C = inv(A')*T*inv(A) <=> C1 '=' T'*C1*T
%C2 = T'*C2*T;
iC1 = T'*iC1*T; % Dual of C, containing the tangent planes
iC2 = T'*iC2*T;

[~,~,V]=svd(iC2);
NN=V(:,1);
tr = R'*c';
v = R(2,:)';

if plots
    % Code to extract points
    [th,y]=meshgrid(0:0.1:(2*pi),-10:1:10);
    X = r*cos(th) + c1; % Rotation och translation
    Y = y;
    Z = r*sin(th) + c3;
    A = cos(th)/r;
    B = zeros(size(y));
    C = sin(th)/r; % zeros(size(z));
    D = - ( A.*X + B.*Y + C.*Z );
    cylinders.pts = T'*[X(:) Y(:) Z(:) ones(size(Z(:)))]'; % Cylinder skrew. jfr x = A*x_t
    cylinders.planes = T'*([A(:) B(:) C(:) D(:)]'); % Planes for every point jfr pi = C*x
    cylinders.Tsave = T;
    
    figure(4); title('Original cameras and cylinders in 3D (setup)')
    hold on;
    plot3(cylinders.pts(1,:),cylinders.pts(2,:),cylinders.pts(3,:),'*','DisplayName',['cylinder ' num2str(cyl_nr)]); hold on;
    plot3(tr(1),tr(2),tr(3),'*k','DisplayName', ['center point cylinder  ' num2str(cyl_nr)]');
    plot3([tr(1) tr(1)+v(1)],[tr(2) tr(2)+v(2)],[tr(3) tr(3)+v(3)],'k', 'DisplayName', ['v/direction cylinder ' num2str(cyl_nr)]);
    legend
end

cylinders.c = c;
cylinders.NN = NN;
cylinders.iC1 = iC1;
cylinders.iC2 = iC2;
cylinders.r = r;
cylinders.R = R;
cylinders.translation = tr;
cylinders.v = v;
end

function lines = lines_from_conic(iC1, iC2, P)

n = length(P);
% Fix lines
lines = zeros(3,2,n);
%Project 3D conics to 2D conics
for i = 1:n
    ic1 = P{i}*iC1*P{i}';
    ic2 = P{i}*iC2*P{i}';
    [~,~,v]=svd(ic2);
    b1 = v(:,2);% b1 and b2 is a base for the nullspace of the vaishing point.
    b2 = v(:,3);
    cc = [b1'*ic1*b1 2*b1'*ic1*b2 b2'*ic1*b2]; % (b1'x + b2)C(b1'x + b2)'=0 for lines that are tangents
    ll = roots(cc);
    lines(:,1,i) = b1*ll(1)+b2;
    lines(:,2,i) = b1*ll(2)+b2;
end


end

function lnew = addnoise(l, noise, plots, figure_no)

if noise~=0
    % normal, direction and point on line
    N = l(1:2,:); V=[-l(2,:); l(1,:)]; P=[l(2,:); -l(1,:)-l(3,:)./l(2,:)];
    Q = (sum(N.*P)./sum(N.*N)).*N; % projektion of p on n
    
    % Add noise
    W1 = Q+1/2*V./sqrt(sum(V.*V)) + normrnd(0,noise,size(Q));
    W2 = Q-1/2*V./sqrt(sum(V.*V)) + normrnd(0,noise,size(Q));
    
    % Construct new line
    k = (W2(2,:)-W1(2,:))./(W2(1,:)-W1(1,:));
    
    lnew = [k; -ones(1,4); -k.*W1(1,:)+W1(2,:)];
    
    if plots
        figure(figure_no);
        hold on;
        rital(l);
        plot(Q(1,:),Q(2,:),'*k')
        plot(W1(1,:),W1(2,:),'*r')
        plot(W2(1,:),W2(2,:),'*g')
        rital(lnew, '--');
    end
else
    lnew = l;
end
end
