function R=randrot(direction);
% RANDROT - creates a rotation matrix that maps vectors in 'direction'
% towards the z-axis. The orientation around this axis is random.
% 'direction' is optional

if nargin==0,
  direction = 2*rand(3,3)-1;
end


[u,d,v]=svd(direction);

if det(u)<0,
 R=(u(:,[3 2 1]))';
else
 R=(u(:,[2 3 1]))';
end;

fi=rand*2*pi;
R=[cos(fi) sin(fi) 0;-sin(fi) cos(fi) 0;0 0 1]*R;
