function R  = Rot(axis, q )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% syms axis
if axis == 'z'
    R=[cos(q) -sin(q) 0;sin(q) cos(q) 0;0 0 1];
end
if axis == 'y'
    R=[cos(q) 0 sin(q);0 1 0;-sin(q) 0 cos(q)];
end
if axis == 'x'
    R=[1 0 0; 0 cos(q) -sin(q);0 sin(q) cos(q)];
end
end

