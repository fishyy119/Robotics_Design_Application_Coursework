clc
clear all

syms rx ry rz wx wy wz real
r=[rx ry rz]'
w=[wx wy wz]'
temp = cross(cross(r,w),r)
temp1 = cross(r,cross(w,r))
temp - temp1
collect(temp, [wx wy wz]')