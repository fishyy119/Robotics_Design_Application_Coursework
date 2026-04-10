if 1
nDegree=6;%”–∂‡…ŸŒ¨

fid=fopen(   'rhs.m','w');
fprintf(fid, 'function zdot=rhs(t,z,flag,m,xG,yG,zG,I1,g)   \n\n');
for ii=1:nDegree
    fprintf(fid, 'q%d = z(%d);   u%d = z(%d);\n',ii,2*ii-1,ii,2*ii);
end
fprintf(fid, 'I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3);I2_1 = I2(1,1);I2_2 = I2(2,2);I2_3 = I2(3,3); \n\n');

for ii=1:nDegree
    for jj=1:nDegree
        fprintf(fid,'M%d%d =;', ii,jj );
    end
    fprintf(fid,'\n\n');
end

fprintf(fid,'MM = [');
for ii=1:nDegree
    for jj=1:nDegree
        fprintf(fid,'M%d%d ',ii,jj);
    end
    if ii~=nDegree
        fprintf(fid,';\n');
    end
end
fprintf(fid,' ];\n\n');

fprintf(fid,'RHS = [');
for ii=1:nDegree
    fprintf(fid,'RHS%d',ii);
    if ii~=nDegree
        fprintf(fid,';');
    end
end
fprintf(fid,' ];\n\n');

fprintf(fid,'X = MM \\ RHS;                                    \n\n');

for ii=1:nDegree 
    fprintf(fid,'ud%d = X(%d); \n',ii,ii);
end

fprintf(fid,'zdot = [');
for ii=1:nDegree 
    fprintf(fid,'u%d ud%d ',ii,ii); 
end
fprintf(fid,' ]'';\n\n');
fclose(fid);
% 
% fid=fopen(   'energy.m','w');
% fprintf(fid, 'function [KE, PE]=energy(t,z,flag,m,xG,yG,zG,I1,g)   \n\n');
%  
% for ii=1:nDegree
%     fprintf(fid, 'q%d = z(%d);   q%d = z(%d);\n',2*ii-1,2*ii-1,2*ii,2*ii);
% end
% fprintf(fid, 'I1_1 = I1(1,1);I1_2 = I1(2,2);I1_3 = I1(3,3);I2_1 = I2(1,1);I2_2 = I2(2,2);I2_3 = I2(3,3); \n');
% fprintf(fid, 'I1_12 = I1(1,2);I1_13 = I1(1,3);I1_23 = I1(2,3);I2_12 = I2(1,2);I2_13 = I2(1,3);I2_23 = I2(2,3); \n');
% fprintf(fid,'KE = %s;\n',char(KE)); 
% fprintf(fid,'PE = %s;\n',char(PE)); 
% fclose(fid);
end