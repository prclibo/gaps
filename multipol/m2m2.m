function row = m2m2(M);


[m,n]=size(M);
row = 'matrix{';
for i=1:m;
 row = [row '{'];
 for j=1:n-1,
  row = [row num2str(M(i,j),8) '*1_R,'];
 end
 if i<m,
  row = [row num2str(M(i,n),8) '*1_R},'];
 else
  row = [row num2str(M(i,n),8) '*1_R}}'];
 end;
end
disp(row);
