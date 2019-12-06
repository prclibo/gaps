function res = runM2(eq,m2path)
%runM2 Generate M2-code, run it and return results to matlab
%[m2dim,m2degree,m2basis] = runM2(eq)

if nargin<2
    m2path = 'M2';
end
if exist('testxx.m2','file')
    delete('testxx.m2');
end
if exist('testxx.m2.out','file')
    delete('testxx.m2.out');
end
fname = 'testxx.m2';
eqs2m2(eq,fname);
eval(['! ' m2path ' testxx.m2']);
%eval(['!m2 testxx.m2.out']);
while ~exist('testxx.m2.out','file')
    pause(1);
end
fid=fopen('testxx.m2.out','r');
l1 = fgetl(fid);
l2 = fgetl(fid);
l3 = fgetl(fid);
fclose(fid);
res.m2dim = str2double(l1);
res.m2degree = str2double(l2);
res.m2basis = l3;
if exist('testxx.m2','file')
    delete('testxx.m2');
end
if exist('testxx.m2.out','file')
    delete('testxx.m2.out');
end

