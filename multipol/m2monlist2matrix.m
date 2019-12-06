function mm = m2monlist2matrix(ltgb,maxdeg);
% mm = m2monlist2matrix(ltgb,maxdeg);
%

if nargin<2
    maxdeg = 7;
end
ltgb = [ltgb ' '];
aa = find(ltgb==' ');
mm = zeros(maxdeg,length(aa));
off = 0;
for kk = 1:length(aa);
    spacepos = aa(kk);
    onemon = ltgb((off+1):(spacepos-1));
    onemm = zeros(maxdeg,1);
    while length(onemon)>0,
        bb = find(onemon=='x');
        if length(bb)>0,
            bb = bb(end);
            onepart = onemon(bb:end);
            powerpos = find(onepart=='^');
            if length(powerpos)==1,
                pow = str2num(onepart((powerpos+1):end));
                onepart = onepart(1:(powerpos-1));
            elseif length(powerpos)==0
                pow = 1;
            else
                error('error in powerpos');
            end
            ii = str2num(onepart(2:end));
            if ii>maxdeg,
                mm( (maxdeg+1):ii,:)=zeros( ii-maxdeg,size(mm,2));
                maxdeg=ii;
            end;
            onemm(ii)=pow;
            onemon = onemon(1:(bb-1));
        end;
    end;
    mm(:,kk)=onemm;
    off = spacepos;
end
% hack to remove all zero columns
badid = find(all(mm==0));
mm(:,badid)=[];

