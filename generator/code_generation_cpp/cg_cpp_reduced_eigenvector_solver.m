function [code_setup_AA,code_extract_sol] = cg_cpp_reduced_eigenvector_solver(red_eig_solver, opt)
    
code_setup_AA = [];
indent = [opt.cg_indentation opt.cg_indentation];
for k = 1:red_eig_solver.AA_sz(2)
    ind = red_eig_solver.ind_AA1{k};
    
    str = sprintf('%sAA.col(%d) = ',indent,k-1);
    for i = 1:size(ind,2);
        if ind(2,i) == 0            
            str = [str sprintf('AMs.col(%d) + ',ind(1,i)-1)];
        else
            str = [str sprintf('zi[%d] * AMs.col(%d) + ',ind(2,i)-1,ind(1,i)-1)];  
        end
    end
    str = str(1:end-3);
    
    code_setup_AA = [code_setup_AA str sprintf(';\n')];
end

for k = 1:red_eig_solver.AA_sz(2)
    ind = red_eig_solver.ind_AA2(:,k);
    
    str = sprintf('%sAA(%d,%d) = AA(%d,%d) - zi[%d];\n',indent,ind(1)-1,k-1,ind(1)-1,k-1,ind(2)-1);    
    code_setup_AA = [code_setup_AA str];
end

code_extract_sol = [];

for k = 1:length(red_eig_solver.ind_var)
    if red_eig_solver.ind_var(k) > 0
        code_extract_sol = [code_extract_sol sprintf('%ssols(%d,i) = s(%d);\n',indent,k-1,red_eig_solver.ind_var(k)-1)];
    else
        code_extract_sol = [code_extract_sol sprintf('%ssols(%d,i) = zi[0];\n',indent,k-1)];
    end
    
end

%     fprintf('AA(:,%d) = %s;\n',k,str(1:end-1));

end

