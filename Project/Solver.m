iter = 0;
residue = 10*min_residue;
while iter<max_iter && residue>min_residue
    for k =2:Nz-1
        dpdx(k) = 1;
    end
    residue = 1;
    iter = iter+1
end