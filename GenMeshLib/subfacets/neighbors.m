function[GenVs] = neighbors(M,F,gamma,genV)


GenVs = [];
l = gamma{genV};
for el = l
    V = setdiff(M.elt(el,:),F(genV));
    for i = 1:length(V)
        Vi = V(i);
        ind = find(F==Vi);
        for j = 1:length(ind)
            genVj = ind(j);
            lj = gamma{genVj};
            if ismember(el,lj)
                GenVs = [GenVs,genVj];
            end
        end
        
    end
end
GenVs = unique(GenVs);





end