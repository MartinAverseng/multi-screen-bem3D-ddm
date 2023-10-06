function [GV,mult] = genVertOfElt(M,l,F,gamma)


GV = [];
mult = [];


for i = 1:length(l)
    el = l(i);
    V = M.elt(el,:);
    for k = 1:M.n+1
        Vk = V(:,k);
        gvs = find(F==Vk);
        mult(end+1) = length(gvs);
        for j = 1:length(gvs)
            if ismember(el,gamma{gvs(j)})
                GV(end+1) = gvs(j);
            end
        end
    end
end

[GV,I] = unique(GV);
mult = mult(I);


end

