function [] = highlight(M,gVs,s)

[F,~,~] = M.generalizedSubfacets(0);
V = F(gVs);
plot3(M.vtx(V,1),M.vtx(V,2),M.vtx(V,3),s);



end

