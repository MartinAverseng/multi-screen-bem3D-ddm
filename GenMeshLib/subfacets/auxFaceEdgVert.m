function [marked] = auxFaceEdgVert(M,F,gamma,GenV,A,marked)

marked(GenV) = 1;
l = neighbors(M,F,gamma,GenV);
l = l(A(F(l))==A(F(GenV)));
for i = 1:length(l)
    if ~marked(l(i))
        marked = auxFaceEdgVert(M,F,gamma,l(i),A,marked);
    end
end


end




