function [Klist] = incidentElements(gamma,dofList)


Klist = [];

for ind = 1:length(dofList)
    genF = dofList(ind);
    Klist = union(Klist,gamma{genF});
end


end

