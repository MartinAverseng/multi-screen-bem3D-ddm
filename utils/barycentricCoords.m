function[bary] = barycentricCoords(VH,Vh)

if size(VH,1) == 2 
    % 1D
    bary = zeros(size(Vh,1),2);
    A = VH(1,:);
    B = VH(2,:);
    AX = Vh - A;
    AB = B - A;
    bary(:,1) = 1 - AX.*AB/(AB.*AB);
    bary(:,2) = 1 - bary(:,1);
else
    A = VH(1,:);
    B = VH(2,:);
    C = VH(3,:);
    P = [A; B; C];
    T = [1 2 3];
    TR = triangulation(T,P);
    id = 0*Vh(:,1)+1;
    bary = cartesianToBarycentric(TR,id,Vh);
end

end