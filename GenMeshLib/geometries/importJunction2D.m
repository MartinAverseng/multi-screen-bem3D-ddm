mOmegaJunction = msh('junction2DOmega.msh');
mGammaJunction = msh('junction2DGamma.msh');

plot(mOmegaJunction);
hold on
plot(mGammaJunction,'r');


Mjunction = fracturedMesh(mOmegaJunction,mGammaJunction);
save('mOmegaJunction','mOmegaJunction');
save('mGammaJunction','mGammaJunction');