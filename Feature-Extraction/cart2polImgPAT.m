function I_polar = cart2polImgPAT(V_rf_sum_smooth,x,z)
[N_z,~] = size(V_rf_sum_smooth);
[X,Z] = meshgrid(x,z);
[~,Rho_query] = cart2pol(X,Z);
rhoRange = linspace(0, max(Rho_query(:)), N_z);
thetaRange = linspace(-pi/2, pi/2, N_z);
[T, R] = meshgrid(thetaRange, rhoRange);
[Z_query, X_query] = pol2cart(T,R);
I_polar = interp2(X, Z, V_rf_sum_smooth, X_query, Z_query,'cubic');
I_polar(isnan(I_polar))=0;
%imagesc(I_polar)
end
