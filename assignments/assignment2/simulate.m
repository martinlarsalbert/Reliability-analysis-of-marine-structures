function result = simulate(C, MEAN, N)

% N : Number of simulations (Number of simulated voyages)
n_seastates = length(MEAN);


%%

X=normrnd(0,1,n_seastates,N);
mean_values = MEAN*ones(1,N);

Z0=real(C'*X);  % (Z return as complex so taking the real part, unsure about this one...)
Z=Z0+mean_values;
Hs = exp(Z);

% Compute the fatigue damages of all simulated sea states by Eq.(3.2).
D=(0.27*Hs.^(2.5) + 0.35*Hs.^2)/11.76e7;
D_voy=sum(D);

%%
result = struct();
result.n_seastates = n_seastates;
result.C = C;
result.Z=Z;
result.Hs=Hs;
result.D=D;
result.D_voy=D_voy;