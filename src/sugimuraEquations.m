%Equations Sugimura - 2013

% N_uv      -> Stress
% Delta_uv  -> Kronecker's delta
% P         -> Pressure
% T         -> Tension???
% A         -> Area
% l         -> Edge length
% lx        -> Component X of l
% ly        -> Component Y of l
% Cell stress:
N_uv = ((-A * P * Delta_uv ) +  T * ( lx * ly ) / l ) / A ;

%% Potencial enerty   -> U
U = U_0 - lambdaX * Tx - lambdaY * Ty;

% first term
U_0 = U_ar + U_lin + U_cor;

U_ar = (K / 2) * (A - A_0)^2; % elasticity of the cell
U_lin = sigma * l; %line tension = sum of cell adhesion and contracting force
U_cor = (upperLambda/ 2)*L^2; % cortical elasticity

% second and third terms add the tissue strecth.
% lambdaX and lambdaY    -> variables parameterizing the system size.
% Tx and Ty              -> stresses applied along the x and y axes.
Lx = lambdaX * Lx_0;
Ly = lambdaY * Ly_0;

% T is derived from the follow differential equation:
dxdt = -lambdaX^2 * Fx + ((Tx - Txx) / lambdaX) * x;
dydt = -lambdaY^2 * Fy + ((Ty - Tyy) / lambdaY) * y;
dLambdaXdt = Tx - Txx;
dLambdaYdt = Ty - Tyy;

F = dU_0/dX;
Txx = sum(K * ( A - A_0) * A) + sum( sigma + upperLambda * (Li + Lj)) * lx*lx/abs(l);
Tyy = sum(K * ( A - A_0) * A) + sum( sigma + upperLambda * (Li + Lj)) * ly*ly/abs(l); %(expecificity T^)

%at final stady state, Fi=0 for all i and T = T^