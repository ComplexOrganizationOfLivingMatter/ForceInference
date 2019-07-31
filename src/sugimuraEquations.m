%Equations Sugimura - 2013
%The mechanical anisotropy in a tissue promotes ordering in hexagonal cell packing
% Kaoru Sugimura and Shuji Ishihara

%% Parameters notation
% N_uv      -> Global stress tensor
% Delta_uv  -> Kronecker's delta
% pressure (P) -> Pressure of the cell
% tension (T)         -> Tension of the edge
% area (A)  -> Area of the cell
% edgeLength (l)         -> Edge length
% lx        -> Component X of edge length l
% ly        -> Component Y of edge length l

%Calculate area of cells
area = calculateAreaOfCells();
pressure = getPressureCells();
tension = getTensionCells();


%Vertices and edges
[lx, ly, edgeLength] = calculateVertices();


%% Potencial enerty   -> U
U = U_0 - lambdaX * Tx - lambdaY * Ty;

% first term
U_0 = U_ar + U_lin + U_cor;

U_ar = (K / 2) * (area - A_0)^2; % elasticity of the cell
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
Txx = sum(K * ( area - A_0) * area) + sum( sigma + upperLambda * (Li + Lj)) * lx*lx/abs(l);
Tyy = sum(K * ( area - A_0) * area) + sum( sigma + upperLambda * (Li + Lj)) * ly*ly/abs(l); %(expecificity T^)

%at final stady state, Fi=0 for all i and T = T^