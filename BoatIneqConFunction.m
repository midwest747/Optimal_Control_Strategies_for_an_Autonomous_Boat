function ceq = BoatIneqConFunction(X,U,e,data)
p = data.PredictionHorizon;
X1 = X(2:p+1,2);
river_length = 10;
ceq = [X1 - (river_length+.01);-X1-.001];

