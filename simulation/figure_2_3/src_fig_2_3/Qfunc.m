function opt = Qfunc(dX,dE,b)
% Objective function for CM estimators
% dX: difference in observable characteristics X between period s and t
% dE: estimated H(E(dy_s-dy_t|X))
% b: beta that is to be estimated (D by 1 vector)
% criterion:
% minsum[i,j,t]{1{E[y_ijt-y_ijs|Xi]>0}1{X_ijt'b<=X_ijs'b}1{X_ikt'b>=X_iks'b}}
% can try strict inequality in above equations

[N,D,J,T0] = size(dX);
%dE1 = 2 * normcdf(subplus(dE)) - 1;

indb = permute(sum(dX .* repmat(b', [N,1,J,T0]), 2),[1,3,4,2]); % N by J by T0
dbm = zeros(N,J,T0);

%dbm(:,j,t) = (indb(:,j,t) <= 0).* prod(double(indb(:,1:end ~=j,t) >= 0),2);
dbm = dE .* double(indb <= 0);

%opt = sum(sum(sum(dE1 .* dbm)));
opt = sum(sum(2*normcdf(subplus(sum(dbm,2))) - 1));

end