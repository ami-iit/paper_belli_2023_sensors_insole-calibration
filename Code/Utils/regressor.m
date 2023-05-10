function Phi = regressor(C, POLYNOMIAL_ORDER, HISTORY_SAMPLES, HISTORY_POLYNOMIAL_ORDER)

% C is m x 1 --> capacitance of one taxel, m are samples in the time
%                or one sample for all taxels

Phi = [];

for ii = POLYNOMIAL_ORDER : (-1) : 0
  Phi = [Phi   C.^ii];
end

C_temp = C;

for kk = 1 : HISTORY_SAMPLES
  C_temp = [C(1);  C_temp];
  
  C_temp(end) = [];

  for ii = 1 :  HISTORY_POLYNOMIAL_ORDER
    Phi = [Phi   C_temp.^ii];
  end
end

end