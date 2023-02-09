function res = meniscusfit(xinput, yinput, coefs)
size_x = size(xinput);
size_y = size(yinput);
if (~all(size_x == size_y))
    error('vector size does not match');
end
len_coefs = length(coefs);
%%
% % f(x) = a + exp(-x/l)(1/x)^(1/2)(poly(x))
l_surf = sqrt(0.072/1e3/9.8)*1e3;
part_polynomial = polyval(coefs(2:len_coefs),xinput);
% l_surf = coefs(2);
% part_polynomial = polyval(coefs(3:len_coefs),xinput);
part_exp = exp(-xinput/l_surf).*sqrt(1./xinput);
fx = coefs(1) + part_exp.*part_polynomial;

res = fx-yinput;
end
