function T = pcc_trans(theta,phi,l)

arguments (Input)
    theta   (1,1) {mustBeReal}
    phi     (1,1) {mustBeReal}
    l       (1,1) {mustBeReal}
end
arguments (Output)
    T       (3,1) {mustBeReal}
end

sinc_mma = @(x) sinc(x/sym(pi));

T = l * [ cos(phi) * sinc_mma(theta/2) * sin(theta/2) ;
          sin(phi) * sinc_mma(theta/2) * sin(theta/2) ;
                          sinc_mma(theta)               ];

end