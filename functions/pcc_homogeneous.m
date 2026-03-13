function H = pcc_homogeneous(theta,phi,l)

arguments (Input)
    theta   (1,1) {mustBeReal}
    phi     (1,1) {mustBeReal}
    l       (1,1) {mustBeReal}
end
arguments (Output)
    H       (4,4) {mustBeReal}
end

H = [pcc_rotmat(theta,phi), pcc_trans(theta,phi,l);
     zeros(1,3), 1];
end