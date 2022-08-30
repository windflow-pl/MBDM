function scew = scew_sym(vec)
    scew(1,1) = 0.0;       scew(1,2) = -vec(3,1);   scew(1,3) = vec(2,1);
    scew(2,1) = vec(3,1);    scew(2,2) = 0.0;       scew(2,3) = -vec(1,1);
    scew(3,1) = -vec(2,1);   scew(3,2) = vec(1,1);    scew(3,3) = 0.0;
end