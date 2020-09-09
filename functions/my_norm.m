function v_out = my_norm(v_in)
minv = min(v_in(:));
maxv = max(v_in(:));
v_out = (v_in-minv) / (maxv-minv);
end