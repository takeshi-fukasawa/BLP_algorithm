    spec=spec_default;
    if method==1 % fixed point iteration
        spec.update_spec=0;
        if skip_contraction_spec==1
            spec.ITER_MAX=1;
       end
    elseif method==2|method==3 % spectral
        if t_dependent_alpha_spec==1
            spec.update_spec=t_dim_id;
        else
            spec.update_spec=[];%%%%%%
        end
        if method==3
            spec.SQUAREM_spec=1;
        end
    end

