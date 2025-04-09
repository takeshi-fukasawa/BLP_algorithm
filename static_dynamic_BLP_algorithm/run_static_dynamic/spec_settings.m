    spec=spec_default;
    if method==1 % fixed point iteration
        spec.update_spec=0;
        if skip_contraction_spec==1
            spec.ITER_MAX=1;
       end
    elseif method==2|method==3 % spectral or SQUAREM
        if t_dependent_spectral_coef_spec==1
            spec.dim_hetero_spectral_coef=t_dim_id;
        else
            spec.update_spec=[];%%%%%%
        end
        if method==3 % SQUAREM
            spec.SQUAREM_spec=1;            %spec.SQUAREM_spec=3;%%%%
        end
    elseif method==4 % Anderson
        spec.Anderson_acceleration=1;
    end

