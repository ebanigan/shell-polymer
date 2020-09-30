
void lamin_pol_binding(int polnum1, int mono1, int mono2)
{

//fprintf(stderr, "lpb for %i %i\n", mono1, mono2);
    int ii, kk;
    int polnum2;
    double dotproduct = 0.;
    double rotation1[DIMENSION];
    double rotation2[DIMENSION];
    double dotproduct_contrib;
    
    
    for(ii = 0; ii < mono_list[mono2].get_lamin_pol_used(); ii++)
    {
        if(mono_list[mono2].get_lamin_pol_index(ii) == mono1)
        {
            polnum2 = ii;
            break;
        }
    }

    for(kk = 0; kk < DIMENSION; kk++)
    {
        dotproduct += mono_list[mono1].get_prev_lamin_pol(polnum1, kk)*mono_list[mono2].get_prev_lamin_pol(polnum2, kk);
    }
    
    if(dotproduct < 0.)
    {
    dotproduct_contrib = integer_power(dotproduct, LAMIN_ROTATION_EXPONENT-1);

    for(kk = 0; kk < DIMENSION; kk++)
    {
        rotation1[kk] = -2.*LAMIN_ROTATION_EXPONENT*LAMIN_ROTATION_EPSILON*dotproduct_contrib*(mono_list[mono2].get_prev_lamin_pol(polnum2, kk) - dotproduct*mono_list[mono1].get_prev_lamin_pol(polnum1, kk))*mono_list[mono1].get_rdiffusion_coeff()*dt;
        rotation2[kk] = -2.*LAMIN_ROTATION_EXPONENT*LAMIN_ROTATION_EPSILON*dotproduct_contrib*(mono_list[mono1].get_prev_lamin_pol(polnum1, kk) - dotproduct*mono_list[mono2].get_prev_lamin_pol(polnum2, kk))*mono_list[mono2].get_rdiffusion_coeff()*dt;
    }
    
    mono_list[mono1].rotate(rotation1);
    mono_list[mono2].rotate(rotation2);
    }//if dotproduct < 0.
//	if dotproduct > 0., do nothing apparently

}

