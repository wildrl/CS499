#max(ic2_mag113$magnitude))

#IC2
mag_plot(ic2_mag11$magnitude, ic2_mag24$magnitude, ic2_mag53$magnitude,
         ic2_mag64$magnitude, ic2_mag113$magnitude, "ic2=(120,0,120,0) L2 Norm", 
         nrow(ic2_mag11), max(ic2_mag64$magnitude))
dot_plot(ic2_mag11$dot_product, ic2_mag24$dot_product, ic2_mag53$dot_product,
         ic2_mag64$dot_product, "ic2=(120,0,120,0) dot product with 113-bit soln", 
         nrow(ic2_mag11), max(ic2_mag64$dot_product))
r_error_plot(ic2_mag11$r_error, ic2_mag24$r_error, ic2_mag53$r_error, ic2_mag64$r_error,
             ic2=(120,0,120,0) relative error from 113-bit soln,
             nrow(ic2_mag11), max(ic2_mag64$r_error))

#IC3
mag_plot(ic3_mag11$magnitude, ic3_mag24$magnitude, ic3_mag53$magnitude,
         ic3_mag64$magnitude, ic3_mag113$magnitude, "ic3=(121,0,120,0) L2 Norm", 
         nrow(ic3_mag11), max(ic3_mag64$magnitude))
dot_plot(ic3_mag11$dot_product, ic3_mag24$dot_product, ic3_mag53$dot_product,
         ic3_mag64$dot_product, "ic3=(121,0,120,0) dot product with 113-bit soln", 
         nrow(ic3_mag11), max(ic3_mag64$dot_product))
r_error_plot(ic3_mag11$r_error, ic3_mag24$r_error, ic3_mag53$r_error, ic3_mag64$r_error,
             ic3=(121,0,120,0) relative error from 113-bit soln,
             nrow(ic3_mag11), max(ic3_mag64$r_error))

#IC4
mag_plot(ic4_mag11$magnitude, ic4_mag24$magnitude, ic4_mag53$magnitude,
         ic4_mag64$magnitude, ic4_mag113$magnitude, "ic4=(120,1,120,0) L2 Norm", 
         nrow(ic4_mag11), max(ic4_mag64$magnitude))
dot_plot(ic4_mag11$dot_product, ic4_mag24$dot_product, ic4_mag53$dot_product,
         ic4_mag64$dot_product, "ic4=(120,1,120,0) dot product with 113-bit soln", 
         nrow(ic4_mag11), max(ic4_mag64$dot_product))
r_error_plot(ic4_mag11$r_error, ic4_mag24$r_error, ic4_mag53$r_error, ic4_mag64$r_error,
             ic4=(120,1,120,0) relative error from 113-bit soln,
             nrow(ic4_mag11), max(ic4_mag64$r_error))



mag_plot113(ic2_mag113$magnitude, ic3_mag113$magnitude, ic4_mag113$magnitude,
            ic5_mag113$magnitude, ic6_mag113$magnitude,
            "L2-norm of f(vi, 113)", 100000,20,
            "f(v0, 113)", "f(v1, 113)", "f(v2, 113)", "f(v3, 113)", "f(v4, 113)")
mag_plot113(ic7_mag113$magnitude, ic8_mag113$magnitude, ic9_mag113$magnitude,
            ic10_mag113$magnitude, ic11_mag113$magnitude,
            "L2-norm of f(ui, 113)", 100000,10,
            "f(u0, 113)", "f(u1, 113)", "f(u2, 113)", "f(u3, 113)", "f(u4, 113)")

mag_plot113(ic2_mag64$magnitude, ic3_mag64$magnitude, ic4_mag64$magnitude,
            ic5_mag64$magnitude, ic6_mag64$magnitude,
            "L2-norm of f(vi, 64)", 200000,20,
            "f(v0, 64)", "f(v1, 64)", "f(v2, 64)", "f(v3, 64)", "f(v4, 64)")
mag_plot113(ic2_mag53$magnitude, ic3_mag53$magnitude, ic4_mag53$magnitude,
            ic5_mag53$magnitude, ic6_mag53$magnitude,
            "L2-norm of f(vi, 53)", 200000,20,
            "f(v0, 53)", "f(v1, 53)", "f(v2, 53)", "f(v3, 53)", "f(v4, 53)")
mag_plot113(ic2_mag24$magnitude, ic3_mag24$magnitude, ic4_mag24$magnitude,
            ic5_mag24$magnitude, ic6_mag24$magnitude,
            "L2-norm of f(vi, 24)", 200000,20,
            "f(v0, 24)", "f(v1, 24)", "f(v2, 24)", "f(v3, 24)", "f(v4, 24)")
mag_plot113(ic2_mag11$magnitude, ic3_mag11$magnitude, ic4_mag11$magnitude,
            ic5_mag11$magnitude, ic6_mag11$magnitude,
            "L2-norm of f(vi, 11)", 200000,20,
            "f(v0, 11)", "f(v1, 11)", "f(v2, 11)", "f(v3, 11)", "f(v4, 11)")


fx_plot(ic3_pol11$th1, ic3_pol11$w1, ic3_pol11$th2, ic3_pol11$w2, ic2_mag11$magnitude,
        "ic2 11 bit function values", 10000, 40)

fx_plot(ic3_pol113$th1, ic3_pol113$w1, ic3_pol113$th2, ic3_pol113$w2, ic2_mag113$magnitude,
        "ic2 11 bit function values", 10000, 40)

