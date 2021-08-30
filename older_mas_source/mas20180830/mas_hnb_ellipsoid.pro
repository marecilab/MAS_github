function hnb_ellipsoid_get_gradients

    common scan_data, project
    
    ci = project.ci
    
    bmatrix = *project.imndarray[ci].b_matrix

    nbvals = n_elements(bmatrix)/6
    
    gradients = fltarr(4,nbvals)

    for b = 0, nbvals-1 do begin
    
        tmp_bmat = [ [ bmatrix[0,b], bmatrix[3,b]/2, bmatrix[4,b]/2 ], $
                     [ bmatrix[3,b]/2, bmatrix[1,b], bmatrix[5,b]/2 ], $
                     [ bmatrix[4,b]/2, bmatrix[5,b]/2, bmatrix[2,b] ] ]
                     
        evals = eigenql(tmp_bmat, eigenvectors=evecs)
        
        gradients[0,b] = total(evals)
        gradients[1:*,b] = evecs[*,0]

    endfor

    return, gradients
    
end

pro hnb_ellipsoid_curvefit, x, y, z

    common scan_data, project

    old_device = !D.name
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        set_plot, old_device
        loadct, 0
        print, !error_state.MSG
        help, /traceback
        void = dialog_message(['An error occurred. Please check to make sure that the', $
                               'Data set loaded has multiple b-values and multiple', $
                               'gradient directions.'], /error, /center)
        return
    endif
    
    ci = project.ci
    gr = hnb_ellipsoid_get_gradients()
    go = 1
        
    !P.multi = [ 0,2,2 ]
    ;;window, 5, xsize=1200, ysize=900
    
    if (n_elements(x) eq 0) then x = project.procpramarray[ci].fdim_start
    if (n_elements(y) eq 0) then y = project.procpramarray[ci].pdim_start
    if (n_elements(z) eq 0) then z = project.procpramarray[ci].sdim_start
        
    sugg_filename = string([x,y,z], format='(%"HnB_%03dX_%03dY_%03dZ_report.ps")')
    
    report = dialog_pickfile(file=sugg_filename, title='Choose a destination for the report:')
    if (report ne '') then begin
        set_plot, 'PS', /copy
        loadct, 6
        device, /landscape, /color, BITS_PER_PIXEL = 8, FONT_SIZE=12, language_level=2
        device, FILENAME=report, /inches, xsiz=11, ysiz=8.5, xoffset=0.25, yoffset=11-(0.25), scale_factor=0.95, /TIMES, /BOLD
        ;;device, decomposed=1
    endif
    
    data   = reform((*project.dataarray[ci].state1)[x,y,z,*])
    min_b  = min(gr[0,*], max=max_b)
    b_axis = (findgen(500)/(500-1)*(max_b-min_b) + min_b)/1000.
    yplot_max = max(data)

    while (go eq 1) do begin
    
        g = gr[1:3, 0]
        b = gr[0, 0]
        
        ;; find the common directions
        sim_ind = where((abs(gr[1,*]-g[0]) lt 0.2) AND $
                        (abs(gr[2,*]-g[1]) lt 0.2) AND $
                        (abs(gr[3,*]-g[2]) lt 0.2), n_found, $
                        complement=other, ncomplement=n_other)
        
        ;; get the bvals and grad and data for the first direction
        cf_grad = gr[1:3,sim_ind]
        cf_bval = gr[0  ,sim_ind]/1000.0
        measurements = reform(data[sim_ind])

        ;; curvefit to get alpha and D
        curvefit_adc_alpha, measurements, cf_bval, A, chisquare, sigma
        yfit_str = A[0]*exp( -((cf_bval)*A[1])^A[2] )
        resid_str = yfit_str - measurements
        yerror_str = sqrt(total(resid_str^2)/3.0)
        fitline_str = A[0]*exp( -((b_axis)*A[1])^A[2] )

        curvefit_adc, measurements, cf_bval, A_o, chisquare, sigma
        A_o = [A_o, 1]
        yfit_sin = A_o[0]*exp( -((cf_bval)*A_o[1])^A_o[2] )
        resid_sin = yfit_sin - measurements
        yerror_sin = sqrt(total(resid_sin^2)/3.0) 
        fitline_sin = A_o[0]*exp( -((b_axis)*A_o[1])^A_o[2] )
        if (report ne '') then begin
            plot, cf_bval, measurements, ystyle=1, xstyle=1, yrange=[0, yplot_max], psym=1, $
                  xtitle='bvalue (*1e-3)', ytitle='Intensity', title=string([x,y,z], format='(%"Stretch Exp Fit, Voxel=[%03d,%03d,%03d]")')
            oplot, b_axis, fitline_str, color=75
            oplot, b_axis, fitline_sin, color=140;'bcbcbc'x
            oplot, b_axis, fitline_str+yerror_str, linestyle=1, color=75;'ff0000'x
            oplot, b_axis, fitline_str-yerror_str, linestyle=1, color=75;'ff0000'x
            
            xyouts, mean(cf_bval), yplot_max-(1*0.05*yplot_max), /data, string(A[2], format='(%"alpha=%0.4f")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(2*0.05*yplot_max), /data, string(A[1]*1e-3, format='(%"D=%0.5f (mm^2/s)")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(3*0.05*yplot_max), /data, string(A[0], format='(%"S0=%0.3f")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(4*0.05*yplot_max), /data, string(g, format='(%"G=[%0.3f,%0.3f,%0.3f]")'), charsize=0.8                                                                    

            plot, cf_bval, measurements, ystyle=1, xstyle=1, yrange=[0, yplot_max], psym=1, $
                  xtitle='bvalue (*1e-3)', ytitle='Intensity', title=string([x,y,z], format='(%"Single Exp Fit, Voxel=[%03d,%03d,%03d,]")')
            oplot, b_axis, fitline_sin, color=140
            oplot, b_axis, fitline_str, color=75;'bcbcbc'x
            oplot, b_axis, fitline_sin+yerror_sin, linestyle=1, color=140;'ff0000'x
            oplot, b_axis, fitline_sin-yerror_sin, linestyle=1, color=140;'ff0000'x
            
            ;;xyouts, mean(cf_bval), yplot_max-(1*0.1*yplot_max), /data, string(A_o[2], format='(%"alpha=%0.4f")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(2*0.05*yplot_max), /data, string(A_o[1]*1e-3, format='(%"D=%0.5f (mm^2/s)")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(3*0.05*yplot_max), /data, string(A_o[0], format='(%"S0=%0.3f")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(4*0.05*yplot_max), /data, string(g, format='(%"G=[%0.3f,%0.3f,%0.3f]")'), charsize=0.8                        
        endif
        print, A[1]
        ;; scale the gradient vector by the result and add it to an array
        if (n_elements(alpha) eq 0) then begin
            diff =  A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ]
            diff = [ diff, -A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
            
            alpha =  A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ]
            alpha = [ alpha, -A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
        endif else begin
            diff = [ diff, A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
            diff = [ diff, -A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
            
            alpha = [ alpha, A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
            alpha = [ alpha, -A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
        endelse
        
        if (n_other eq 0) then begin
            go = 0
            break
        endif
        
        gr = gr[*,other]
        data = data[other]
                 
    endwhile

    if (report ne '') then begin
        device, /close_file
        set_plot, old_device
        loadct, 0
    endif

    alpha = reform(alpha, 3, n_elements(alpha)/3)
    diff  = reform(diff, 3, n_elements(diff)/3)/max(diff)*0.9
    
    if (n_elements(diff)/3 lt 4) then begin
        void = dialog_message(['Not enough gradient measurements to create ellipsoids.'], /error, /center)
        return
    endif
    
    ;; create the best-fit ellipse to the scaled gradients
    hnb_ellipsoid_pointset, alpha, title="Alpha"
    
;    qhull, diff, tr,  /delaunay
;    poly = tetra_surface(diff, tr)
;    o1 = obj_new('idlgrpolygon', diff, color=[128,128,128], polygon=poly, alpha=0.5, reject=1)
;    o2 = obj_new('idlgrpolygon', diff, style=0, thick=5, color=[128,0,0])
;    an0 = obj_new('idlgrpolyline', [[0,0,0],[0,0,1], [0,0,0], [0,0,-1]], color=[128,128,128], linestyle=2)
;    an1 = obj_new('idlgrpolyline', [[0,0,0],[0,1,0], [0,0,0], [0,-1,0]], color=[128,128,128], linestyle=2)
;    an2 = obj_new('idlgrpolyline', [[0,0,0],[1,0,0], [0,0,0], [-1,0,0]], color=[128,128,128], linestyle=2)
;    xobjview, [o1,o2, an0, an1, an2], background=[0,0,0], title='D'
    
    ;;hnb_ellipsoid_pointset, diff, title="D"
    
    do_adt_regress_standard_singlevoxel, [x,y,z], eval=eval, evec=evec
    hnb_ellipsoid_pointset, eval=eval, evec=evec, title='DTI'
    
    
    !P.multi = 0
    
end

pro hnb_ellipsoid_curvefit_simulation, noise_sig=noise_sig

    common scan_data, project

    ci = project.ci
        
    !P.multi = [ 0,2,2 ]
        
    sugg_filename = 'HnB_simulation_report.ps'
    
    report = dialog_pickfile(file=sugg_filename, title='Choose a destination for the report:')
    old_device = !D.name
    if (report ne '') then begin
        set_plot, 'PS'
        device, /landscape, /color, BITS_PER_PIXEL = 8, FONT_SIZE=12, language_level=2
        device, FILENAME = report, /inches, xsiz=11, ysiz=8.5, xoffset=0.25, yoffset=11-(0.25), scale_factor=0.95, /TIMES, /BOLD
        device, decomposed=1
    endif
    
    ;; sim
    cf_bval = [102.200, 1005.59, 2007.64, 3009.22, $
               4010.55, 5011.72, 6012.78, 7013.75, $
               8014.66, 9015.51, 12017.8, 15019.8, $
               18021.7, 21023.4, 24024.9]/1000
    S0     = findgen(21)*2.5 + 25.
    D      = (findgen(21)/(21-1)/2.0) * 1e-3 + 0.0003
    i = 0 & seed = 21312L
    is_sim = 1
    noise_sig = (n_elements(noise_sig) eq 0) ? 0.0 : noise_sig
    ;;;
    min_b  = min(cf_bval, max=max_b)
    b_axis = (findgen(500)/(500-1)*(max_b-min_b) + min_b);/1000.

    while (i lt n_elements(S0)) do begin

        yplot_max = S0[i] + 0.1*S0[i]
        ;; get the bvals and grad and data for the first direction
        measurements = S0[i]*exp( -((cf_bval)*D[i]*1e3) ) +$
                       abs(complex(randomn(22, n_elements(cf_bval))*noise_sig, $
                                   randomn(22, n_elements(cf_bval))*noise_sig))
        x = (y = (z = 0))
        ;; curvefit to get alpha and D
        curvefit_adc_alpha, measurements, cf_bval, A, chisquare, sigma
        yfit_str = A[0]*exp( -((cf_bval)*A[1])^A[2] )
        resid_str = yfit_str - measurements
        yerror_str = sqrt(total(resid_str^2)/3.0)
        fitline_str = A[0]*exp( -((b_axis)*A[1])^A[2] )

        curvefit_adc, measurements, cf_bval, A_o, chisquare, sigma
        A_o = [A_o, 1]
        yfit_sin = A_o[0]*exp( -((cf_bval)*A_o[1])^A_o[2] )
        resid_sin = yfit_sin - measurements
        yerror_sin = sqrt(total(resid_sin^2)/3.0) 
        fitline_sin = A_o[0]*exp( -((b_axis)*A_o[1])^A_o[2] )
        if (report ne '') then begin
            plot, cf_bval, measurements, ystyle=1, xstyle=1, yrange=[0, yplot_max], psym=1, $
                  xtitle='bvalue (*1e-3)', ytitle='Intensity', title=string(noise_sig, format='(%"Simulation, sigma=%0.2f")')
            oplot, b_axis, fitline_str
            oplot, b_axis, fitline_sin, color='bcbcbc'x
            oplot, b_axis, fitline_str+yerror_str, linestyle=1, color='ff0000'x
            oplot, b_axis, fitline_str-yerror_str, linestyle=1, color='ff0000'x
            
            xyouts, mean(cf_bval), yplot_max-(1*0.05*yplot_max), /data, string(A[2], format='(%"alpha=%0.4f")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(2*0.05*yplot_max), /data, string(A[1]*1e-3, D[i], format='(%"D=%0.5f (mm^2/s) (True D=%0.5f)")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(3*0.05*yplot_max), /data, string(A[0], S0[i], format='(%"S0=%0.3f (True S0=%0.3f)")'), charsize=0.8

            plot, cf_bval, measurements, ystyle=1, xstyle=1, yrange=[0, yplot_max], psym=1, $
                  xtitle='bvalue (*1e-3)', ytitle='Intensity', title=string(noise_sig, format='(%"Simulation, sigma=%0.2f")')
            oplot, b_axis, fitline_sin
            oplot, b_axis, fitline_str, color='bcbcbc'x
            oplot, b_axis, fitline_sin+yerror_sin, linestyle=1, color='ff0000'x
            oplot, b_axis, fitline_sin-yerror_sin, linestyle=1, color='ff0000'x
            
            xyouts, mean(cf_bval), yplot_max-(2*0.05*yplot_max), /data, string(A_o[1]*1e-3, D[i], format='(%"D=%0.5f (mm^2/s) (True D=%0.5f)")'), charsize=0.8
            xyouts, mean(cf_bval), yplot_max-(3*0.05*yplot_max), /data, string(A_o[0], S0[i], format='(%"S0=%0.3f (True S0=%0.4f")'), charsize=0.8                        
            i++
            
        endif
        
        ;; scale the gradient vector by the result and add it to an array
;        if (n_elements(alpha) eq 0) then begin
;            diff =  A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ]
;            diff = [ diff, -A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
;            
;            alpha =  A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ]
;            alpha = [ alpha, -A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
;        endif else begin
;            diff = [ diff, A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
;            diff = [ diff, -A[1]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
;            
;            alpha = [ alpha, A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
;            alpha = [ alpha, -A[2]*[ cf_grad[0,0], cf_grad[1,0], cf_grad[2,0] ] ]
;        endelse
        
;        if (n_other eq 0) then begin
;            go = 0
;            break
;        endif
        
;        gr = gr[*,other]
;        data = data[other]
                 
    endwhile

 ;   alpha = reform(alpha, 3, n_elements(alpha)/3)
 ;   diff  = reform(diff, 3, n_elements(diff)/3)
    
    ;; create the best-fit ellipse to the scaled gradients
 ;   hnb_ellipsoid_pointset, alpha, title="Alpha"
 ;   hnb_ellipsoid_pointset, diff, title="D"
    
    if (report ne '') then begin
        device, /close_file
        set_plot, old_device
    endif
    
    !P.multi = 0
    
end

pro hnb_convex_hull, matrix, title=title

    qhull, matrix, tr,  /delaunay
    poly = tetra_surface(matrix, tr)
    hull = obj_new('idlgrpolygon', matrix, color=[128,128,128], polygon=poly, alpha=0.5, reject=1)
    dots = obj_new('idlgrpolygon', matrix, style=0, thick=5, color=[128,0,0])
    card_z = obj_new('idlgrpolyline', [[0,0,0],[0,0,1], [0,0,0], [0,0,-1]]*1.2, color=[128,128,128], linestyle=2)
    card_y = obj_new('idlgrpolyline', [[0,0,0],[0,1,0], [0,0,0], [0,-1,0]]*1.2, color=[128,128,128], linestyle=2)
    card_x = obj_new('idlgrpolyline', [[0,0,0],[1,0,0], [0,0,0], [-1,0,0]]*1.2, color=[128,128,128], linestyle=2)
    axlbl = obj_new('idlgrtext', ['X','Y','Z'], locations=[[1,0,0], [0,1,0], [0,0,1]]*1.2, color=[255,255,255], /onglass)
    
    xobjview, [hull, dots, card_x, card_y, card_z, axlbl], $
              background=[0,0,0], $
              title=(n_elements(title) eq 0) ? 'HNB_CVX' : title
end

pro hnb_ellipsoid_fit_from_pointset, matrix, eval=eval, evec=evec

    n = n_elements(matrix)/3
    
    X = reform(matrix[0,*])
    Y = reform(matrix[1,*])
    Z = reform(matrix[2,*])
    
    A = [ [X^2  ], [Y^2  ], [Z^2  ], $
          [2*X*Y], [2*X*Z], [2*Y*X], $
          [2*X  ], [2*Y  ], [2*Z  ] ]
          
    A = transpose(reform(A, n, 9))
    
    ones = replicate(1.0, n)
    
    LS = invert((transpose(A) ## A)) ## transpose(A) ## ones
    
    A = [ [ LS[0], LS[3], LS[4], LS[6] ], $
          [ LS[3], LS[1], LS[5], LS[7] ], $
          [ LS[4], LS[5], LS[2], LS[8] ], $
          [ LS[6], LS[7], LS[8], -1.0  ] ]
        
    center  = transpose(-invert(A[0:2, 0:2]) ## [LS[6:8]])
    
    print, "--------------------"
    print, "center:", center
    
    eval = eigenql(A[0:2,0:2], eigenvectors=evec)
    eval = 1.0/sqrt(eval)

end

pro hnb_ellipsoid_pointset, matrix, eval=eval, evec=evec, title=title

    common scan_data, project

    if (n_elements(matrix) ne 0) then begin
        hnb_ellipsoid_fit_from_pointset, matrix, eval=eval, evec=evec
    endif else begin
        eval /= max(eval)
    endelse
    
    eval_sort = reverse(sort(eval))
    eval = eval[eval_sort]
    evec = evec[*,eval_sort]
    
    if (total(finite(eval, /nan)) ne 0) then begin
        hnb_convex_hull, matrix, title=title
        return
    endif
    
    print, "--------------------"
    print, eval
    print, "--------------------"
    print, evec
    print, "--------------------"

    evec = transpose(evec)
    
    surf = sp_make_sphere(30, vert=vert, poly=poly)
    sz = (size(vert))[2]
    shift = transpose([ [replicate(-(30/2.), sz)], $
                        [replicate(-(30/2.), sz)], $
                        [replicate(-(30/2.), sz)] ])
    vert += shift
    for v = 0, n_elements(vert)/3 - 1 do begin
        n = norm(vert[*,v])
        if (n ne 0) then begin
            vert[*,v] /= norm(vert[*,v])
        endif
    endfor

    ;; shift, rotate, de-shift    
    evec = [ [evec[*,0], 0], $
             [evec[*,1], 0], $
             [evec[*,2], 0], $
             [0, 0, 0, 1] ]  
             
    tx = evec ## diag_matrix([eval, 1.0])
    vert = vert_t3d(vert, /no_copy, matrix=tx)
    
    o1 = obj_new('idlgrpolygon', vert, color=[128,128,128], polygon=poly, alpha=0.5, reject=1)
    if (n_elements(matrix) ne 0) then begin
        o2 = obj_new('idlgrpolygon', matrix, style=0, thick=5, color=[128,0,0])
    endif
    ln0 = obj_new('idlgrpolyline', [[0,0,0],[reform(evec[0,0:2])], [0,0,0], [-reform(evec[0,0:2])]]*1.2, color=[0,128,0])
    ln1 = obj_new('idlgrpolyline', [[0,0,0],[reform(evec[1,0:2])], [0,0,0], [-reform(evec[1,0:2])]]*1.2, color=[0,128,0])
    ln2 = obj_new('idlgrpolyline', [[0,0,0],[reform(evec[2,0:2])], [0,0,0], [-reform(evec[2,0:2])]]*1.2, color=[0,128,0])
    axlbl1 = obj_new('idlgrtext', ['E1','E2','E3'], locations=[[reform(evec[0,0:2])], [reform(evec[1,0:2])], [reform(evec[2,0:2])]]*1.2, color=[0,128,0], /onglass)
    an0 = obj_new('idlgrpolyline', [[0,0,0],[0,0,1], [0,0,0], [0,0,-1]]*1.2, color=[128,128,128], linestyle=2)
    an1 = obj_new('idlgrpolyline', [[0,0,0],[0,1,0], [0,0,0], [0,-1,0]]*1.2, color=[128,128,128], linestyle=2)
    an2 = obj_new('idlgrpolyline', [[0,0,0],[1,0,0], [0,0,0], [-1,0,0]]*1.2, color=[128,128,128], linestyle=2)
    axlbl2 = obj_new('idlgrtext', ['X','Y','Z'], locations=[[1,0,0], [0,1,0], [0,0,1]]*1.2, color=[255,255,255], /onglass)
    
    if (n_elements(matrix) ne 0) then begin
        xobjview, [o1,o2, ln0, ln1, ln2, an0, an1, an2, axlbl1, axlbl2], background=[0,0,0], title=(n_elements(title) ne 0) ? title : "Ellipse"
    endif else begin
        xobjview, [o1, ln0, ln1, ln2, an0, an1, an2, axlbl1, axlbl2], background=[0,0,0], title=(n_elements(title) ne 0) ? title : "Ellipse"
    endelse
    
end