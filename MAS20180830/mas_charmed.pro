pro UCB_double_up, coord, state, simdata=simdata

    common scan_data
    
    bmat = *project.imndarray[0].b_matrix
    
    do_adt_regress_standard_singlevoxel, [37,49,0], $
        fa=fa, eval=eval, evec=evec, ad=ad, fit=fit
        
    evals = eval[sort(eval)] * 1;e1; ; [ 1.0, 0.3, 0.3 ]*1e-4
    evals = [ 1.0, 0.1, 0.1 ] * 1e-3
    
    vec1 = 1.0/sqrt(2) * [1,1,0]
    vec1_a = vec1 ## [ [ cos(!PI/2), -sin(!PI/2), 0 ], [ sin(!PI/2), cos(!PI/2), 0], [0,0,1] ]
    vec1_b = crossp(vec1_a, vec1)
    evecs1 = [ [vec1], [vec1_a], [vec1_b] ]
    tens1 = transpose(evecs1) ## diag_matrix(evals) ## evecs1
    tens1 = [tens1[0,0], tens1[1,1], tens1[2,2], tens1[0,1], tens1[0,2], tens1[1,2]]
    
    vec2 = 1.0/sqrt(2) * [1,-1,0]
    vec2_a = vec2 ## [ [ cos(!PI/2), -sin(!PI/2), 0 ], [ sin(!PI/2), cos(!PI/2), 0], [0,0,1] ]
    vec2_b = crossp(vec2_a, vec2)    
    evecs2 = [ [vec2], [vec2_a], [vec2_b] ]
    tens2 = transpose(evecs2) ## diag_matrix(evals) ## evecs2
    tens2 = [tens2[0,0], tens2[1,1], tens2[2,2], tens2[0,1], tens2[0,2], tens2[1,2]]
    
    frac = 0.25
    
    simdata = 30.0 * ( frac*exp(- total(bmat*rebin(tens1, 6, n_elements(bmat)/6), 1) ) + $
                       (1-frac)*exp(- total(bmat*rebin(tens2, 6, n_elements(bmat)/6), 1) ) )
    
    (*project.dataarray[0].state1)[30,30,0,*] = simdata
    message, /info, "Done."

end

function charmed_oh_my_god, dfit, grad_dirs, maxG, delta, smalldel, te, nhin, nres, eta

    common charmed_fit, data
    
    gamma = 4257.0
;    
    maxq = gamma * smalldel * maxG/10e6
;    
;    grad_dirs = grad_dirs * maxq
;    
    grad_sph = cv_coord(from_rect=grad_dirs, /to_sphere)
    
    theta_q = reform(grad_sph[0,*])
    phi_q   = !PI/2 - reform(grad_sph[1,*])
    R_q     = reform(grad_sph[2,*])
 
    nerve= [ [0,1,0], [1,1,0] ]
    nnerves = 2L
    
    f_h = 0.7
    f_r = replicate(1-f_h, nnerves)/nnerves ;,0.12]
    
    temp = cv_coord(from_rect=Nerve, /to_sphere)
    theta_N = transpose([temp[0, 0:nnerves-1]]) ;;[ 0.785398, -0.785398 ]
    phi_N = !PI/2 - transpose([temp[1, 0:nnerves-1]]) ;;[ 1.57080, 1.57080 ]
    D_r = replicate(1.9, nnerves) ;, 1.9 ]
    
    D_vec = float([ 1, 1, 1, 0, 0, 0 ]) * 1e-3

    axon = findgen(6)+1.5
    ;axon /= 2
    a = axon
    R = a
    
    weight = [1.826, 9.2374, 16.7562, 22.986, 18.525, 16.8636]
    weight /= total(weight)
    
    tau_coeff = te/2.0
    
    data = { grad_dirs: grad_dirs, $
             theta_q: theta_q, $
             phi_q: phi_q, $
             R_q: R_q, $
             ydata: replicate(0, n_elements(grad_Dirs)/3), $
             nhin: nhin, $
             nres: nres, $
             delta: delta, $
             smalldel: smalldel, $
             weight: weight, $
             tau_coeff: tau_coeff, $
             maxQ: maxQ, $
             R: R }
    
    X0 = [ f_h, f_r, theta_N, phi_N, D_vec, D_r, eta ]
    E = charmed_model(grad_Dirs, X0)
    return, E ;;sqrt(E^2 + eta^2); + randomn(1L, n_elements(grad_dirs)/3)*eta
    
    
end


function charmed_regrid, dim

    arr   = float(array_indices([dim,dim,dim], lindgen(dim^3), /dimensions)) - dim/2.0
    norms = sqrt(total(arr^2, 1))
    norms[*] = max(norms)
    arr  /= rebin(transpose(norms), 3, n_elements(norms))
    
    return, arr

;    grid = fltarr(3,dim^3)
;    n = 0L
;    norms = fltarr(dim^3)
;    for x = 0, dim-1 do begin
;        for y = 0, dim-1 do begin
;            for z = 0, dim-1 do begin
;                grid[*,n] = ([x,y,z] - (dim-1)/2.0);/(dim/2.0)
;                norms[n] = sqrt(total(grid[*,n]^2))
;                ++n
;            endfor
;        endfor
;    endfor
;
;    return, grid / rebin(replicate(max(norms), 3), 3, dim^3) ;(grid - dim/2.0)/((dim-1)/2)
    
end

function charmed_model, grad, x, E_h=E_h, E_r=E_r

    common charmed_fit, data
    
    theta_q   = data.theta_q
    phi_q     = data.phi_q
    R_q       = data.R_q
    ydata     = data.ydata
    nhin      = data.nhin
    nres      = data.nres
    delta     = data.delta
    smalldel  = data.smalldel
    weight    = data.weight
    tau_coeff = data.tau_coeff
    R         = data.R
    if (0 and n_elements(grad) ne 0) then begin
        grad_dirs = grad * data.maxq
        grad_sph = cv_coord(from_rect=grad_dirs, /to_sphere)
        theta_q  = reform(grad_sph[0,*])
        phi_q    = !PI/2 - reform(grad_sph[1,*])
        R_q      = reform(grad_sph[2,*])
    endif else begin
        grad_dirs = data.grad_dirs
    endelse
    
    if (nhin gt 0) then begin
        f_h      = x[0:nhin-1]
        D_vec    = x[nhin+3*nres:nhin+3*nres+6-1];
        eta      = x[nhin+4*nres+6];
    endif
    if (nres gt 0) then begin
        f_r      = x[nhin:nhin+nres-1]
        theta_N  = x[nhin+nres:nhin+2*nres-1]
        phi_N    = x[nhin+2*nres:nhin+3*nres-1];
        D_r      = x[nhin+3*nres+6:nhin+4*nres+6-1];
    endif
    
    Dxx      = D_vec[0];
    Dyy      = D_vec[1];
    Dzz      = D_vec[2];
    Dxy      = D_vec[3];
    Dxz      = D_vec[4];
    Dyz      = D_vec[5];
    D_mat    = [[Dxx, Dxy, Dxz], $
                [Dxy, Dyy, Dyz], $
                [Dxz, Dyz, Dzz]]
    a = R
    
    l_q = n_elements(R_q)
    l_a = n_elements(a)
    
    R_matrix = rebin(R, l_a, l_q, /sample)
    gamma = weight
    gamma_matrix = rebin(gamma, l_a, l_q, /sample)
    
    E_h = fltarr(l_q, nhin)
    E_r = fltarr(l_q, nres)
    E_tot = fltarr(l_q)
    
    diag = lindgen(n_elements(E_tot))
    for i = 0, nhin-1 do begin
        ;gT_D_g = (grad_dirs ## D_mat ## transpose(grad_dirs))[diag,diag]
        gT_D_g = (total( (grad_dirs ## D_mat) * grad_dirs, 1))
        E_h[*,i] = f_h[i]*exp(-4*!PI^2*(delta-smalldel/3.0)*gT_D_g)
        ;E_tot = E_tot + f_h[i]*E_h
    endfor
    E_tot = E_tot + total(reform(E_h, l_q, nhin), 2)
    
    for j = 0, nres-1 do begin
        factor_angle_term_par = abs(sin(theta_q)*sin(theta_N[j])*cos(phi_q-phi_N[j])+cos(theta_q)*cos(theta_N[j]))
        factor_angle_term_perp = sqrt(1-factor_angle_term_par^2)
        
        q_par_sq = (R_q * factor_angle_term_par)^2
        q_par_sq_matrix = rebin(q_par_sq, l_q, l_a, /sample)
        q_perp_sq = (R_q * factor_angle_term_par)^2
        
        E = exp(-4 * !PI^2 * transpose(q_par_sq_matrix) * (delta-smalldel/3.0) * D_r[j]) * $
            exp(-4 * !PI^2 * q_perp_sq##R^4/(D_r[j]*tau_coeff)*(7.0/96.0)*(2-(99.0/112.0)*transpose(R_matrix)^2/(D_r[j]*tau_coeff)))
        E_r[*,j] = total(E * gamma_matrix,1) * f_r[j]
        ;E_tot = E_tot + f_r[j]*E_r
    endfor
    E_tot = E_tot + total(reform(E_r, l_q, nres), 2)
    
    E_tot = sqrt(E_tot^2 + eta^2)

    return, E_tot
    
end

function charmed, ydata, dfit, grad_dirs, maxG, delta, smalldel, te, nhin, nres, $
                  yfit=yfit, bestnorm=bestnorm

    common charmed_fit, data
    
    gamma = 4257.0
    
    maxq = gamma * smalldel * maxG/10e6
    
    grad_dirs = grad_dirs * maxq
    
    grad_sph = cv_coord(from_rect=grad_dirs, /to_sphere)
    
    theta_q = reform(grad_sph[0,*])
    phi_q   = !PI/2 - reform(grad_sph[1,*])
    R_q     = reform(grad_sph[2,*])
    
    if (nhin gt 0) then begin
        f_h = 0.7 * replicate(1, nhin)/nhin
    endif else begin
        f_h = 0.0
    endelse
    
    if (nres gt 0) then begin
        f_r = 0.3 * replicate(1, nres)/nres
    endif else begin
        f_r = 0.0
    endelse
    
    ;dfit =[ 0.0348 ,   0.1236 ,   0.0336 ,  -0.0041   ,-0.0010 ,   0.0000]*1e-3
    Dxx = dfit[0] * 1e3
    Dyy = dfit[1] * 1e3
    Dzz = dfit[2] * 1e3
    Dxy = dfit[3] * 1e3
    Dxz = dfit[4] * 1e3
    Dyz = dfit[5] * 1e3
    
    D_vec = dfit * 1e3
    
    D_mat = [ [Dxx, Dxy, Dxz], [Dxy, Dyy, Dyz], [Dxz, Dyz, Dzz] ]
    
    eta = 0.08
    
    D = eigenql(D_mat, eigenvectors=V)
    void = max(D, max_eigval_ind)
    
    if (nres gt 0) then begin
        D_r = replicate(1.9, nres) ;; orig: 1.9
        N_sph = cv_coord(from_rect=reform(V[*, max_eigval_ind]), /to_sphere)
        theta_N = replicate(reform(N_sph[0]), nres)
        phi_N   = replicate(!PI/2.0 - reform(N_sph[1]), nres)
        X0 = [ f_h, f_r, theta_N, phi_N, D_vec, D_r, eta ]
        min_val = [replicate(0,nhin), replicate(0,nres), $
            replicate(-4*!PI, nres), $
            replicate(-4*!PI, nres), $
            0, 0, 0, -1, -1, -1, $
            replicate(0.9,nres),  0.0001 ]
            
        max_val = [replicate(1,nhin), replicate(1,nres), $
            replicate(4*!PI, nres), $
            replicate(4*!PI, nres), $
            2, 2, 2, 1, 1, 1, $
            replicate(2.0,nres),  0.10 ]
    endif else begin
        X0 = [ f_h, D_vec, eta ]
        min_val = [replicate(0,nhin), $
                   0, 0, 0, -1, -1, -1, 0.0001 ]
            
        max_val = [replicate(1,nhin), $
                   2, 2, 2, 1, 1, 1, 0.10 ]
    endelse
    
    little = where(abs(X0) lt 1e-4, nlittle)
    if (nlittle gt 0) then X0[little] = 0.0
    
    parinfo = replicate({ value:0.D, fixed:0, limited:[0,0], limits:[0.D,0] }, n_elements(X0))  
    
    for i = 0, n_elements(X0)-1 do begin
        parinfo[i].limited = [1,1]
        parinfo[i].limits = [min_val[i], max_val[i]]
        ;parinfo[i].value = X0[i]
    endfor
    ;parinfo[n_elements(parinfo)-1].fixed = 1B
    
    axon = findgen(6)+1.5
    ;axon /= 2
    a = axon
    R = a
    
    weight = [1.826, 9.2374, 16.7562, 22.986, 18.525, 16.8636]
    weight /= total(weight)
    
    tau_coeff = te/2.0
    
    data = { grad_dirs: grad_dirs, $
             theta_q: theta_q, $
             phi_q: phi_q, $
             R_q: R_q, $
             ydata: ydata, $
             nhin: nhin, $
             nres: nres, $
             delta: delta, $
             smalldel: smalldel, $
             weight: weight, $
             tau_coeff: tau_coeff, $
             maxQ: maxQ, $
             R: R }
    
    fitted = MPFITFUN('charmed_model', grad_dirs, ydata, replicate(1.0, n_elements(ydata)), x0, $
                      yfit=yfit, quiet=0, parinfo=parinfo, bestnorm=bestnorm) 
                      
    return, fitted

end

pro UCB_charmed_fit, coord, state

    common scan_data
    common charmed_fit, data
    
    if (ptr_valid(state)) then begin
        ci = (*state).proj_index
    endif else begin
        ci = project.ci
    endelse
    
    ;; Obtain the gradient amplitude values
    method_str = project.imndarray[ci].imnd_file
    tmp_amps = stregex(method_str, 'PVM_DwGradAmp=\( ([0-9]+) \)([0-9 \.]+)', /extract, /subex)
    if (n_elements(tmp_amps) eq 3) then begin
        amps = float(strsplit(tmp_amps[2], ' ', /extract))
    endif
    namps = n_elements(amps)
    
    ;; Obtain the b-values to indicate how many scans there are
    bvals = *project.imndarray[ci].bval_array
    amps_matrix = transpose((rebin(amps, namps, n_elements(bvals)/namps, /sample))[*])
    amps_matrix = rebin(amps_matrix, 3, n_elements(amps_matrix), /sample)
    maxG = max(amps)
    
    ;; obtain echo time
    te = project.imndarray[ci].echo_time
    
    ;; convert gradient directions from sphereical to rectangular
    grad_dirs = cv_coord(from_sphere=transpose([ [transpose([*project.imndarray[0].angle_phi])], $
                                    [transpose([90-*project.imndarray[0].angle_theta])], $
                                    [replicate(1, n_elements(bvals))] ]), $
                         /to_rect, /degrees)
    grad_dirs *= (amps_matrix/maxG)

    ;; obtain delta
    smalldel = (*project.imndarray[0].small_delta)[0]
    delta    = (*project.imndarray[0].big_delta)[0]
    
    ;; Specify the number of hindered and number of restricted compartments
    nhin = 1
    nres = 1
    
    xdim = 142L
    ydim = 94L

    x = coord[0]
    y = coord[1]

    ydata = reform((*project.dataarray[ci].state1)[x, y, 0, *])
    ;; I am replacing the ydata with simulated y data
    ;ydata = charmed_oh_my_god(dfit, grad_dirs, maxG, delta, smalldel, te, nhin, nres, 0.09) * 30.0
    ;(*project.dataarray[ci].state1)[0,0,0,*] = ydata

    ;; do DTI fitting for the current voxel
    do_adt_regress_standard_singlevoxel, [x,y,0], eval=eval, evec=evec, fit=dfit, bval_threshold=2500

    ;; I am replacing it with generic information.
    ;;dfit = [1.0,1.0,1.0,0,0,0] * 1e-3
    dfit = dfit[1:*]
    ;; obtain the ydata from data array 
    
    ;; Obtain an estimate for S0 using the many gradient directions.
    S0 = fltarr(n_elements(ydata)/namps)
    for i = 0, n_elements(ydata)/namps-1 do begin 
        j = ydata[i*namps:i*namps+3]
        curvefit_adc, j, bvals[0:3]*1e-3, A
        S0[i] = a[0]
    endfor

    ;; compute measured E(q)
    E = ydata/mean(S0)
    
    ;; perform the fit to charmed model
    fit = charmed(E, dfit, grad_dirs, maxG, delta, smalldel, te, nhin, nres, $
                  bestnorm=bestnorm, yfit=yfit)

    ;; unpack the fitted parameters
    if (nres gt 0) then begin
        f_r      = fit[nhin:nhin+nres-1]
        theta_N  = fit[nhin+nres:nhin+2*nres-1]
        phi_N    = fit[nhin+2*nres:nhin+3*nres-1]
        D_r      = fit[nhin+3*nres+6:nhin+4*nres+6-1]

        n_rect = cv_coord(from_sphere=[ transpose(theta_N), $
                                         transpose(!PI/2-phi_n), $
                                         transpose(replicate(1,nres)) ], $
                          /to_rect)
        
        print, "f_r.......: ", strjoin(strtrim(f_r,2), ', ')
        print, "theta_N...: ", strjoin(strtrim(theta_N*1.0/!DTOR,2), ', ')
        print, "phi_N.....: ", strjoin(strtrim(phi_N*1.0/!DTOR,2), ', ')
        print, "D_r.......: ", strtrim(D_r, 2)
        print, "N_rect....: ", '[ '+strjoin(strtrim(n_rect,2), ', ')+' ]'
        
    endif
    f_h      = fit[0:nhin-1]
    D_vec    = fit[nhin+3*nres:nhin+3*nres+6-1]
    eta      = fit[nhin+4*nres+6]
    print, "D_vec.....: ", strjoin(strtrim(D_vec,2), ', ')
    print, "Evec_in...: ", '[ '+strjoin(strtrim(reform(evec[*,0]),2), ',')+' ]'
    
    D_ten = [ [ D_vec[0], D_vec[3], D_vec[4] ], $
              [ D_vec[3], D_vec[1], D_vec[5] ], $
              [ D_vec[4], D_vec[5], D_vec[2] ] ]*1e-3
    eval_out = eigenql(D_ten, eigenvectors=evec_out)
    print, "Evec_ot...: ", '[ '+strjoin(strtrim(reform(evec[*,0]),2), ',')+' ]'
    print, "f_h.......: ", strjoin(strtrim(f_h,2), ', ')
    print, "eta.......: ", strtrim(eta,2)
    print, "bestnrm...: ", strtrim(bestnorm, 2)
    print, "ssqresd...: ", strtrim(sqrt(total( (E-yfit)^2 )), 2)
    print, "angledif..: ", strtrim(1.0/!DTOR*acos(abs(total(n_rect*reform(evec[*,0])))),2)
    ;; do a plot
    do_plot = 'No'
    do_plot = dialog_message("Show Plot?", /question, /center)
    if (do_plot eq 'Yes') then begin
        !P.multi = [0, 4, 6]
        window, 5, xsize=1700, ysize=1000
        device, decomposed=1
        for i = 0, 21-1 do begin
            plot, data.R_q, yfit[i*15:i*15+15-1], yrange=[0,1], xstyle=1, ystyle=1,$
                  xtitle='q', ytitle='E_tot(q)'
            oplot, data.R_q, E[i*15:i*15+15-1], psym=1, color='0000FF'x
            
        endfor
        !P.multi = 0
    endif
    
    ;;return
    
    if (1) then begin
        dim = long(275)
        grid = charmed_regrid(dim)
        gamma = 4257.0
        maxq = gamma * smalldel * maxG/10e6
        grad_dirs = grid * 5.0 * maxq
        grad_sph = cv_coord(from_rect=grad_dirs, /to_sphere)
        theta_q = reform(grad_sph[0,*])
        phi_q   = !PI/2 - reform(grad_sph[1,*])
        R_q     = reform(grad_sph[2,*])
        data = { grad_dirs: grad_dirs, $
                 theta_q: theta_q, $
                 phi_q: phi_q, $
                 R_q: R_q, $
                 ydata: ydata, $
                 nhin: data.nhin, $
                 nres: data.nres, $
                 delta: data.delta, $
                 smalldel: data.smalldel, $
                 weight: data.weight, $
                 tau_coeff: data.tau_coeff, $
                 maxQ: data.maxQ, $
                 R: data.R }

        E_regrid = charmed_model(grad_dirs, fit, E_h=E_h, E_r=E_r)

        int_factor = 1
        
        volh = reform(E_h[*,0], dim, dim, dim)
        volh = shift(abs(fft(volh, 1, /overwrite)), int_factor*dim/2, int_factor*dim/2, int_factor*dim/2)
        shade_volume, volh, max(volh)/4, v, p
        if (n_elements(v) gt 0) then begin
            o_h = obj_new('idlgrpolygon', v, polygons=p, color=[128,255,128], shading=1)
            xobjview, [o_h]
        endif

        volr = reform(E_r[*,0], dim, dim, dim)
        volr = shift(abs(fft(volr, 1, /overwrite)),int_factor*dim/2, int_factor*dim/2, int_factor*dim/2)
        shade_volume, volr, max(volr)/4, v, p
        if (n_elements(v) ne 0) then begin
            o_r = obj_new('idlgrpolygon', v, polygons=p, color=[128,128,255], shading=1)
            xobjview, [o_r]
        endif
        
        if (obj_valid(o_r) and obj_valid(o_h)) then begin
            xobjview, [o_h, o_r]
        endif
        
    endif
    
    heap_gc
    
    print, "Done."
    
end



