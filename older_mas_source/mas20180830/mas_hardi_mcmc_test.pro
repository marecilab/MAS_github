;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; The following code is proof-of-concept / work in progress / experimental
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; EXPRIMENTAL EXPRIMENTAL EXPRIMENTAL EXPRIMENTAL EXPRIMENTAL
; EXPRIMENTAL

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mas_HARDI_MCMC_simulator::init, $
   sigma=sigma, $
   data=data, $
   weights=weights, $
   seed1=seed1, $
   seed2=seed2, $
   burn_time=burn_time, $
   jump_interval=jump_interval, $
   start_direction=start_direction, $
   reco_object=reco_object
  
    if (keyword_set(sigma)) then begin
       self.sigma = float(sigma)
    endif

    if (keyword_set(seed1)) then begin
       temp = long(seed1)
    endif else begin
       temp = 100L
    endelse
    void = randomu(temp)
    self.seed1 = temp

    if (keyword_set(seed2)) then begin
       self.seed2 = long(seed2)
    endif

    if (keyword_set(burn_time)) then begin
       self.burn_time = long(burn_time)
    endif

    if (keyword_set(jump_interval)) then begin
       self.jump_interval = long(jump_interval)
    endif

    if (keyword_set(start_direction)) then begin
       self.start_direction = start_direction
    endif else begin
       self.start_direction = float([0,0,0])
    endelse

    if (keyword_set(reco_object)) then begin
       if (obj_valid(reco_object)) then begin
          self.reco_obj = reco_object
       endif
    endif else begin
       message, 'This object requires an existing diffusion reconstructor.'
       return, 0
    endelse

    if keyword_set(weights) then begin
       self.curr_weights = ptr_new(weights)
    endif else if keyword_set(data) then begin
       reco_object->getProperty, can_use_weights=can_use_weights
       if (can_use_weights) then begin
          w = reco_object->computeWeights(data)
          self.curr_weights = ptr_new(w)
       endif else begin
          self.curr_data = ptr_new(data)
       endelse
    endif else begin
       message, 'Data or Weights required.'
       return, 0
    endelse

    mas_tessellator_make_icos, level=3, vertlist=vlist, polylist=plist
    pr = (self.reco_obj)->getDisplacementProbability(weights=*self.curr_weights, $
                                                     direction=vlist)
    max_pr = max(pr, min=min_pr)
    self.max_pr = max_pr
    self.min_pr = min_pr
    
    return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_hardi_MCMC_simulator::reset


end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mas_HARDI_MCMC_simulator::jump, n_jumps, accept_ratio=accept_ratio

    forward_function mas_mcmc_studt_n
    forward_function mas_mcmc_beta

    old_x = (self.last_direction)[0]
    old_y = (self.last_direction)[1]
    old_z = (self.last_direction)[2]
    seed = self.seed1

;; Student's t Dist
;     n = 5    
;     new_x = mas_mcmc_studt_n(n, seed=seed) * self.sigma + old_x
;     new_y = mas_mcmc_studt_n(n, seed=seed) * self.sigma + old_y
;     new_z = mas_mcmc_studt_n(n, seed=seed) * self.sigma + old_z
;;

;; Beta Dist
     a = 8 & b = 12
     new_x = mas_mcmc_beta(a,b, seed=seed, prob=p_new_x) + old_x
     new_y = mas_mcmc_beta(a,b, seed=seed, prob=p_new_y) + old_y
     new_z = mas_mcmc_beta(a,b, seed=seed, prob=p_new_z) + old_z
;;

;; Normal Dist    
;     new_x = randomu(seed, /normal) * self.sigma + old_x
;     new_y = randomu(seed, /normal) * self.sigma + old_y
;     new_z = randomu(seed, /normal) * self.sigma + old_z
;;

    direction = [ new_x, new_y, new_z ]/sqrt(new_x^2 + new_y^2 + new_z^2)
  
    pr_prop = (self.reco_obj)->getDisplacementProbability(weights=*self.curr_weights, $
                                                          direction=direction)

    pr_prop = (pr_prop[0] - self.min_pr)/(self.max_pr - self.min_pr)
    
    ;; proposal density assumed to be symmetric. note: beta may not be symm.!
    arg = min([1., pr_prop/(self.pr_last > 1e-8)]) 

    self.nth_sample++

    if (arg eq 1. or arg lt randomu(self.seed1)) then begin

       self.n_accept++
       self.accept_ratio = float(self.n_accept)/float(self.nth_sample)

       if (self.nth_sample ge self.burn_time) then begin
          
       endif 
       
       self.last_direction = direction
       self.pr_last = pr_prop
       
    endif
    
    self.seed1 = seed
    accept_ratio = self.accept_ratio

    return, direction

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_HARDI_MCMC_simulator::getProperty, $
   sigma=sigma, $
   accept_ratio=accept_ratio, $
   naccept=naccept

     sigma = self.sigma
     accept_ratio = self.accept_ratio
     naccept = self.n_accept

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_HARDI_MCMC_simulator::setProperty, sigma=sigma
   
    if keyword_set(sigma) then self.sigma = float(sigma)

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_HARDI_MCMC_simulator::cleanup

    ptr_free, self.curr_weights

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_HARDI_MCMC_simulator__define

    struct = { MAS_HARDI_MCMC_SIMULATOR, $
               sigma: 0.15, $
               curr_weights: ptr_new(), $
               curr_data: ptr_new(), $
               start_direction: fltarr(3), $
               last_direction: fltarr(3), $
               pr_last: double(1e-5), $
               max_pr: 0D, $
               min_pr: 0D, $
               seed1: lonarr(36), $
               seed2: 100L, $
               curr_t: 0L, $
               accept_ratio: 0.0, $
               n_accept: 0L, $
               nth_sample: 0L, $
               burn_time: 1000L, $
               jump_interval: 10L, $
               reco_obj: obj_new() }

end

;; randomly sample from student's t
function mas_mcmc_studt_n, n, seed=seed

    if (not keyword_set(seed)) then seed = 100L; systime(1)

    R1 = randomu(seed)

    if (n gt 1) then begin

       R2 = randomu(seed)
       cos_2pir2 = cos(2*!PI*R2)
       
       T_n = (sqrt(n) * cos_2pir2) / (sqrt(1.0/(1.-R1^(2./(n-1))) - cos_2pir2^2))
       
    endif else begin

       T_n = tan ((R1 - 0.5) * !PI)

    endelse

    return, T_n
end

;; randomly sample beta(a,b) dist
function mas_mcmc_beta, a, b, seed=seed, prob=prob, compute_prob=compute_prob

    center = float(a)/(a+b)

    G1 = randomu(seed, gamma=a)
    G2 = randomu(seed, gamma=b)
       
    B_ab = G1/(G1 + G2)

    if (keyword_set(compute_prob)) then begin
       prob = ( B_ab^(a) * (1.0-B_ab)^(b) )/BETA(a,b)
    endif

    ;print, B_ab, prob

    ;; center it at the mean for our purposes
    return, B_ab - center

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mas_hardi_mcmc_test, directions=directions, $
                         sample_every=sample_every, $
                         comp_method=comp_method, $
                         decon_level=decon_level, $
                         coord=coord, $
                         A_arg=A_arg, $
                         weights=weights, $
                         pr_vlist=pr_vlist, $
                         xobj_tlb=xobj_tlb
    
    common mas_hardi_mcmc_test, $
       seed1, seed2,   $
       ntimes,         $
       sigma1, sigma2, $
       theta, phi,     $
       data, b_mat, time 

    common probtrack_objective, Pr, vlist, pr_interp_tol

    common scan_data, project

    b_mat = *project.imndarray[0].b_matrix * 1e-3
    time  = ((*project.imndarray[0].big_delta)-(*project.imndarray[0].small_delta/3D))[0]
    if (not keyword_set(coord)) then begin
       coord = [ project.procpramarray[0].fdim_start, $
                 project.procpramarray[0].pdim_start, $
                 project.procpramarray[0].sdim_start ]
    endif

    comp_method  = keyword_set(comp_method) ? comp_method : 0
    decon_level  = keyword_set(decon_level) ? decon_level : 2
    data         = reform((*project.dataarray[0].state1)[coord[0], coord[1], coord[2], *])
    sample_every = keyword_set(sample_every) ? long(sample_every) : 1L
    nsamples     = keyword_set(ntimes) ? ntimes : 50000L
    seed1        = keyword_set(seed1) ? long(seed1) : 100L
    sigma1       = keyword_set(sigma1) ? sigma1 : 0.09
    burn_time    = 1L; 1000L

    mm = obj_new('mas_mow', $
                 b_matrix=b_mat, $
                 diff_time=time, $
                 radius=15.0, $
                 deconvolution_level=decon_level,$
                 computation_method=comp_method, $
                 debug=1)

    weights = mm->computeWeights(data)

    mas_tessellator_make_icos, level=3, vertlist=vlist, polylist=plist
    Pr = mm->getDisplacementProbability(weights=weights, direction=vlist)
    max_pr = max(Pr, min=min_pr)
    Pr_vlist = (Pr - min_pr)/(max_pr - min_pr)
    
    mcmc = obj_new('mas_HARDI_MCMC_simulator', $
                   weights=weights, $
                   reco_object=mm, $
                   seed1=seed1, $
                   sigma=sigma1, $
                   burn_time=burn_time)

    directions = fltarr(3, (nsamples)/sample_every)

    time_start = systime(1)
    nth        = 0L
    ac_ratios  = fltarr(nsamples+burn_time)
    sigmas     = fltarr(nsamples+burn_time)
    sigma_tune = sigma1

    print, "mas_hardi_mcmc_test: using coord.....: ", coord

    for i = 0L, (nsamples+burn_time)-1 do begin

       if (i mod sample_every eq 0 and i gt burn_time) then begin
          directions[*,nth] = mcmc->jump(1, accept_ratio=ar)
          nth++
       endif else begin
          junk = mcmc->jump(1, accept_ratio=ar)
       endelse

       ac_ratios[i] = ar
       sigmas[i]    = sigma_tune

       if (0 and i mod 100 eq 0) then begin
          if (ar lt 0.21) then begin
;             sigma_tune = (sigma_tune + 0.01) < 0.25
             sigma_tune = (sigma_tune * (1. + 0.008)) < 0.2
             mcmc->setProperty, sigma = sigma_tune
          endif else if (ar gt 0.22) then begin
;             sigma_tune = (sigma_tune - 0.01) > 0.001
             sigma_tune = (sigma_tune * (1. - 0.008)) > 0.001
             mcmc->setProperty, sigma = sigma_tune
          endif
       endif

    endfor

    mcmc->getProperty, naccept=acceptance, $
                       accept_ratio=accept_ratio, $
                       sigma=sigma_tune

    print, "mas_hardi_mcmc_test: Elapsed Time....: ", systime(1) - time_start
    print, "mas_hardi_mcmc_test: Total (w/burn)..: ", i
    print, "mas_hardi_mcmc_test: # Samples.......: ", nth
    print, "mas_hardi_mcmc_test: Sigma (start)...: ", sigma1
    print, "mas_hardi_mcmc_test: Sigma (end).....: ", sigma_tune
    print, "mas_hardi_mcmc_test: Acceptance......: ", acceptance
    print, "mas_hardi_mcmc_test: Acceptance Rate.: ", accept_ratio

    omodel = obj_new('idlgrmodel', depth_test_disable=2)

    o = obj_new('idlgrpolygon', directions, $
                style=0, color=[255,0,0], $
                alpha_channel=0.50)

    vlist_p = fltarr(size(vlist, /dimensions))
    vlist_p[0,*] = vlist[0,*]*pr_vlist
    vlist_p[1,*] = vlist[1,*]*pr_vlist
    vlist_p[2,*] = vlist[2,*]*pr_vlist

    sp = obj_new('idlgrpolygon', vlist_p, polygons=plist, $
                 shininess=30., $
                 shading=0, texture_interp=0, $
                 style=1, color=[200,200,120], $
                 depth_test_disable=0, $
                 alpha_channel=0.15)

    omodel->add, [o, sp];;model

    xobjview, omodel, background=[0,0,0], renderer=0, tlb=xobj_tlb
    
    x_samp = reform(directions[0,*])
    x_hist = histogram(x_samp, locations=x_loc, nbins=100)

    y_samp = reform(directions[1,*])
    y_hist = histogram(y_samp, locations=y_loc, nbins=100)

    z_samp = reform(directions[2,*])
    z_hist = histogram(z_samp, locations=z_loc, nbins=100)

    !P.multi = [0, 2, 3]
    window, /free, xsize=400*2, ysize=600
    plot, x_loc, x_hist/total(x_hist), TITLE='X probability'
    plot, x_samp, TITLE='X Time Series', xrange=[0,n_elements(x_samp)], xstyle=1
    plot, y_loc, y_hist/total(y_hist), TITLE='Y probability'
    plot, y_samp, TITLE='Y Time Series', xrange=[0,n_elements(y_samp)], xstyle=1
    plot, z_loc, z_hist/total(z_hist), TITLE='Z probability'
    plot, z_samp, TITLE='Z Time Series', xrange=[0,n_elements(z_samp)], xstyle=1
    !P.multi = 0

    !P.multi = [0, 1, 2]
    window, /free, xsize=800, ysize=2*200
    plot, reform(ac_ratios), $
          TITLE='Acceptance Ratio', xstyle=1, $
          xrange=[0L,n_elements(ac_ratios)], $
          yrange=[min(ac_ratios),max(ac_ratios)]
    plot, reform(sigmas), $
          TITLE='Sigma', xstyle=1, $
          xrange=[0L,n_elements(sigmas)], $
          yrange=[min(sigmas),max(sigmas)]    
    !P.multi = 0

    obj_destroy, mcmc
    obj_destroy, mm
    
end
