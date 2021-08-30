;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; NOTE: This is not part of MAS, this code exists as a project which did not
;       complete due to lack of time and/or interest by the lab. 
;       
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; $Id$
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;function mas_probtrack_make_seeds_3d, mask_ptr, density=density
;
;    sz_mask = size(*mask_ptr, /dimensions)
;    density = (keyword_set(density)) ? density : 1
;
;    active = where(*mask_ptr ne 0, n_active)
;    
;    active_ind = array_indices(*mask_ptr, active)
;
;    if (density eq 1) then return, active_ind
;    
;    work = bytarr(density,density,density)
;
;    seeds = 0D
;    
;    t3d, /reset, scale=double([1,1,1])/density
;    
;    for n = 0, n_elements(active_ind)/3 - 1 do begin
;    
;        vox = array_indices(work, indgen(density^2-1))
;        vox = vert_t3d(vox, /no_copy)
;        
;        vox[0,*] += active_ind[0,n]
;        vox[1,*] += active_ind[1,n]
;        vox[2,*] += active_ind[2,n]
;
;        seeds = (n eq 0) ? vox : [ [seeds], [vox] ]
;        
;    endfor    
;
;    t3d, /reset
;    
;    return, seeds
;    
;end
;
;function mas_probtrack_make_seeds_2d, ptr_seeds, density=density
;
;    density = keyword_set(density) ? float(density) : 1.0
;    
;    plane_sz = size(*ptr_seeds[1])
;    num_planes = plane_sz[0] eq 1 ? 1 : plane_sz[2]
;    pa = *ptr_seeds[0]
;    planes = *ptr_seeds[1]
;    
;    sq = 1./sqrt(density)
;    ct = 0
;    cum_inc = 0
;    for i = 0, num_planes-1 do begin
;
;       if (not pa[i]) then continue
;       
;       inc=0
;                
;       if planes[0,i] then w = 'x' else  $
;          if planes[1,i] then w = 'y' else $
;             if planes[2,i] then w = 'z' else w='o'
;
;       seedsp = seeder(pa[i], w, sq)
;       sz_seedsp = size(seedsp)
;       inc = sz_seedsp[1]
;
;       if cum_inc eq 0 then begin
;          seeds = seedsp 
;       endif else begin
;          seeds = [seeds, seedsp]
;       endelse
;
;       cum_inc = cum_inc+inc
;
;    endfor
;    
;    ptr_free, pa
;    
;    return, transpose(seeds)
;
;end

;
;function mas_probtrack_get_seeds_from_ftroi, res_factor, density=density, $
;              roi_objects=roi_objects_new
;
;    common share_fibers_roiobj, roi_objects
;    common scan_data
;    
;    res_x = [res_factor[0]]
;    res_y = [res_factor[1]]
;    res_z = [res_factor[2]]
;
;    density = keyword_set(density) ? float(density) : 1.
;    
;    roi_objects = ptrarr(1) ;; Part of fibers common block
;
;    ptr_seeds = ft_roi_seeder(project.dataarray[project.ci].frac_ani, $
;                               res_x, res_y, res_z, 0)
;    
;    ;; roi_objects is modified by ft_roi_seeder to contain the
;    ;; roi pointers
;    
;    num_rois = n_elements(roi_objects)
;    num_rois_valid = total(ptr_valid(roi_objects))
;    
;    if (num_rois_valid eq 0) then return, 0
;    
;    tmp_rois = ptrarr(num_rois_valid)
;
;    ct = 0
;    for i=0, num_rois-1 do begin
;        if ptr_valid(roi_objects[i]) then begin
;            tmp_rois[ct++] = roi_objects[i]
;        endif
;    endfor
;    roi_objects = tmp_rois
;    roi_objects_new = roi_objects
;    
;    return, ptr_seeds
;
;end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_get_maxima_yy
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Editing Information:
;;

;pro mas_probtrack_get_maxima_yy, $
;   Pr_in, $
;   vlist_restart, $
;   n_rsverts, $
;   pow_tol, $
;   min_probability=min_probability, $
;   angle_sep=angle_sep, $
;   degrees=degrees, $
;   relative_dir=relative_dir, $
;   directions=directions, probs=probs, $
;   ncross=ncross, nthru=nthru, $
;   thru_idx=thru_idx, cross_idx=cross_idx, $
;   pr_lmt=pr_lmt, $
;   vs=vs, $
;   ps=ps
;
;    common probtrack_objective, Pr, vlist, pr_interp_tol
;
;    ;;st = systime(1)
;    Pr = Pr_in
;    ps = fltarr(n_rsverts)
;    vs = fltarr(3, n_rsverts)
;   
;    if (keyword_set(angle_sep)) then begin
;       if keyword_set(degrees) then begin
;          angle_cos = cos(angle_sep * !DTOR)
;       endif else begin
;          angle_cos = cos(angle_sep)
;       endelse
;    endif else begin
;       angle_cos = cos(!PI/6.)
;    endelse
;
;    if not keyword_set(min_probability) then begin
;       min_probability = 0.65
;    endif
;
;    for rv = 0, n_rsverts-1 do begin
;       
;       ;; Select a vertex from the "restart" vertices; call it "test"
;       test = reform(vlist_restart[*,rv])
;       
;       ;; find neighbors of "test" (within pr_interp_tol)
;       ;; call these the "potential coordinates"
;       pot_coords = where(abs(vlist[0,*] - test[0]) lt pr_interp_tol and $
;                          abs(vlist[1,*] - test[1]) lt pr_interp_tol and $
;                          abs(vlist[2,*] - test[2]) lt pr_interp_tol, count)
;       
;       ;; proceed in the direction of greatest increase in Pr
;       ;; assumes that the probability function is reasonably smooth
;       ;; (which is it in this case)
;       curr_pr = 0.0
;       iter = 0L
;       while (1) do begin
;          if (count eq 0) then begin
;             print, "No vertices returned"
;             break
;          endif else begin
;             
;             ;; get the probabilites associated with the pot_coords
;             ;; and choose the max out of all of them
;             possible_pr = Pr[pot_coords]
;             max_pr = max(possible_pr)
;             possible_coord = (pot_coords[where(possible_pr eq max_pr)])[0]
;             
;             ;; if this probabilty is greater than current max, save
;             ;; it + the coords
;             if (max_pr gt curr_pr) then begin
;                test = vlist[*, possible_coord]
;                curr_pr = max_pr
;                iter++
;             endif else begin
;                ;; break out of the loop when the probabilites begin
;                ;; to decrease
;                iter++
;                break
;             endelse
;          endelse
;
;          ;; now find the neighbors for the next iteration
;          pot_coords = where(abs(vlist[0,*] - test[0]) lt pr_interp_tol and $
;                             abs(vlist[1,*] - test[1]) lt pr_interp_tol and $
;                             abs(vlist[2,*] - test[2]) lt pr_interp_tol, count)
;          
;       endwhile
;
;       ;; at the end of the while loop, "test" will have the vector
;       ;; associated witht he maximum probability 
;
;       ps[rv] = max_pr
;       vs[*,rv] = test
;       
;    endfor
;
;    ;; Stage 2: Armed with the local maxima for each restart vertex,
;    ;; we have to weed out duplicates and vertices that are close
;    ;; enough to each other that they could be considered the same.
;    
;     max_ps = (where(ps eq max(ps)))[0]
;
;     dp       = matrix_multiply(vs, vs, /atranspose)
;     dirtrak  = bytarr(3,n_elements(ps))
;     register = bytarr(n_elements(ps))
;     dir      = fltarr(3,3)
;     probs    = fltarr(3)
;     dp_rel   = fltarr(3)
;     test_ver = reform(dp[max_ps, *])
;     
;     for i = 0, 2 do begin
;
;        close = (abs(test_ver) gt angle_cos and register eq 0)
;        
;        dirtrak[i,*] = close
;        
;        close_idx = where(close eq 1, nclose)
;        
;        dirs = reform(vs[*, close_idx])
;        
;        neg_dp = where(test_ver[close_idx] lt 0, neg_ct)
;        if (neg_ct gt 0) then begin
;           dirs[*,neg_dp] = - dirs[*,neg_dp]
;        endif
;        
;        prob = reform(ps[close_idx])
;        if (nclose gt 1) then begin
;           dir[*,i] = total(dirs, 2)/nclose
;           dir[*,i] /= sqrt(total(dir[*,i] * dir[*,i]))
;           probs[i] = total(prob)/nclose
;           dirs = reform(dir[*,i])
;        endif else begin
;           dir[*,i] = dirs
;           probs[i] = prob
;        endelse
;        
;        dp_rel[i] = abs(total(dirs * relative_dir))
;        
;        register += close
;        far = where(register eq 0, nfar)
;        if (nfar eq 0) then begin
;           i++
;           break
;        end
;        test_ver = reform(dp[far[0], *])
;
;     endfor
;     
;     thru_indicator  = (dp_rel eq max(dp_rel))
;     cross_indicator = (dp_rel lt max(dp_rel))
;
;     invalid = where(probs lt min_probability, inv_count)
;     if (inv_count gt 0) then begin
;        thru_indicator[invalid] = 0
;        cross_indicator[invalid] = 0
;     endif
;     
;     thru_idx   = where(thru_indicator eq 1, nthru)
;     cross_idx  = where(cross_indicator eq 1, ncross)
;     directions = dir
;     probs      = probs
;
;     ;;print, "Time:", systime(1) - st
;end

pro mas_probtrack_bootstrap_data, method_spec, $
        threshold=threshold, $
        nbootstrap=nbootstrap
        
    common scan_data
    common probtrack_objective, Pr, vlist, pr_interp_tol
    
    ci = project.ci
    
    if (not obj_valid(method_spec)) then return
    
    if (n_elements(pr_interp_tol) eq 0) then pr_interp_tol = 0.250
    
    if (not keyword_set(normalize_probabilities)) then begin
        normalize_probabilities = 0
    endif
    
    if (not keyword_set(nbootstrap)) then begin
       nbootstrap = 100
    endif

    threshold  = n_elements(threshold) ne 0  ? float(threshold) : 1.5
    
    data = project.dataarray[ci].state1
    
    data_sz = size((*data), /dimensions)
    
    fdim = data_sz[0]
    pdim = data_sz[1]
    sdim = data_sz[2]
    adim = data_sz[3]
       
    mas_tessellator_make_icos, level=5, vertlist=tmp;;, polylist=polylist
    vlist = tmp[*,where(tmp[2,*] ge 0)]
    
    method_spec->setReconstructionDirections, vlist
    tx = *project.imndarray[project.ci].acq_matrix
    if (obj_class(method_spec) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
        vlist = vert_t3d(vlist, matrix=tx, /double)
    endif
    
    ;; prepare the "restart vertices" for finding max Pr
    mas_tessellator_make_icos, level=0, vertlist=tmp
    vlist_restart = tmp[*,where(tmp[2,*] ge 0)]
    
    if (obj_class(method_spec) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
        vlist_restart = vert_t3d(vlist_restart, matrix=tx, /double)
    endif
    
    n_rsverts = n_elements(vlist_restart)/3
    ;;;;;;;;;
    
    pbar = obj_new('progressbar', $
        text='Processing Directionality', $
        title='Processing', /fast_loop)
    pbar->start
    
    mask = mas_roi_get_current_mask([fdim, pdim], crop=crop_dims, /no_transform)
    
    pbar->setProperty, text='Processing Threshold...'
    max_signal = max(*data)
    pct_threshold = max_signal * threshold/100.0
    
    pbar->setProperty, text='Processing Diffusion Data'
    
    th_final_1 = ptr_new(fltarr(fdim,pdim,sdim,nbootstrap))
    ph_final_1 = ptr_new(fltarr(fdim,pdim,sdim,nbootstrap))
    th_final_2 = ptr_new(fltarr(fdim,pdim,sdim,nbootstrap))
    ph_final_2 = ptr_new(fltarr(fdim,pdim,sdim,nbootstrap))
    
    anis_ptr = project.dataarray[ci].frac_ani
    
    iter = 0.
    n  = n_elements(vlist)/3
    n_used_max = 0L
    n_used_multi = 0L
    st = systime(1)
    for ss = 0, sdim-1 do begin
        print, "Starting Slice", ss
        for rr = 0, fdim-1 do begin
            iter++
            for pp = 0, pdim-1 do begin
            
                voxel_data = reform((*data)[rr,pp,ss,*])
                
                if (total(voxel_data) eq 0) then continue
                ;;if (voxel_data[0] lt pct_threshold) then continue
                if (ptr_valid(mask)) then begin
                    if (*mask)[rr,pp] eq 0 then continue
                endif
                
                Pr = method_spec->getDisplacementProbability(data=voxel_data)
                
                nans = finite(Pr, /nan)
                if (total(nans) ne 0) then continue
                seed = systime(1)
                void = method_spec->makeBootstrap(data=voxel_data)
                
                for b = 0, nbootstrap-1 do begin
                
                    Pr = method_spec->sampleBootstrap(seed, data=voxel_data)
                    mx = max(Pr, min=mn)
                    if (mn ne mx) then Pr = (Pr - mn)/(mx - mn)
                    
                    if (0 && 1 or (*anis_ptr)[rr,pp,ss] gt 0.5) then begin
                        n_used_max++
                        max_Pr = 1.0;; max(Pr)
                        max_Pr_ind = where(Pr eq max_Pr, count)
                        if (count eq 0) then begin
                            print, "OOF!"
                            break
                        endif
                        max_dir = vlist[*,max_Pr_ind[0]]
                    
                        vec_sp = cv_coord(from_rect=max_dir, /to_sphere) ;, /degrees)
                        if (vec_sp[0] lt 0) then vec_sp[0] += 2*!PI
                        (*ph_final_1)[rr,pp,ss,b] = vec_sp[0]
                        (*th_final_1)[rr,pp,ss,b] = vec_sp[1]
                        
                    endif else begin
                        ;;;
                        n_used_multi++
                        mas_hardi_tractography_get_maxima, Pr, vlist_restart, 0.65, $
                            vs=vs, ps=ps, nvs=nvs
                            
                        prim_dir = reform(vs[*,0])
                        ;; theta is polar angle, phi is equatorial angle
                        vec_sp = cv_coord(from_rect=prim_dir, /to_sphere) ;, /degrees)
                        if (vec_sp[0] lt 0) then vec_sp[0] += 2*!PI
                        (*ph_final_1)[rr,pp,ss,b] = vec_sp[0]
                        (*th_final_1)[rr,pp,ss,b] = !DPI/2 - vec_sp[1]
                        
                        if (nvs gt 1) then begin
                            secd_dir = reform(vs[*,1])
                            vec_sp = cv_coord(from_rect=secd_dir, /to_sphere) ;, /degrees)
                            if (vec_sp[0] lt 0) then vec_sp[0] += 2*!PI
                            (*ph_final_2)[rr,pp,ss,b] = vec_sp[0]
                            (*th_final_2)[rr,pp,ss,b] = !DPI/2 - vec_sp[1]
                        endif else begin
                            (*ph_final_2)[rr,pp,ss,b] = vec_sp[0]
                            (*th_final_2)[rr,pp,ss,b] = !DPI/2 - vec_sp[1]
                        endelse
                        
                    endelse
;;;
;;;
;;                     mas_probtrack_get_maxima_yy, $
;;                        Pr, vlist_restart, n_rsverts, 0.01, $
;;                        relative_dir=max_dir, $
;;                        directions=directions, $
;;                        thru_idx=thru_idx, cross_idx=cross_idx, $
;;                        ncross=ncross, nthru=nthru, probs=probs, $
;;                        ps=ps, vs=vs
                    
;;                     pr_total = 0.

;;                     if (nthru gt 0) then begin
;;                        prim_dir = reform(directions[*,thru_idx])
;;                        if (prim_dir[2] lt 0) then prim_dir *= -1.0
;;                        vec_sp = cv_coord(from_rect=prim_dir, /to_sphere) ;, /degrees)
;;                        if (vec_sp[0] lt 0) then vec_sp[0] += 2*!PI
;;                        (*ph_final_1)[rr,pp,ss,b] = vec_sp[0]
;;                        (*th_final_1)[rr,pp,ss,b] = vec_sp[1]
;;                     endif
                    
;;                     if (ncross gt 0) then begin
;;                        secd_dir = reform(directions[*,cross_idx[0]])
;;                        if (secd_dir[2] lt 0) then secd_dir *= -1.0
;;                        vec_sp = cv_coord(from_rect=secd_dir, /to_sphere) ;, /degrees)
;;                        if (vec_sp[0] lt 0) then vec_sp[0] += 2*!PI
;;                        (*ph_final_2)[rr,pp,ss,b] = vec_sp[0]
;;                        (*th_final_2)[rr,pp,ss,b] = vec_sp[1]
;;                     endif
                    
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                    
                    
                 endfor
                
             endfor
            
            if rr mod 2 eq 0 then begin
                pbar->update, iter/(fdim*sdim) * 100.0
                if (pbar->checkCancel()) then begin
                    pbar->destroy
                    anis = 0
                    goto, CLEANUP
                endif
            endif
            
        endfor
        print, "Done Slice (max, multi)", ss, n_used_max, n_used_multi
    endfor
    
    pbar->destroy
    
    print, "Used Max:", n_used_max
    print, "Used_multi:", n_used_multi
    print, "Time: ", systime(1)-st
    
    path = project.current_path
    print, "Exporting to: "+path
    mas_export_nifti, data_ptr=th_final_1, file_name=path+'/qbi_sing_th1_final.nii'
    mas_export_nifti, data_ptr=ph_final_1, file_name=path+'/qbi_sing_ph1_final.nii'
    mas_export_nifti, data_ptr=th_final_2, file_name=path+'/qbi_sing_th2_final.nii'
    mas_export_nifti, data_ptr=ph_final_2, file_name=path+'/qbi_sing_ph2_final.nii'
    
    CLEANUP:
    
    ptr_free, th_final_1
    ptr_free, ph_final_1
    ptr_free, th_final_2
    ptr_free, ph_final_2    
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_show_voxel
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; ** THIS IS A WORK IN PROGRESS, THE FINAL FORM 
;;    HAS NOT BEEN DETERMINED **
;;
;; Purpose of subroutine:
;;
;; Editing Information:
;;
pro mas_probtrack_show_voxel, reco_obj, coord, $
        thru_angle_lmt=thru_angle_lmt, $
        branch_angle_lmt=branch_angle_lmt, $
        pr_int_tol=pr_int_tol, $
        powell_tol=powell_tol, $
        restart_tess=restart_tess, $
        use_state_ptr=use_state_ptr
        
    common scan_data
    common probtrack_objective, Pr, vlist, pr_interp_tol
    
    ci = project.ci
    
    if (not keyword_set(use_state_ptr)) then begin
    
        ;; arguments that affect fiber tracking
        pr_interp_tol  = keyword_set(pr_int_tol)      ? float(pr_int_tol)     : 0.085
        powell_tol     = arg_present(powell_tol)      ? float(powell_tol)     : 0.001
        thru_ang_lmt   = keyword_set(thru_ang_lmt)    ? float(thru_ang_lmt)   : 45.0
        branch_ang_lmt = keyword_set(branch_ang_lmt)  ? float(branch_ang_lmt) : 45.0
        prob_surf_tess = 32
        prob_branch_tol= 0.75
        
    endif else begin
    
        pt_state = project.procpramarray[ci].probtrack_state
        
        if (not ptr_valid(pt_state)) then begin
            junk = dialog_message('Cannot find state pointer.', /error, /center)
            return
        endif
        
        pr_interp_tol  = (*pt_state).prob_interp_tol
        powell_tol     = (*pt_state).prob_opt_tol
        thru_ang_lmt   = (*pt_state).thru_ang_lmt
        branch_ang_lmt = (*pt_state).branch_ang_lmt
        prob_branch_tol= (*pt_state).prob_branch_tol
        prob_surf_tess = (*pt_state).prob_surf_tess
        
    endelse
        
    mas_tessellator_make_icos, level=(prob_surf_tess<5)>0, $
        vertlist=vlist, $
        polylist=polylist
        
    reco_obj->setReconstructionDirections, vlist
    
    restart_tess = (keyword_set(restart_tess)) ? floor(restart_tess) : 1
    
    cos_thru_lmt   = cos(thru_ang_lmt*!DTOR)
    cos_branch_lmt = cos(branch_ang_lmt*!DTOR)
    step_size      = 2.0
    
    ;;;
    ;junk = get_orient_rgb_axis(tx_matrix=tx)
    ;tx = *project.imndarray[project.ci].acq_matrix
    ;if (obj_class(reco_obj) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
    ;    vlist = vert_t3d(vlist, matrix=tx, /double)
    ;endif
    
    ;;;
    mas_tessellator_make_icos, level=restart_tess, vertlist=vlist_restart
    
    ;if (obj_class(reco_obj) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
    ;    vlist_restart = vert_t3d(vlist_restart, matrix=tx, /double)
    ;endif
    
    n_rsverts = n_elements(vlist_restart)/3
    
    ;;;
    data = reform((*project.dataarray[ci].state1)[coord[0], coord[1], coord[2], *])
    
    Pr = reco_obj->getDisplacementProbability(data=data)
    Pr = (Pr - min(Pr))/(max(Pr)-min(Pr))
    
    if (total(finite(Pr, /nan)) gt 0) then begin
        junk = dialog_message('Voxel location has NAN or INF probability.', /error, /center)
        return
    endif
    
    max = where(Pr eq max(Pr), count)
    max_dir = vlist[*,max[0]]

    mas_hardi_tractography_get_maxima, Pr, vlist_restart, prob_branch_tol, $
                                 vs=vs, ps=ps, nvs=nvs, maxima_mask=maxima_mask
    directions=vs        
    tmp = vlist
    
    print, 'surface_area_before: ' + $
        strcompress(string(mesh_surfacearea(tmp, polylist)), /remove_all)
        
    tmp[0,*] = tmp[0,*]*step_size*0.5*Pr
    tmp[1,*] = tmp[1,*]*step_size*0.5*Pr
    tmp[2,*] = tmp[2,*]*step_size*0.5*Pr
        
    print, 'surface_area_after: ' + $
        strcompress(string(mesh_surfacearea(tmp, polylist)), /remove_all)
        
    psurface = obj_new('idlgrpolygon', tmp, $
        name=strcompress("Anis=", /remove_all),$
        polygons=polylist, alpha_channel=0.65, $
        color=[180,180,180], style=1, shading=0)
        
    tmp = vs
    tmp[0,*] = tmp[0,*]*ps*step_size*0.5
    tmp[1,*] = tmp[1,*]*ps*step_size*0.5
    tmp[2,*] = tmp[2,*]*ps*step_size*0.5
        
    nbootstrap = 1500
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if (nbootstrap gt 0) then begin
        seed = systime(1)
        void = reco_obj->makeBootstrap(data=data)
        
        for b = 0, nbootstrap-1 do begin
        
            Pr = reco_obj->sampleBootstrap(seed, data=data)
            mx = max(Pr, min=mn)
            if (mn ne mx) then Pr = (Pr - mn)/(mx - mn)
            max = where(Pr eq max(Pr), count)
            max_dir = vlist[*,max[0]]

            mas_hardi_tractography_get_maxima, Pr, vlist_restart, prob_branch_tol, $
                                         vs=vs, ps=ps, nvs=nvs, maxima_mask=maxima_mask
            directions=vs        
            tmp = vlist
            
;            mas_probtrack_get_maxima_yy, $
;                Pr, vlist_restart, n_rsverts, powell_tol, $
;                min_probability=prob_branch_tol, $
;                relative_dir=max_dir, $
;                directions=directions, $
;                thru_idx=thru_idx, cross_idx=cross_idx, $
;                ncross=ncross, nthru=nthru, $
;                ps=ps, vs=vs         
            ;;;
            odirs = objarr(n_elements(directions)/3)
            
            for j = 0, nvs-1 do begin
                d = reform(directions[*,j])
                col = [0,0,255]
                odirs[j] = obj_new('idlgrpolyline', $
                    [ [-(d*step_size*0.65)], [(d*step_size*0.65)] ], $
                    color = col, $
                    thick=1)
            endfor
            
            ooodirs = (n_elements(ooodirs) eq 0) ? odirs : [ooodirs, odirs ]
            
        endfor
        
        odirs = ooodirs
        
    endif else begin
    
        odirs = objarr(n_elements(directions)/3)
        
        for j = 0, n_elements(directions)/3-1 do begin ;;;(nthru+ncross)-1 do begin
            d = reform(directions[*,j])
            ;;if j eq thru_idx then col = [255,0,0] else col = [0,0,255]
            col = [0,0,255]
            odirs[j] = obj_new('idlgrpolyline', $
                [ [-(d*step_size*0.65)], [(d*step_size*0.65)] ], $
                color = col, $
                thick=1)
        endfor
    endelse

    xobjview, [ psurface, odirs[where(obj_valid(odirs))] ], background=[0,0,0], renderer=0
    
    HEAP_GC
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_track
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;; 
;; ** THIS IS A WORK IN PROGRESS, THE FINAL FORM 
;;    HAS NOT BEEN DETERMINED **
;;
;; Purpose of subroutine:
;;
;; Editing Information:
;;
pro mas_probtrack_track, track_method_spec, use_state_ptr=use_state_ptr, $
                         anis_halt_thr=anis_halt_thr, $
                         anis_cross_thr=anis_cross_thr, $
                         anis_scl=anis_scl, $
                         thru_ang_lmt=thru_ang_lmt, $
                         branch_ang_lmt=branch_ang_lmt, $
                         seed_density=seed_dentsity, $
                         step_size=step_size,$
                         use_ft_roi=use_ft_roi,$
                         use_roi_file=use_roi_file, $
                         iter_max=iter_max, $
                         pr_int_tol=pr_int_tol, $
                         powell_tol=powell_tol, $
                         cross_depth=cross_depth, $
                         branch_per_vxl=branch_per_vxl, $
                         cross_show=cross_show, $
                         output_dir=output_dir, $
                         trs_filename=trs_filename, $
                         accounting=accounting
                         
    return
end
;pro mas_probtrack_track, track_method_spec, use_state_ptr=use_state_ptr, $
;                         anis_halt_thr=anis_halt_thr, $
;                         anis_cross_thr=anis_cross_thr, $
;                         anis_scl=anis_scl, $
;                         thru_ang_lmt=thru_ang_lmt, $
;                         branch_ang_lmt=branch_ang_lmt, $
;                         seed_density=seed_dentsity, $
;                         step_size=step_size,$
;                         use_ft_roi=use_ft_roi,$
;                         use_roi_file=use_roi_file, $
;                         iter_max=iter_max, $
;                         pr_int_tol=pr_int_tol, $
;                         powell_tol=powell_tol, $
;                         cross_depth=cross_depth, $
;                         branch_per_vxl=branch_per_vxl, $
;                         cross_show=cross_show, $
;                         output_dir=output_dir, $
;                         trs_filename=trs_filename, $
;                         accounting=accounting
;
;    common scan_data
;    common share_fibers_roiobj, roi_objects
;    common probtrack_objective, Pr, vlist, pr_interp_tol
;    forward_function mas_tricubic_4d
;
;    if (not obj_valid(track_method_spec)) then return
;
;    ci = project.ci
;
;    if (not keyword_set(use_state_ptr)) then begin
;
;        ;; arguments that affect fiber tracking
;        pr_interp_tol  = keyword_set(pr_int_tol)      ? float(pr_int_tol)     : 0.085
;        pr_thresh      = keyword_set(prob_branch_thr) ? float(prob_branch_thr): 0.75
;        powell_tol     = arg_present(powell_tol)      ? float(powell_tol)     : 0.001
;        anis_scl       = keyword_set(anis_scl)        ? float(anis_scl)       : 1.0
;        anis_halt_thr  = arg_present(anis_halt_thr)   ? float(anis_halt_thr)  : 0.17
;        anis_cross_thr = arg_present(anis_cross_thr)  ? float(anis_cross_thr) : 0.3
;        vdim = 24
;        thru_ang_lmt   = keyword_set(thru_ang_lmt)    ? float(thru_ang_lmt)   : 45.0
;        branch_ang_lmt = keyword_set(branch_ang_lmt)  ? float(branch_ang_lmt) : 45.0
;        
;        step_size      = keyword_set(step_size)       ? float(step_size)      : 0.3
;        seed_density   = keyword_set(seed_density)    ? fix(seed_density)     : 1
;        max_level      = keyword_set(cross_depth)     ? cross_depth           : 0
;        cross_show     = keyword_set(cross_show)      ? 1                     : 0
;        branch_per_vxl = keyword_set(branch_per_vxl)  ? branch_per_vxl        : 5
;        iter_max       = keyword_set(iter_max)        ? iter_max              : 1000
;        iter_min       = keyword_set(iter_min)        ? iter_min              : 20
;        output_stream  = 1
;        output_surface = 0
;        output_voxeldata = 0
;
;    endif else begin
;
;        pt_state = project.procpramarray[ci].probtrack_state
;
;        if (not ptr_valid(pt_state)) then begin
;            junk = dialog_message('Cannot find state pointer.', /error)
;            return
;        endif
;
;        pr_interp_tol  = (*pt_state).prob_interp_tol
;        powell_tol     = (*pt_state).prob_opt_tol
;        pr_thresh      = (*pt_state).prob_branch_tol
;        vdim           = (*pt_state).prob_surf_tess
;
;        anis_scl       = 0.175; 1.0 
;        anis_halt_thr  = (*pt_state).anis_halt_thr
;        anis_cross_thr = (*pt_state).anis_branch_thr
;        
;        thru_ang_lmt   = (*pt_state).thru_ang_lmt
;        branch_ang_lmt = (*pt_state).branch_ang_lmt
;        
;        step_size      = (*pt_state).step_size
;        seed_density   = (*pt_state).seed_density
;        max_level      = (*pt_state).max_branch_lvl
;        branch_per_vxl = (*pt_state).max_branch_vxl
;        cross_show     = (*pt_state).show_crossings
;
;        iter_max       = (*pt_state).max_length
;        iter_min       = 0.; (*pt_state).min_length
;        output_dir     = (*pt_state).output_directory
;        trs_filename   = (*pt_state).output_filename
;        output_surface = (*pt_state).output_surface
;        output_stream  = (*pt_state).output_stream
;        output_voxeldata = (*pt_state).output_voxeldata
;        seed_mask_files = (*pt_state).roi_filelist_3d
;        seed_mask_types = (*pt_state).roi_filelist_3d_types
;        if ((*pt_state).res_crossings eq 0) then begin
;            max_level = 0
;            anis_cross_thr = 0.0
;        endif
;
;    endelse
;
;    ;; need to apply data-set specific transformations to the
;    ;; vertices
;    ;;mas_tessellator_make_sphere, vdim=vdim, /normalize, $
;    ;;                            vertlist=vlist, polylist=polylist
;    mas_tessellator_make_icos, level=3, vertlist=vlist, polylist=polylist
;
;    track_method_spec->setReconstructionDirections, vlist
;    junk = get_orient_rgb_axis(tx_matrix=tx)
;    tx = *project.imndarray[project.ci].acq_matrix # tx
;    if (obj_class(track_method_spec) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
;       vlist = vert_t3d(vlist, matrix=tx, /double)
;    endif
;
;    ;; prepare the "restart vertices" for finding max Pr
;    mas_tessellator_make_icos, level=0, vertlist=vlist_restart
;    if (obj_class(track_method_spec) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
;       vlist_restart = vert_t3d(vlist_restart, matrix=tx, /double)
;    endif
;    
;    n_rsverts = n_elements(vlist_restart)/3
;
;    ;; data set size and related parameters
;    sz   = size(*project.dataarray[project.ci].state1, /dimensions)
;    xdim = sz[0]
;    ydim = sz[1]
;    zdim = sz[2]
;
;    interp_factors = [ project.procpramarray[ci].freq_interp, $
;                      project.procpramarray[ci].phase_interp, $
;                      project.procpramarray[ci].slice_interp ]
;
;    voxel_dims = [ project.imndarray[ci].f_voxsz, $
;                   project.imndarray[ci].p_voxsz, $
;                   project.imndarray[ci].s_voxsz ]
;
;    res_factors = voxel_dims/interp_factors    
;    res_x = res_factors[0]
;    res_y = res_factors[1]
;    res_z = res_factors[2]
;
;    paths_to_rawdata = diag_matrix(double([1.0/interp_factors, 1.0]))
;    paths_to_rawdata = diag_matrix(replicate(1D,4))
;    
;    sl   = project.procpramarray[project.ci].sdim_start
;    n    = n_elements(vlist)/3
;    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; ROI Selection
;
;    terminator_mask = ptr_new(bytarr(xdim,ydim,zdim))
;
;    if keyword_set(trs_filename) then begin
;        junk = stregex(trs_filename, '^(.+)(\.(trs|ROI))*$')
;        if (junk ne -1) then begin
;            junk = stregex(trs_filename, '^(.+)(\.(trs|ROI))*$', /subexpr, /extract)
;            trs_filename = junk[1]
;        endif
;    endif else begin
;        trs_filename = 'tracts_prob'
;    endelse
;
;    if not keyword_set(output_dir) then begin
;        output_dir = project.current_path
;    endif
;    junk = strmid(output_dir, strlen(output_dir)-1)
;    if (junk eq (get_dir_slash())[0]) then begin
;        output_dir = strmid(output_dir, 0, strlen(output_dir)-1)
;    endif
;
;    ;; this omodel gets passed to the view_fibers routines and can
;    ;; contain anything.
;    omodel = obj_new('idlgrmodel', depth_test_disable=2)
;    
;    ;; we need to prepare the seed points using either an existing ROI
;    ;; file or generate new ROIs
;    have_roi = 0
;    
;    if (keyword_set(use_roi_file)) then begin
;
;        filename = dialog_pickfile(filter='*.ROI', path=output_dir)
;        if not file_test(filename, /read) then return
;
;        output_dir = file_dirname(filename)
;        trs_filename = file_basename(filename, '.ROI')
;
;        restore, filename
;        
;        if (n_elements(ptr_seeds) ne 0) then begin
;            roi_sdpts = mas_probtrack_make_seeds_2d(ptr_seeds, density=seed_density)
;            ptr_free, ptr_seeds
;        endif else begin
;            obj_destroy, track_method_spec
;            return
;        endelse
;        
;        have_roi = 1
;        
;    endif else if (keyword_set(seed_mask_files)) then begin
;        
;        roi_mask = ptr_new(bytarr(xdim,ydim,zdim))
;        roi_num = 0
;        for i = 0, n_elements(seed_mask_files)-1 do begin
;        
;            if (seed_mask_files[i] eq '' or not file_test(seed_mask_files[i], /read)) then continue
;            
;            nif = mas_read_nifti(nifti_filename=seed_mask_files[i], read_status=rs)
;            ;*nif.voxel_data = byte(*nif.voxel_data)
;            
;            for lab = 1, 10 do begin
;
;                labeled_voxels = where(*nif.voxel_data eq lab, n_labeled)
;                if (n_labeled eq 0) then continue
;                
;                tmproi = array_indices(*nif.voxel_data, labeled_voxels)
;                roi_objects = roi_num++ eq 0 ? ptr_new(tmproi) : [ roi_objects, ptr_new(tmproi) ]
;
;                if (seed_mask_types[i] eq 0) then begin
;                   *roi_mask = (*roi_mask or *nif.voxel_data)
;                endif else if (seed_mask_types[i] eq 2) then begin
;                   print, 'terminator'
;                   *terminator_mask = (*terminator_mask or *nif.voxel_data)
;                endif
;
;            endfor
;            ptr_free, nif.voxel_data
;            
;        endfor
;        
;        mas_indices = where(*roi_mask ne 0, n_voxels)
;        if (n_voxels ne 0) then begin
;            roi_sdpts = mas_probtrack_make_seeds_3d(roi_mask, density=seed_density)            
;            save, roi_objects, filename=output_dir+(get_dir_slash())[0]+trs_filename+'.ROI'
;            have_roi = 1
;        endif
;        
;    endif
;    
;    if (have_roi eq 0) then begin
;        ptr_seeds = mas_probtrack_get_seeds_from_ftroi([res_x,res_y,res_z], $
;                         roi_objects=roi_objects)
;                         
;        ;; int type returned means no roi selected
;
;        if (size(ptr_seeds, /type) eq 2) then begin
;            void = dialog_message('No ROI selected.', /center, /error)
;            return
;        endif
;        
;        save, ptr_seeds, roi_objects, $
;              filename=output_dir+(get_dir_slash())[0]+trs_filename+'.ROI'
;        
;        roi_sdpts = mas_probtrack_make_seeds_2d(ptr_seeds, density=seed_density)
;        
;        ;;ptr_free, ptr_seeds
;        
;    endif
;
;    n_roi_sdpts = (size(roi_sdpts, /dimensions))[1]
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Initialize Tractography Data Structures
;    
;    ;; these all grow dynamically
;    start_points  = ulon64arr(1) ;; Seed Points
;    mid_points    = ulon64arr(1) ;; Seed Points
;    branch_levels = lonarr(max_level+1) 
;    pts_x         = fltarr(1) ;; X coordinates of streamlines
;    pts_y         = fltarr(1) ;; Y      "       "      "
;    pts_z         = fltarr(1) ;; Z      "       "      "
;    sdpts         = fltarr(3,1) ;; list of ROI seed points
;    sdvel         = fltarr(3,1) ;; list of initial directions for seed points
;    iter          = 0
;
;    res_fac_x=min([res_x, res_y, res_z])/res_x
;    res_fac_y=min([res_x, res_y, res_z])/res_y
;    res_fac_z=min([res_x, res_y, res_z])/res_z
;    
;    res_fac = [res_fac_x, res_fac_y, res_fac_z]*step_size
;
;    cos_branch_lmt = cos(branch_ang_lmt*!dtor)
;    cos_thru_lmt   = cos(thru_ang_lmt*!dtor)
;    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Initialize seed points and look for any seed points that have
;;; crossings
;
;    for j = 0L, n_roi_sdpts-1 do begin
;
;         pt = roi_sdpts[*,j]
;         
;         ;;data = reform(mas_odf_interpolate(pt))
;         ;;data = mas_4d_data_average(pt)
;         data = mas_tricubic_4d(pt)
;         ;;data = mas_cubic_convolution(pt)
;         
;         Pr = track_method_spec->getDisplacementProbability(data=data)
;         max_pr = max(Pr, min=min_pr)
;         Pr = (Pr - min_pr)/(max_pr - min_pr)
;
;         if total(finite(Pr,/nan)) then begin
;            terminate = 7
;            continue
;         endif
;
;         anis = (*project.dataarray[project.ci].frac_ani)[pt[0],pt[1],pt[2]]
;        
;         ;; get the main direction
;         max = where(Pr eq max(Pr), count)
;         max_dir = vlist[*,max[0]]
;
;         if (anis lt anis_halt_thr) then continue
;         if (anis le anis_cross_thr) then begin
;            
;            mas_probtrack_get_maxima_yy, $
;               Pr, vlist_restart, n_rsverts, powell_tol, $
;               min_probability=pr_thresh, $
;               relative_dir=max_dir, $
;               directions=directions, $
;               thru_idx=thru_idx, cross_idx=cross_idx, $
;               ncross=cross, nthru=thru, $
;               ps=ps, vs=vs
;
;            if (thru gt 0) then begin
;               sdpts = [ [sdpts],  [pt] ]
;               sdvel = [ [sdvel],  [reform(directions[*, thru_idx])] ]
;            endif
;               
;            if (cross gt 0) then begin
;               sdpts = [ [sdpts],  [pt] ]
;               sdvel = [ [sdvel],  [reform(directions[*, cross_idx[0]])] ]
;            endif 
;
;            if (cross gt 1) then begin
;               sdpts = [ [sdpts],  [pt] ]
;               sdvel = [ [sdvel],  [reform(directions[*, cross_idx[1]])] ]
;            endif 
;                        
;         endif else begin
;             sdpts = [ [sdpts],  [pt] ]
;             sdvel = [ [sdvel],  [max_dir] ]
;         endelse
;
;      endfor
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Begin Tracking
;
;    start_time  = systime(1)
;    seedpt_num  = 0L
;    n_seedpts   = long((size(sdpts, /dimensions))[1])
;    sdpts       = reform(sdpts[*, 1:n_seedpts-1])
;    sdvel       = reform(sdvel[*, 1:n_seedpts-1])
;    n_seedpts--
;    orig_ct     = n_seedpts
;
;    data_size  = size(*project.dataarray[ci].state1, /dimensions)
;    accounting = bytarr(data_size[0],data_size[1], data_size[2],max_level+1)
;    cr_accounting = bytarr(data_size[0],data_size[1], data_size[2],max_level+1)
;
;    level      = 0
;    nfibers    = 0
;    branch_levels[0] = n_seedpts
;    curr_pt    = { pos: fltarr(3), vel: fltarr(3) }
;    parent_seeds = 0
;    
;    pbar = obj_new('progressbar', text='HARDI Fiber Tracking...', /fast_loop, /nocancel)
;    pbar->start
;    
;    while seedpt_num lt n_seedpts-1 do begin
;
;        curr_pt.pos = sdpts[*,seedpt_num]
;        curr_pt.vel = sdvel[*,seedpt_num]
;        pts       = [curr_pt.pos]
;        terminate = 0
;        fiter     = 0
;        iter      = 0
;        ncrossings= 0 
;
;        pts_fwd_x = fltarr(1)
;        pts_fwd_y = fltarr(1)
;        pts_fwd_z = fltarr(1)
;
;        txd_pt = ([curr_pt.pos, 1.] # paths_to_rawdata)[0:2]
;        pts_fwd_x[0] = txd_pt[0] ;curr_pt.pos[0]
;        pts_fwd_y[0] = txd_pt[1] ;curr_pt.pos[1]
;        pts_fwd_z[0] = txd_pt[2] ;curr_pt.pos[2]
;
;        ;; Keep track of branch levels
;        if seedpt_num ge orig_ct then begin
;            level++
;            orig_ct = n_seedpts
;            branch_levels[level] = orig_ct
;        endif
;
;        ;; prevents looking for crossings while on the last depth
;        ;; level. commented out because i don't think it is appropriate
;        ;if level eq max_level then anis_cross_thr = 0.0
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;        cross = 0
;        thru  = 0
;        cross_indicator = bytarr(n_rsverts)
;        thru_indicator = bytarr(n_rsverts)
;
;        ;; Forward
;        while (terminate eq 0) do begin
;            
;            thru               = 0
;            cross              = 0
;
;            ;;data = reform(mas_odf_interpolate(curr_pt.pos))
;            ;;data = mas_4d_data_average(curr_pt.pos)
;            data = mas_tricubic_4d(curr_pt.pos)
;            ;;data = mas_cubic_convolution(curr_pt.pos)
;            
;            Pr = track_method_spec->getDisplacementProbability(data=data)
;            Pr = (Pr - min(Pr))/(max(Pr)-min(Pr))
;
;            if total(finite(Pr,/nan)) then begin
;               terminate = 7
;               Pr[*] = 1
;            endif
;            
;            anis = (*project.dataarray[project.ci].frac_ani)[curr_pt.pos[0], $
;                                                             curr_pt.pos[1], $
;                                                             curr_pt.pos[2]]
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;            if (anis le anis_halt_thr) then begin
;
;                 terminate = 2
;
;             endif else if (anis le anis_cross_thr and fiter ne 0) then begin
;
;                ;; resolve crossings
;                mas_probtrack_get_maxima_yy, $
;                   Pr, vlist_restart, n_rsverts, powell_tol, $
;                   min_probability=pr_thresh, $
;                   relative_dir=last_dir, $
;                   directions=directions, $
;                   thru_idx=thru_idx, cross_idx=cross_idx, $
;                   ncross=cross, nthru=thru, $
;                   ps=ps, vs=vs
;                
;                if (thru eq 0) then begin
;                   thru = cross
;                   cross = 0
;                   thru_idx = cross_idx
;                endif
;
;                direction = reform(directions[*, thru_idx[0]])
;                
;                ;; don't record any more crossings if we've reached
;                ;; the limit
;                if (level lt max_level) then begin
;
;                   ;; don't record any more crossings for this
;                   ;; voxel if we've reached that limit
;                   if (cr_accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level] lt branch_per_vxl) then begin
;
;                      if (cross gt 0) then begin
;                         sdpts = [ [sdpts],  [curr_pt.pos] ]
;                         sdvel = [ [sdvel],  [reform(directions[*, cross_idx[0]])] ]
;                         n_seedpts++
;                      endif 
;                      
;                      if (cross gt 1) then begin
;                         sdpts = [ [sdpts],  [curr_pt.pos] ]
;                         sdvel = [ [sdvel],  [reform(directions[*, cross_idx[1]])] ]
;                         n_seedpts++
;                      endif 
;
;                      cr_accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level]++
;
;                   endif else begin
;                      cross = 0
;                   endelse
;
;                endif else begin
;                   cross = 0
;                endelse
;                
;             endif else if (anis gt anis_cross_thr and fiter ne 0) then begin
;                 ;; The anisotropy is high enough to assume
;                 ;; uni-directional model.
;
;                 max = (where(Pr eq max(Pr)))[0]
;                 direction = vlist[*,max]
;
;             endif
;
;             if fiter eq 0 then begin
;                 ;; if this is the launch...
;                 direction = curr_pt.vel
;                 direction = -direction
;                 
;             endif else begin
;                 ;; Angle test
;                 ;; We need to allow for a larger cos_thru_lmt when
;                 ;; the unidirectional model is in effect. Right now,
;                 ;; the limits are the same for uni/multi directional.
;
;                 dp = total(direction*last_dir)
;                
;                 if (dp lt 0) then begin
;                     direction = -direction
;                     dp        = -dp
;                 endif
;                
;                 cos_theta = dp
;
;                 if (cos_theta lt cos_thru_lmt) then begin
;                     terminate = 1
;                     continue
;                 endif
;
;             endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;            ;;; Debug step
;            if (cross_show and cross gt 0 and level lt max_level) then begin
;                
;                tmp = vlist
;                tmp[0,*] = tmp[0,*]*Pr*res_fac[0]*0.5 + curr_pt.pos[0]
;                tmp[1,*] = tmp[1,*]*Pr*res_fac[1]*0.5 + curr_pt.pos[1]
;                tmp[2,*] = tmp[2,*]*Pr*res_fac[2]*0.5 + curr_pt.pos[2]
;                junk1 = obj_new('idlgrpolygon', tmp, $
;                                name=strcompress("Anis="+string(anis, format='(1F0.3)')+', Pos=('+$
;                                                 string(curr_pt.pos[0], format='(1F0.3)') + ',' +$
;                                                 string(curr_pt.pos[1], format='(1F0.3)') + ',' +$
;                                                 string(curr_pt.pos[2], format='(1F0.3)') + ',' +$
;                                                 ')'), $
;                                polygons=polylist, alpha_channel=0.75, $
;                                color=[180,180,180], style=1, shading=0)
;                
;                tmp = vs
;                tmp[0,*] = tmp[0,*]*ps*res_fac[0]*0.5 + curr_pt.pos[0]
;                tmp[1,*] = tmp[1,*]*ps*res_fac[1]*0.5 + curr_pt.pos[1]
;                tmp[2,*] = tmp[2,*]*ps*res_fac[2]*0.5 + curr_pt.pos[2]
;
;                ttt = obj_new('idlgrpolygon',$
;                              tmp,$
;                              style=0, $
;                              thick=3,$
;                              color=[0,255,0])
;
;                odirs = objarr(3)
;                
;                for j = 0, (thru+cross)-1 do begin
;                   d = reform(directions[*,j])
;                   if j eq thru_idx then col = [255,0,0] else col = [0,0,255]
;                   odirs[j] = obj_new('idlgrpolyline', $
;                                      [ [curr_pt.pos], [curr_pt.pos+(d*res_fac*0.5)] ], $
;                                      color = col, $
;                                      thick=3)
;                endfor
;
;                ;omodel->add, [ junk1, ttt, odirs[where(obj_valid(odirs) eq 1)] ]
;                omodel->add, [ odirs[where(obj_valid(odirs) eq 1)] ]
;
;             endif else if (0 && cross_show) then begin
;                
;                d = direction
;                odirs = obj_new('idlgrpolyline', $
;                                [ [curr_pt.pos], [curr_pt.pos+(d*res_fac*0.5)] ], $
;                                color = [255,0,0], $
;                                thick=3)
;                
;                omodel->add, [ odirs ]
;
;             endif
;
;            accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level] = $
;              accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level] + (cross gt 0) ? cross+thru : 1
;
;            curr_pt.pos += ([direction[0], direction[1], direction[2]])*(res_fac)
;
;            if (curr_pt.pos[0] le 0 or curr_pt.pos[1] le 0 or curr_pt.pos[2] le 0) then begin
;                terminate = 4
;                continue
;            endif else if (curr_pt.pos[0] ge sz[0] or $
;                             curr_pt.pos[1] ge sz[1] or $
;                             curr_pt.pos[2] ge sz[2]) then begin
;               terminate = 4
;               continue
;            endif else if ((*terminator_mask)[curr_pt.pos[0], curr_pt.pos[1], curr_pt.pos[2]] ne 0) then begin
;               terminate = 6
;            endif
;             
;            txd_pt = ([curr_pt.pos,1.] # paths_to_rawdata)[0:2] 
;            
;            pts_fwd_x = [pts_fwd_x, txd_pt[0]]
;            pts_fwd_y = [pts_fwd_y, txd_pt[1]]
;            pts_fwd_z = [pts_fwd_z, txd_pt[2]]
;
;;            pts_fwd_x = [pts_fwd_x, curr_pt.pos[0]]
;;            pts_fwd_y = [pts_fwd_y, curr_pt.pos[1]]
;;            pts_fwd_z = [pts_fwd_z, curr_pt.pos[2]]
;            
;            if (fiter++ gt iter_max) then terminate = 5
;
;            last_dir = direction
;
;        endwhile
;        
;;;        if (1 or fiter ge iter_min) then  begin
;        if (1 or fiter ge 0) then  begin
;            case terminate of
;                1: term_cond = 'ANGLE: '+strcompress(string(cos_theta)+'('+string(cos_thru_lmt)+')', /remove_all)
;                2: term_cond = 'FA_THR:'+strcompress(string(anis)     +'('+string(anis_halt_thr)+')', /remove_all)
;                3: term_cond = 'NO_MAX'
;                4: term_cond = 'LOWER_BOUND:'+string(curr_pt.pos[0])+','+string(curr_pt.pos[1])+','+string(curr_pt.pos[2])
;                5: term_cond = 'ITERATIONS:'+string(fiter)
;                6: term_cond = 'TERMINATOR!'
;                7: term_cond = 'NAN/INF PROBABILITY'
;                else: term_cond = 'UNKNOWN'
;            endcase
;            print, "Pot. Cross: ", cross, '(F)'
;            print, "Iterations: ", fiter, '(F)'
;            print, "Point     : ", seedpt_num, " of ", n_seedpts
;            print, "Level     : ", level, strcompress('('+string(max_level)+')', /remove_all)
;            print, "Next Level: ", orig_ct
;            print, "Terminate : ", term_cond
;            print, ""
;        endif
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;        curr_pt.pos = sdpts[*,seedpt_num]
;        curr_pt.vel = sdvel[*,seedpt_num]
;
;        terminate = 0
;        biter     = 0
;
;        pts_bwd_x = fltarr(1)
;        pts_bwd_y = fltarr(1)
;        pts_bwd_z = fltarr(1)
;        
;        txd_pt = ([curr_pt.pos, 1.] # paths_to_rawdata)[0:2]
;        pts_bwd_x[0] = txd_pt[0] ;curr_pt.pos[0]
;        pts_bwd_y[0] = txd_pt[1] ;curr_pt.pos[1]
;        pts_bwd_z[0] = txd_pt[2] ;curr_pt.pos[2]
;
;        ;; Backwards
;        while (terminate eq 0) do begin
;            
;            thru = 0
;            cross = 0
;
;            ;;data = reform(mas_odf_interpolate(curr_pt.pos))
;            ;;data = mas_4d_data_average(curr_pt.pos)
;            data = mas_tricubic_4d(curr_pt.pos)
;            ;;data = mas_cubic_convolution(curr_pt.pos)
;            
;            Pr = track_method_spec->getDisplacementProbability(data=data)
;            Pr = (Pr - min(Pr))/(max(Pr)-min(Pr))
;
;            if total(finite(Pr,/nan)) then begin
;               terminate = 9
;               Pr[*] = 1
;            endif
;
;            anis = (*project.dataarray[project.ci].frac_ani)[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2]]
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;            if (anis le anis_halt_thr) then begin
;
;                 terminate = 2
;
;             endif else if (anis le anis_cross_thr and biter ne 0) then begin
;
;                mas_probtrack_get_maxima_yy, $
;                   Pr, vlist_restart, n_rsverts, powell_tol, $
;                   min_probability=pr_thresh, $
;                   relative_dir=last_dir, $
;                   directions=directions, $
;                   thru_idx=thru_idx, cross_idx=cross_idx, $
;                   ncross=cross, nthru=thru, $
;                   ps=ps, vs=vs
;               
;                if (thru eq 0) then begin
;                   thru = cross
;                   cross = 0
;                   thru_idx = cross_idx
;                endif
;
;                direction = reform(directions[*, thru_idx[0]])
;                
;                if (level lt max_level) then begin
;
;                   if (cr_accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level] lt branch_per_vxl) then begin
;
;                      if (cross gt 0) then begin
;                         sdpts = [ [sdpts],  [curr_pt.pos] ]
;                         sdvel = [ [sdvel],  [reform(directions[*, cross_idx[0]])] ]
;                         n_seedpts++
;                      endif 
;                      
;                      if (cross gt 1) then begin
;                         sdpts = [ [sdpts],  [curr_pt.pos] ]
;                         sdvel = [ [sdvel],  [reform(directions[*, cross_idx[1]])] ]
;                         n_seedpts++
;                      endif 
;
;                      cr_accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level]++
;
;                   endif else begin
;                      cross = 0
;                   endelse
;
;                endif else begin
;                   cross = 0
;                endelse
;
;             endif else if (anis gt anis_cross_thr and biter ne 0) then begin
;
;                 max = (where(Pr eq max(Pr)))[0]
;                 direction = vlist[*,max]
;
;             endif
;
;             if biter eq 0 then begin
;                 ;; if this is the launch...
;                 direction = curr_pt.vel
;                 ;;direction = direction
;                 
;             endif else begin
;
;                 dp = total(direction*last_dir)
;                
;                 if (dp lt 0) then begin
;                     direction = -direction
;                     dp        = -dp
;                 endif
;                
;                 cos_theta = dp
;
;                 if (cos_theta lt cos_thru_lmt) then begin
;                     terminate = 1
;                     continue
;                 endif
;
;             endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;            ;;; Debug step
;            if (cross_show and cross gt 0 and level lt max_level) then begin
;                
;                tmp = vlist
;                tmp[0,*] = tmp[0,*]*Pr *res_fac[0]*0.5 + curr_pt.pos[0]
;                tmp[1,*] = tmp[1,*]*Pr *res_fac[1]*0.5 + curr_pt.pos[1]
;                tmp[2,*] = tmp[2,*]*Pr *res_fac[2]*0.5 + curr_pt.pos[2]
;                junk1 = obj_new('idlgrpolygon', tmp, $
;                                name=strcompress("Anis="+string(anis, format='(1F0.3)')+', Pos=('+$
;                                                 string(curr_pt.pos[0], format='(1F0.3)') + ',' +$
;                                                 string(curr_pt.pos[1], format='(1F0.3)') + ',' +$
;                                                 string(curr_pt.pos[2], format='(1F0.3)') + ',' +$
;                                                 ')'), $
;                                polygons=polylist, alpha_channel=0.75, $
;                                color=[180,180,180], style=1, shading=0)
;                                
;                tmp = vs; vs[*,where(cross_indicator eq 1)]
;                tmp[0,*] = tmp[0,*]*ps *res_fac[0]*0.5 + curr_pt.pos[0]
;                tmp[1,*] = tmp[1,*]*ps *res_fac[1]*0.5 + curr_pt.pos[1]
;                tmp[2,*] = tmp[2,*]*ps *res_fac[2]*0.5 + curr_pt.pos[2]
;
;                ttt = obj_new('idlgrpolygon',$
;                              tmp,$
;                              style=0, $
;                              thick=3,$
;                              color=[0,255,0])
;
;                odirs = objarr(3)
;                
;                for j = 0, (thru+cross)-1 do begin
;                   d = reform(directions[*,j])
;                   if j eq thru_idx then col = [255,0,0] else col = [0,0,255]
;                   odirs[j] = obj_new('idlgrpolyline', $
;                                      [ [curr_pt.pos], [curr_pt.pos+(d*res_fac*0.5)] ], $
;                                      color = col, $
;                                      thick=3)
;                endfor
;
;                ;omodel->add, [ junk1, ttt, odirs[where(obj_valid(odirs) eq 1)] ]
;                omodel->add, [ odirs[where(obj_valid(odirs) eq 1)] ]
;
;             endif else if (0 && cross_show) then begin
;                
;                d = direction
;                odirs = obj_new('idlgrpolyline', $
;                                [ [curr_pt.pos], [curr_pt.pos+(d*res_fac*0.5)] ], $
;                                color = [255,0,0], $
;                                thick=3)
;                
;                omodel->add, [ odirs ]
;
;             endif
;
;            accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level] = $
;              accounting[curr_pt.pos[0],curr_pt.pos[1],curr_pt.pos[2],level] + (cross gt 0) ? 10 : 5
;
;            curr_pt.pos += ([direction[0], direction[1], direction[2]])*(res_fac)
;
;            if (curr_pt.pos[0] le 0 or curr_pt.pos[1] le 0 or curr_pt.pos[2] le 0) then begin
;                terminate = 4
;                continue
;            endif else if (curr_pt.pos[0] ge sz[0] or curr_pt.pos[1] ge sz[1] or curr_pt.pos[2] ge sz[2]) then begin
;                terminate = 4
;                continue
;            endif else if ((*terminator_mask)[curr_pt.pos[0], curr_pt.pos[1], curr_pt.pos[2]] ne 0) then begin
;               terminate = 6
;            endif
;            
;            txd_pt = ([curr_pt.pos,1.] # paths_to_rawdata)[0:2]
;            
;            pts_bwd_x = [pts_bwd_x, txd_pt[0]]
;            pts_bwd_y = [pts_bwd_y, txd_pt[1]]
;            pts_bwd_z = [pts_bwd_z, txd_pt[2]]
; 
;;            pts_bwd_x = [pts_bwd_x, curr_pt.pos[0]]
;;            pts_bwd_y = [pts_bwd_y, curr_pt.pos[1]]
;;            pts_bwd_z = [pts_bwd_z, curr_pt.pos[2]]
;            
;            if (biter++ gt iter_max) then terminate = 5
;
;            last_dir = direction
;
;        endwhile
;        
;;;        if (1 or biter ge iter_min) then  begin
;        if (1 or biter ge 0) then  begin
;            case terminate of
;                1: term_cond = 'ANGLE: '+strcompress(string(cos_theta)+'('+string(cos_thru_lmt)+')', /remove_all)
;                2: term_cond = 'FA_THR:'+strcompress(string(anis)        +'('+string(anis_halt_thr)+')', /remove_all)
;                3: term_cond = 'NO_MAX'
;                4: term_cond = 'LOWER_BOUND:'+string(curr_pt.pos[0])+','+string(curr_pt.pos[1])+','+string(curr_pt.pos[2])
;                5: term_cond = 'ITERATIONS:'+string(biter)
;                6: term_cond = 'TERMINATOR!'
;                7: term_cond = 'NAN/INF PROBABILITY'
;                else: term_cond = 'UNKNOWN'
;            endcase
;            print, "Pot. Cross: ", cross, '(B)'
;            print, "Iterations: ", biter, '(B)'
;            print, "Point     : ", seedpt_num, " of ", n_seedpts
;            print, "Level     : ", level, strcompress('('+string(max_level)+')', /remove_all)
;            print, "Next Level: ", orig_ct
;            print, "Terminate : ", term_cond
;            print, ""
;        endif
;
;;; There is a bug in here somewhere that sometimes occurrs when
;;; iter_min > 0 -- that it why it is forced to 0 at the beginning for
;;;                 now.
;
;;;        if (n_elements(pts_fwd_x) gt max([10,iter_min]) and
;;;        n_elements(pts_bwd_x) gt max([10,iter_min]) ) then begin
;        if (0 and fiter gt min([1, iter_min]) and biter gt min([1, iter_min])) then begin
;           if seedpt_num eq 0 then begin
;              pts_x = [pts_fwd_x, pts_bwd_x]
;              pts_y = [pts_fwd_y, pts_bwd_y]
;              pts_z = [pts_fwd_z, pts_bwd_z]
;              start_points[0] = 0
;              mid_points[0] = fiter+1
;              start_points = [start_points, fiter+biter+2]
;           endif else begin
;              pts_x = [pts_x, pts_fwd_x, pts_bwd_x]
;              pts_y = [pts_y, pts_fwd_y, pts_bwd_y]
;              pts_z = [pts_z, pts_fwd_z, pts_bwd_z]
;              mid_points = [mid_points, start_points[n_elements(start_points)-1]+fiter+1]
;              start_points=[start_points, mid_points[n_elements(mid_points)-1]+biter+1]
;           endelse
;        endif else begin
;            if (seedpt_num eq 0) then begin
;               pts_x = [ reverse(pts_fwd_x,1), pts_bwd_x[1:*] ]
;               pts_y = [ reverse(pts_fwd_y,1), pts_bwd_y[1:*] ]
;               pts_z = [ reverse(pts_fwd_z,1), pts_bwd_z[1:*] ]
;               start_points[0] = 0L
;               start_points = [start_points, fiter+biter]
;               next_start_pt = ulong(fiter+biter+1)
;               mid_points = fiter
;            endif else begin
;               if (fiter ne 0) then begin
;                   pts_x = [ pts_x, reverse(pts_fwd_x,1) ]
;                   pts_y = [ pts_y, reverse(pts_fwd_y,1) ]
;                   pts_z = [ pts_z, reverse(pts_fwd_z,1) ]
;               endif
;               if (biter ne 0) then begin
;                   pts_x = [ pts_x, pts_bwd_x[1:*] ]
;                   pts_y = [ pts_y, pts_bwd_y[1:*] ]
;                   pts_z = [ pts_z, pts_bwd_z[1:*] ]
;               endif
;               start_points = [start_points, next_start_pt, next_start_pt+fiter+biter]
;               next_start_pt += fiter+biter+1
;               mid_points = [ mid_points, fiter ]
;            endelse
;        
;         endelse
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;        pbar->update, float(seedpt_num)/n_seedpts * 100.0
;        seedpt_num++
;
;    endwhile
;
;    pbar->Destroy
;    
;    print, "Time: ", systime(1) - start_time
;
;    print, "Anis. Scaling.....: ", anis_scl
;    print, "Anis. Halt Thr....: ", anis_halt_thr
;    print, "Anis. Branch Thr..: ", anis_cross_thr
;    print, "Thru Angle Limit..: ", thru_ang_lmt
;    print, "Branch Angle Limit: ", branch_ang_lmt
;    print, "Step Size.........: ", res_fac
;    print, "Max SL Length.....: ", iter_max
;    print, "Pr Interp. Tol....: ", pr_interp_tol
;    print, "Powell Opt. Tol...: ", powell_tol
;    print, "Max Crossing Depth: ", max_level
;
;    fname = output_dir+(get_dir_slash())[0]+trs_filename
;    fname = fname + '.trs'
;    matrix = [xdim, ydim, zdim]
;    rrx = temporary(pts_x)
;    rry = temporary(pts_y)
;    rrz = temporary(pts_z)
;
;    save, filename=fname, $
;      rrx, rry, rrz, branch_levels, $
;      start_points, mid_points, $
;      matrix
;    
;    if ((size(accounting))[0] gt 3) then begin
;       accounting = total(accounting, 4)
;    endif
;
;    ;; output the voxel accounting file if requested
;    if (output_voxeldata eq 1) then begin
;       tab = string(09B)
;       sz_acct = size(accounting, /dimensions)
;       vx_fname = output_dir+(get_dir_slash())[0]+trs_filename
;       vx_fname = vx_fname + '.voxels'
;       openw, lun, vx_fname, /get_lun
;       printf, lun, 'x'+tab+'y'+tab+'z'+tab+'n'
;       fmt = '(%"%d\t%d\t%d\t%d")'
;       for zz = 0, sz_acct[2]-1 do begin
;          for xx = 0, sz_acct[0]-1 do begin
;             for yy = 0, sz_acct[1]-1 do begin
;                if (accounting[xx,yy,zz] eq 0) then continue
;                str = string(xx,yy,zz, 1, format=fmt)
;                printf, lun, str
;             endfor
;          endfor
;       endfor
;       close, lun
;       free_lun, lun
;       
;       p_accounting = ptr_new(smooth(accounting, 1))
;       (*p_accounting)[where(*p_accounting ne 0.)] = 1.0
;       mas_export_nifti, file_name=output_dir+(get_dir_slash())[0]+trs_filename+'_voxelmask.nii', $
;                         data_type='FLOAT', $
;                         data_ptr=p_accounting
;       ptr_free, p_accounting
;       
;    endif
;
;    ;; add the tract surface to the 3d view
;    if (output_surface eq 1) then begin
;       osurface_model = obj_new('idlgrmodel', depth_test_disable=2)
;       shade_volume, accounting, 1, vert, poly
;       osurf = obj_new('idlgrpolygon', vert, polygons=poly, $
;                       color=[100,100,150], alpha_channel=0.90, $
;                       shading=1, texture_interp=1)
;       osurface_model->add, osurf
;       omodel->add, osurface_model
;    endif
;    
;    vol = fltarr(size(accounting, /dimensions))
;    vol[where(accounting ne 0)] = 1.0
;    print, 'Total number of voxels included in tracts: ', total(vol)
;    vx_dim = [res_x, res_y, res_z] * 10.
;    print, 'Voxel Dimensons (mm)  : ', vx_dim
;    tmp = vx_dim / min(vx_dim)
;    print, 'Estmated Volume (mm^3): ', (total(vol)*product(tmp))*product(vx_dim)
;
;    make_fibers, fname, fiber_data=fibers, /use_new_method;;, transform=invert(paths_to_rawdata)
;    view_fibers, fibers, /show_axis, alpha=.5, thick=1, /in_mas, supplemental_omodel=omodel
;   
;    return
;    
;end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_gui_cleanup
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; 
;; Called by xmanager when window is destroyed. 
;;
;; Editing Information:
;;

pro mas_probtrack_gui_cleanup, tlb

    widget_control, tlb, get_uvalue=state
    if (ptr_valid(state)) then ptr_free, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_make_default_state
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Creates a probtrack state structure with default parameters
;;
;; Editing Information:
;;

function mas_probtrack_make_default_state

    common scan_data

    state = ptr_new({ PROBTRACK_STATE, $
                      step_size: 0.5, $
                      seed_density: 1, $
                      int_method: 0, $ ; 0 = euler, 1 = rk4
                      max_length: 5000L, $
                      min_length: 0, $
                      res_crossings: 0, $
                      show_crossings: 0, $
                      anis_halt_thr: 0.175, $
                      anis_branch_thr: 0.3, $
                      thru_ang_lmt: 45., $
                      branch_ang_lmt: 45., $
                      prob_comp_method: 0, $
                      prob_surf_tess: 3, $    ;;26, $
                      prob_interp_tol: 0.25, $ ;;0.1, $
                      prob_opt_tol: 0.001, $
                      prob_branch_tol: 0.25, $
                      max_branch_lvl: 1, $
                      max_branch_vxl: 2, $
                      roi_filelist_3D: strarr(5), $
                      roi_filelist_3d_types: bytarr(5), $
                      output_directory: project.current_path, $
                      output_filename: 'probtrack_tracts', $
                      output_stream: 1, $
                      output_surface: 0, $
                      output_voxeldata: 0 })
    
    return, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_gui2state
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Transfers parameters from the input fields on the GUI window
;; into the probtrack state structure.
;;
;; Editing Information:
;;

pro mas_probtrack_gui2state, tlb

   common scan_data

   widget_control, tlb, get_uvalue=wid_state
   
   if not ptr_valid(wid_state) then begin
      junk = dialog_message('Widget State is not consistent.', /error, /center)
      return
   endif

   state = project.procpramarray[project.ci].probtrack_state

   if not ptr_valid(state) then begin
      state = mas_probtrack_make_default_state()
      project.procpramarray[project.ci].probtrack_state = state
   endif

   widget_control, (*wid_state).txt_stepsize, get_value=step_size
   (*state).step_size = (float(strcompress(step_size, /remove_all)))[0]

   widget_control, (*wid_state).txt_density, get_value=seed_density
   (*state).seed_density = (fix(strcompress(seed_density, /remove_all)))[0]
   
   widget_control, (*wid_state).txt_maxlength, get_value=max_length
   (*state).max_length = (floor(float(strcompress(max_length, /remove_all))))[0]
   
   widget_control, (*wid_state).txt_minlength, get_value=min_length
   (*state).min_length = (floor(float(strcompress(min_length, /remove_all))))[0]
   
   (*state).res_crossings  = widget_info((*wid_state).btn_res_crossings, /button_set)
   (*state).show_crossings = widget_info((*wid_state).btn_show_crossings, /button_set)
   
   widget_control, (*wid_state).txt_anis_halt_thr, get_value=anis_halt_thr
   (*state).anis_halt_thr = (float(strcompress(anis_halt_thr, /remove_all)))[0]
   
   widget_control, (*wid_state).txt_anis_branch_thr, get_value=anis_branch_thr
   (*state).anis_branch_thr = (float(strcompress(anis_branch_thr, /remove_all)))[0]
   
   widget_control, (*wid_state).txt_thru_ang_lmt, get_value=thru_ang_lmt
   (*state).thru_ang_lmt = (float(strcompress(thru_ang_lmt, /remove_all)))[0]
   
   widget_control, (*wid_state).txt_branch_ang_lmt, get_value=branch_ang_lmt
   (*state).branch_ang_lmt = (float(strcompress(branch_ang_lmt, /remove_all)))[0]
   
   widget_control, (*wid_state).txt_prob_surf_tess, get_value=prob_surf_tess
   (*state).prob_surf_tess = (floor(float(strcompress(prob_surf_tess, /remove_all))))[0]

   prob_method = widget_info((*wid_state).dl_prob_method, /droplist_select)
   (*state).prob_comp_method = prob_method
   
   widget_control, (*wid_state).txt_prob_interp_tol, get_value=prob_interp_tol
   (*state).prob_interp_tol = (float(strcompress(prob_interp_tol, /remove_all)))[0]
   
;   widget_control, (*wid_state).txt_prob_opt_tol, get_value=prob_opt_tol
;   (*state).prob_opt_tol = (float(strcompress(prob_opt_tol, /remove_all)))[0]
   
   widget_control, (*wid_state).txt_prob_branch_thr, get_value=prob_branch_thr
   (*state).prob_branch_tol = (float(strcompress(prob_branch_thr, /remove_all)))[0]
   
   widget_control, (*wid_state).txt_max_branch_lvl, get_value=max_branch_lvl
   (*state).max_branch_lvl = (floor(float(strcompress(max_branch_lvl, /remove_all))))[0]
   
   widget_control, (*wid_state).txt_max_branch_voxel, get_value=max_branch_vxl
   (*state).max_branch_vxl = (floor(float(strcompress(max_branch_vxl, /remove_all))))[0]
   
   widget_control, (*wid_state).txt_output_directory, get_value=output_directory
   (*state).output_directory = strcompress(output_directory, /remove_all)
   
   widget_control, (*wid_state).txt_output_trs_filename, get_value=output_trs_filename
   (*state).output_filename = strcompress(output_trs_filename, /remove_all)

   output_surface = widget_info((*wid_state).btn_surface, /button_set)
   (*state).output_surface = output_surface

   output_stream = widget_info((*wid_state).btn_stream, /button_set)
   (*state).output_stream = output_stream
 
   output_voxeldata = widget_info((*wid_state).btn_voxeldata, /button_set)
   (*state).output_voxeldata = output_voxeldata
  
   for i = 0, 4 do begin
       roi_num = string(i+1, format='(I0)')

       name = 'btn_roi_enable_'+roi_num
       roi_enable_wid = widget_info(tlb, find_by_uname=name)
       if (not widget_info(roi_enable_wid, /button_set)) then begin
           (*state).roi_filelist_3d[i] = ''
           continue
       endif
       
       name = 'txt_roi_path_'+roi_num
       roi_path_wid = widget_info(tlb, find_by_uname=name)
       widget_control, roi_path_wid, get_value=val
       (*state).roi_filelist_3d[i] = val
       
       name = 'dl_roi_type_'+roi_num
       roi_type_wid = widget_info(tlb, find_by_uname=name)
       sel = widget_info(roi_type_wid, /droplist_select)
       (*state).roi_filelist_3d_types[i] = sel
       print, (*state).roi_filelist_3d_types

   endfor
end

pro mas_probtrack_gui_select_roi_event, event

    common scan_data
    
    widget_control, event.top, get_uvalue=wid_state
    
    name = widget_info(event.id, /uname)
    
    roi_info   = stregex(name, '^btn_roi_(choose|enable)_([0-9]+)$', /subexpr, /extract)
    roi_action = roi_info[1]
    roi_num    = roi_info[2]
    
    case roi_action of
    
        'enable': begin
        
            button = widget_info(event.top, find_by_uname='btn_roi_choose_'+roi_num)
            widget_control, button, sensitive=event.select
            
            textbox = widget_info(event.top, find_by_uname='txt_roi_path_'+roi_num)
            widget_control, textbox, sensitive=event.select
            
            droplist = widget_info(event.top, find_by_uname='dl_roi_type_'+roi_num)
            widget_control, droplist, sensitive=event.select
        
        end
        
        'choose': begin

            file = dialog_pickfile(path=project.current_path, /read, /file, filter='*.nii')
            if (file eq '') then return
            
            textbox = widget_info(event.top, find_by_uname='txt_roi_path_'+roi_num)
            if (widget_info(textbox, /valid_id)) then begin
                widget_control, textbox, set_value=file
            endif
                
        end
        
        else: print, "wot?"
    
    endcase
    

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_gui_event
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Event handler for probtrack GUI window
;;
;; Editing Information:
;;

pro mas_probtrack_gui_event, event

    common scan_data

    widget_control, event.top, get_uvalue=wid_state

    pt_state = project.procpramarray[project.ci].probtrack_state

    if (not ptr_valid(pt_state)) then begin
       mas_probtrack_gui2state, event.top
       pt_state = project.procpramarray[project.ci].probtrack_state
    endif

    name = widget_info(event.id, /uname)

    case name of 

        'btn_specify_roi': begin
            
            mas_probtrack_gui2state, event.top

            case (*pt_state).prob_comp_method of 

               0: begin
                  method_state = project.procpramarray[project.ci].odf_state

                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_odf_gui_make_default_state()
                     project.procpramarray[project.ci].odf_state = method_state
                  endif

                  mas_odf_init_obj, $
                     lmax=(*method_state).odf_lmax, $
                     odf_obj=method, /bulk
               end

               1: begin
                  method_state = project.procpramarray[project.ci].dot_state
                  
                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_dot_gui_make_default_state()
                     project.procpramarray[project.ci].dot_state = method_state
                  endif
                  
                  mas_dot_init_obj, $
                     lmax=(*method_state).dot_lmax, $
                     radius=(*method_state).dot_rzero, $
                     dot_obj=method, /bulk
               end

               2: begin
                  method_state = project.procpramarray[project.ci].mow_state
                  
                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_mow_gui_make_default_state()
                     project.procpramarray[project.ci].mow_state = method_state
                  endif
                  
                  mas_mow_init_obj, $
                     deco_level=(*method_state).mow_deco_order, $
                     radius=(*method_state).mow_radius, $
                     comp_meth=(*method_state).mow_comp_method, $
                     mow_obj=method, /bulk
               end

            endcase

            mas_probtrack_track, method, /use_state_ptr, /use_ft_roi   
            obj_destroy, method

         end

 
        'btn_repeat_track': begin

            mas_probtrack_gui2state, event.top

            case (*pt_state).prob_comp_method of 

                0: begin
                  method_state = project.procpramarray[project.ci].odf_state

                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_odf_gui_make_default_state()
                     project.procpramarray[project.ci].odf_state = method_state
                  endif

                  mas_odf_init_obj, $
                     lmax=(*method_state).odf_lmax, $
                     odf_obj=method, /bulk
               end

               1: begin
                  method_state = project.procpramarray[project.ci].dot_state
                  
                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_dot_gui_make_default_state()
                     project.procpramarray[project.ci].dot_state = method_state
                  endif
                  
                  mas_dot_init_obj, $
                     lmax=(*method_state).dot_lmax, $
                     radius=(*method_state).dot_rzero, $
                     dot_obj=method, /bulk
               end

               2: begin
                  method_state = project.procpramarray[project.ci].mow_state
                  
                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_mow_gui_make_default_state()
                     project.procpramarray[project.ci].mow_state = method_state
                  endif
                  
                  mas_mow_init_obj, $
                     deco_level=(*method_state).mow_deco_order, $
                     radius=(*method_state).mow_radius, $
                     comp_meth=(*method_state).mow_comp_method, $
                     mow_obj=method, /bulk
               end

            endcase

            mas_probtrack_track, method, /use_state_ptr, /use_roi_file
            obj_destroy, method

        end

        'btn_select_dir': begin

            widget_control, (*wid_state).txt_output_directory, get_value=curr_path
            print, curr_path
            dirname = dialog_pickfile(path=curr_path, /directory)

            if (dirname eq '') then return

            widget_control, (*wid_state).txt_output_directory, set_value=dirname            

        end

        'btn_revisualize': begin
            
            widget_control, (*wid_state).txt_output_directory, get_value=output_dir
            widget_control, (*wid_state).txt_output_trs_filename, get_value=output_file

            junk = strmid(output_dir, strlen(output_dir)-1)
            if (junk eq (get_dir_slash())[0]) then begin
                output_dir = strmid(output_dir, 0, strlen(output_dir)-1)
            endif

            tract_file = output_dir + (get_dir_slash())[0] + output_file + '.trs'

            if not file_test(tract_file, /regular, /read) then begin
                junk = dialog_message('Tract file: '+tract_file+ $
                                      ' not found or cannot be read.', /error, /center)
                return
            endif
            
            make_fibers, tract_file, fiber_data=fiber_data, /use_new_method
            
            datasize = size(*project.dataarray[project.ci].state1)
            
            result = total([ datasize[1], datasize[2], datasize[3] ] - (*fiber_data).matrix)
            
            if (result ne 0) then begin
                junk = dialog_message('Project dimensions do not match '+$
                                      'tract file. Aborting', /error, /center)
                ptr_free, (*fiber_data).fibers
                ptr_free, fiber_data
                return
            endif
            
            view_fibers, fiber_data, thick=1, alpha=.75, /show_axis, /in_mas
            
        end
        
        'btn_config_prob': begin
            base = widget_base(title="Configure Probability", group_leader=event.top, /modal)

            mas_probtrack_gui2state, event.top

            case (*pt_state).prob_comp_method of 
               0: begin
                  mas_odf_gui_make, base, 0
                  widget_control, base, /realize
                  xmanager, 'mas_odf_gui', base, cleanup='mas_odf_gui_cleanup'
               end
               1: begin
                  mas_dot_gui_make, base, 0
                  widget_control, base, /realize
                  xmanager, 'mas_dot_gui', base, cleanup='mas_dot_gui_cleanup'
               end
               2: begin
                  mas_mow_gui_make, base, 0
                  widget_control, base, /realize
                  xmanager, 'mas_mow_gui', base, cleanup='mas_mow_gui_cleanup'
               end

            endcase
        end

        'btn_res_crossings': begin
            
            if (event.select eq 0) then begin
                
                widget_control, (*wid_state).txt_anis_branch_thr, sensitive = 0
                widget_control, (*wid_state).txt_branch_ang_lmt, sensitive = 0
                widget_control, (*wid_state).txt_prob_interp_tol, sensitive = 0
                ;widget_control, (*wid_state).txt_prob_opt_tol, sensitive = 0
                widget_control, (*wid_state).txt_prob_branch_thr, sensitive = 0
                widget_control, (*wid_state).txt_max_branch_lvl, sensitive = 0
                widget_control, (*wid_state).txt_max_branch_voxel, sensitive = 0
                widget_control, (*wid_state).btn_show_crossings, sensitive = 0
                
            endif else begin

                widget_control, (*wid_state).txt_anis_branch_thr, sensitive = 1
                widget_control, (*wid_state).txt_branch_ang_lmt, sensitive = 1
                widget_control, (*wid_state).txt_prob_interp_tol, sensitive = 1
                ;widget_control, (*wid_state).txt_prob_opt_tol, sensitive = 1
                widget_control, (*wid_state).txt_prob_branch_thr, sensitive = 1
                widget_control, (*wid_state).txt_max_branch_lvl, sensitive = 1
                widget_control, (*wid_state).txt_max_branch_voxel, sensitive = 1
                widget_control, (*wid_state).btn_show_crossings, sensitive = 1

            endelse

        end

        'btn_examine_vxl': begin
            
            data_sz = size(*project.dataarray[project.ci].state1, /dimensions)

            mas_probtrack_gui2state, event.top
            widget_control, (*wid_state).txt_show_vxl_x, get_value=vxl_x
            if (vxl_x eq '') then begin
               vxl_x = project.procpramarray[project.ci].fdim_start
            endif else begin
               vxl_x = (float(strcompress(vxl_x, /remove_all)))[0]
            endelse

            widget_control, (*wid_state).txt_show_vxl_y, get_value=vxl_y
            if (vxl_y eq '') then begin
               vxl_y = project.procpramarray[project.ci].pdim_start
            endif else begin
               vxl_y = (float(strcompress(vxl_y, /remove_all)))[0]
            endelse

            widget_control, (*wid_state).txt_show_vxl_z, get_value=vxl_z
            if (vxl_z eq '') then begin
               vxl_z = project.procpramarray[project.ci].sdim_start
            endif else begin
               vxl_z = (float(strcompress(vxl_z, /remove_all)))[0]
            endelse

            vxl_x = ( vxl_x > 0 ) < (data_sz[0]-1)
            vxl_y = ( vxl_y > 0 ) < (data_sz[1]-1)
            vxl_z = ( vxl_z > 0 ) < (data_sz[2]-1)

            case (*pt_state).prob_comp_method of 

               0: begin
                  method_state = project.procpramarray[project.ci].odf_state

                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_odf_gui_make_default_state()
                     project.procpramarray[project.ci].odf_state = method_state
                  endif

                  mas_odf_init_obj, $
                     lmax=(*method_state).odf_lmax, $
                     odf_obj=method, /bulk
               end

               1: begin
                  method_state = project.procpramarray[project.ci].dot_state
                  
                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_dot_gui_make_default_state()
                     project.procpramarray[project.ci].dot_state = method_state
                  endif
                  
                  mas_dot_init_obj, $
                     lmax=(*method_state).dot_lmax, $
                     radius=(*method_state).dot_rzero, $
                     dot_obj=method, /bulk
               end

               2: begin
                  method_state = project.procpramarray[project.ci].mow_state
                  
                  if (not ptr_valid(method_state)) then begin
                     method_state = mas_mow_gui_make_default_state()
                     project.procpramarray[project.ci].mow_state = method_state
                  endif
                  
                  mas_mow_init_obj, $
                     deco_level=(*method_state).mow_deco_order, $
                     radius=(*method_state).mow_radius, $
                     comp_meth=(*method_state).mow_comp_method, $
                     mow_obj=method, /bulk
               end

            endcase

            mas_probtrack_show_voxel, method, [vxl_x, vxl_y, vxl_z], /use_state_ptr
            
            obj_destroy, method ;ptr_free, method

        end

        else: begin
            print, "mas_probtrack_gui_event: Unknown event from: "+name
        end
        
    endcase

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack_gui
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;      Creates the probtrack GUI window and input fields
;;
;; Editing Information:
;;

pro mas_probtrack_gui, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

  common scan_data
  common common_widgets

  ci = project.ci

  if ptr_valid(project.procpramarray[ci].probtrack_state) then begin
      state = project.procpramarray[ci].probtrack_state
  endif else begin
      state = mas_probtrack_make_default_state()
      project.procpramarray[ci].probtrack_state = state
  endelse

  main = Widget_Base( GROUP_LEADER=WID_BASE_MAIN, UNAME='main'  $
      ,XOFFSET=5 ,YOFFSET=5 ,TITLE='ProbTrack - Probability-based Tracking' $
      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  tab_base = widget_tab(main)
  
  param_tab = widget_base(tab_base, row=1, title="Tracking Parameters")
  roi_tab   = widget_base(tab_base, title="ROI Selection")
  
  WID_BASE_0 = Widget_Base(param_tab, UNAME='WID_BASE_0' ,XOFFSET=3  $
      ,YOFFSET=3 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,COLUMN=1)

  ;;;;;;;;;;;;;;;;;;;[ streamline parameters ];;;;;;;;;;;;;;;;;;;;;;;;

  base_streamline = Widget_Base(WID_BASE_0, UNAME='base_streamline'  $
      ,FRAME=1 ,XOFFSET=3 ,YOFFSET=3 ,TITLE='Streamlining' ,SPACE=3  $
      ,XPAD=3 ,YPAD=3 ,COLUMN=1)

  
  base_int_method = Widget_Base(base_streamline,  $
      UNAME='base_int_method' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_int_method = Widget_Label(base_int_method,  $
      UNAME='lbl_int_method' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Integration Method:')

  
  dl_int_method = Widget_Droplist(base_int_method,  $
      UNAME='dl_int_method' ,XOFFSET=100 ,YOFFSET=3 ,VALUE=[ 'Eulers'+ $
      ' Method' ])

  base_stepsize = Widget_Base(base_streamline, UNAME='base_stepsize'  $
      ,XOFFSET=3 ,YOFFSET=34 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,ROW=1)

  
  lbl_stepsize = Widget_Label(base_stepsize, UNAME='lbl_stepsize'  $
      ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Step Size:')

  
  txt_stepsize = Widget_Text(base_stepsize, UNAME='txt_stepsize'  $
      ,XOFFSET=56 ,YOFFSET=3 ,/EDITABLE $
      ,VALUE=string((*state).step_size, format='(1F0.2)') $
      ,XSIZE=6 ,YSIZE=1)

;;;;
  lbl_density = Widget_Label(base_stepsize, UNAME='lbl_density' $
      ,VALUE="Seed Density:", /align_left)
  txt_density = Widget_Text(base_stepsize, UNAME='txt_density'  $
      ,XOFFSET=56 ,YOFFSET=3 ,/EDITABLE $
      ,VALUE=string((*state).seed_density, format='(2I0)') $
      ,XSIZE=4 ,YSIZE=1)
;;;;
  
  base_maxlength = Widget_Base(base_streamline, UNAME='base_maxlength'  $
      ,XOFFSET=3 ,YOFFSET=64 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,ROW=1)

  
  lbl_maxlength = Widget_Label(base_maxlength,  $
      UNAME='lbl_maxlength' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Max. Stream Length:')

  
  txt_maxlength = Widget_Text(base_maxlength,  $
      UNAME='txt_maxlength' ,XOFFSET=127 ,YOFFSET=3 ,/EDITABLE  $
      ,VALUE=string((*state).max_length, format='(8I0)') $
      ,XSIZE=8 ,YSIZE=1)

  
  base_minlength = Widget_Base(base_streamline, UNAME='base_minlength'  $
      ,XOFFSET=3 ,YOFFSET=94 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,ROW=1)

  
  lbl_minlength = Widget_Label(base_minlength, UNAME='lbl_minlength'  $
      ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Min. Stream'+ $
      ' Length:')

  
  txt_minlength = Widget_Text(base_minlength, UNAME='txt_minlength'  $
      ,XOFFSET=124 ,YOFFSET=3 ,/EDITABLE ,XSIZE=8  $
      ,VALUE=string((*state).min_length, format='(8I0)') $
      ,YSIZE=1)

  base_res_crossings = Widget_Base(base_streamline,  $
      UNAME='base_res_crossings' ,XOFFSET=3 ,YOFFSET=124 ,TITLE='IDL'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  WID_BASE_1 = Widget_Base(base_res_crossings, UNAME='WID_BASE_1'  $
      ,XOFFSET=3 ,YOFFSET=3 ,TITLE='IDL' ,ROW=1 ,/NONEXCLUSIVE)

  
  btn_res_crossings = Widget_Button(WID_BASE_1,  $
      UNAME='btn_res_crossings' ,/ALIGN_LEFT ,VALUE='Resolve'+ $
      ' Crossings')

  btn_show_crossings = Widget_Button(WID_BASE_1,  $
      UNAME='btn_show_crossings' ,/ALIGN_LEFT ,VALUE='Show'+ $
      ' Crossings')

  widget_control, btn_res_crossings, set_button=(*state).res_crossings
  widget_control, btn_show_crossings, set_button=(*state).show_crossings

  ;;;;;;;;;;;;;;;;;;;[ Anis. and Angle parameters ];;;;;;;;;;;;;;;;;;;;;;;;

  
  WID_BASE_6 = Widget_Base(WID_BASE_0, UNAME='WID_BASE_6' ,FRAME=1  $
      ,XOFFSET=3 ,YOFFSET=161 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,COLUMN=1)

  base_int_method_2 = Widget_Base(WID_BASE_6,  $
      UNAME='base_int_method_2' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_anis_measure = Widget_Label(base_int_method_2,  $
      UNAME='lbl_anis_measure' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Anisotropy Measure:')

  
  dl_anis_measure = Widget_Droplist(base_int_method_2,  $
      UNAME='dl_anis_measure' ,XOFFSET=104 ,YOFFSET=3 ,VALUE=[  $
      'Fractional Anis.'])

  
  base_anis_halt_thr = Widget_Base(WID_BASE_6,  $
      UNAME='base_anis_halt_thr' ,XOFFSET=3 ,YOFFSET=34 ,TITLE='IDL'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_anis_halt_thr = Widget_Label(base_anis_halt_thr,  $
      UNAME='lbl_anis_halt_thr' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Halting Anisotropy Threshold:')

  
  txt_anis_halt_thr = Widget_Text(base_anis_halt_thr,  $
      UNAME='txt_anis_halt_thr' ,XOFFSET=146 ,YOFFSET=3 ,/EDITABLE  $
      ,VALUE=string((*state).anis_halt_thr, format='(1F0.3)') ,XSIZE=8 ,YSIZE=1)

  
  base_anis_branch_thr = Widget_Base(WID_BASE_6,  $
      UNAME='base_anis_branch_thr' ,XOFFSET=3 ,YOFFSET=64  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_anis_branch_thr = Widget_Label(base_anis_branch_thr,  $
      UNAME='lbl_anis_branch_thr' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Branching Anisotropy Threshold:')

  
  txt_anis_branch_thr = Widget_Text(base_anis_branch_thr,  $
      UNAME='txt_anis_branch_thr' ,XOFFSET=161 ,YOFFSET=3 ,/EDITABLE  $
      ,SENSITIVE=(*state).res_crossings $
      ,VALUE=string((*state).anis_branch_thr, format='(1F0.3)') ,XSIZE=8 ,YSIZE=1)

  
  base_thru_ang_lmt = Widget_Base(WID_BASE_6,  $
      UNAME='base_thru_ang_lmt' ,XOFFSET=3 ,YOFFSET=94 ,TITLE='IDL'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_thru_angle_lmt = Widget_Label(base_thru_ang_lmt,  $
      UNAME='lbl_thru_angle_lmt' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Thru Angle Limit:')

  txt_thru_ang_lmt = Widget_Text(base_thru_ang_lmt,  $
      UNAME='txt_thru_ang_lmt' ,XOFFSET=87 ,YOFFSET=3 ,/EDITABLE  $
      ,VALUE=string((*state).thru_ang_lmt, format='(2F0.2)') ,XSIZE=5 ,YSIZE=1)

  
  WID_LABEL_0 = Widget_Label(base_thru_ang_lmt, UNAME='WID_LABEL_0'  $
      ,XOFFSET=142 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='(Deg.)')

  
  base_branch_ang_lmt = Widget_Base(WID_BASE_6,  $
      UNAME='base_branch_ang_lmt' ,XOFFSET=3 ,YOFFSET=124  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_branch_angle_lmt = Widget_Label(base_branch_ang_lmt,  $
      UNAME='lbl_branch_angle_lmt' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Branch Angle Limit:')

  
  txt_branch_ang_lmt = Widget_Text(base_branch_ang_lmt,  $
      UNAME='txt_branch_ang_lmt' ,XOFFSET=99 ,YOFFSET=3 ,/EDITABLE  $
      ,SENSITIVE=(*state).res_crossings $
      ,VALUE=string((*state).branch_ang_lmt, format='(2F0.2)') ,XSIZE=5 ,YSIZE=1)

  
  WID_LABEL_1 = Widget_Label(base_branch_ang_lmt, UNAME='WID_LABEL_1'  $
      ,XOFFSET=154 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='(Deg.)')

  
  ;;;;;;;;;;;;;;;;;;;[ Probability Thresholds  ];;;;;;;;;;;;;;;;;;;;;;;;

  WID_BASE_4 = Widget_Base(param_tab, UNAME='WID_BASE_4'  $
      ,XOFFSET=234 ,YOFFSET=3 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,COLUMN=1)

  
  WID_BASE_3 = Widget_Base(WID_BASE_4, UNAME='WID_BASE_3' ,FRAME=1  $
      ,XOFFSET=3 ,YOFFSET=3 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,COLUMN=1)

  base_prob_method = Widget_Base(WID_BASE_3, $
      UNAME='base_prob_method' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_prob_method = Widget_Label(base_prob_method,  $
      UNAME='lbl_prob_method' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Probability Computation:')

  dl_prob_method = Widget_Droplist(base_prob_method,  $
      UNAME='dl_prob_method' ,XOFFSET=100 ,YOFFSET=3 ,VALUE=[ 'ODF', 'DOT', 'MOW' ])

  btn_config_prob = widget_button(base_prob_method, $
      UNAME='btn_config_prob', VALUE='Configure', SCR_XSIZE=80, SCR_YSIZE=10)

  base_prob_surf_tess = Widget_Base(WID_BASE_3,  $
      UNAME='base_prob_surf_tess' ,XOFFSET=3 ,YOFFSET=3 ,TITLE='IDL'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_prob_surf_tess = Widget_Label(base_prob_surf_tess,  $
      UNAME='lbl_prob_surf_tess' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Probability Surface Tesselation:')

  
  txt_prob_surf_tess = Widget_Text(base_prob_surf_tess,  $
      UNAME='txt_prob_surf_tess' ,XOFFSET=138 ,YOFFSET=3 ,/EDITABLE  $
      ,VALUE=string((*state).prob_surf_tess, format='(3I0)') ,XSIZE=4 ,YSIZE=1)

  
  base_prob_interp_tol = Widget_Base(WID_BASE_3,  $
      UNAME='base_prob_interp_tol' ,XOFFSET=3 ,YOFFSET=3 ,TITLE='IDL'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_prob_interp_tol = Widget_Label(base_prob_interp_tol,  $
      UNAME='lbl_prob_interp_tol' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Probability Interpolation Tol:')

  
  txt_prob_interp_tol = Widget_Text(base_prob_interp_tol,  $
      UNAME='txt_prob_interp_tol' ,XOFFSET=138 ,YOFFSET=3 ,/EDITABLE  $
      ,SENSITIVE=(*state).res_crossings $
      ,VALUE=string((*state).prob_interp_tol, format='(1F0.3)') ,XSIZE=8 ,YSIZE=1)

  
;  base_prob_opt_tol = Widget_Base(WID_BASE_3,  $
;      UNAME='base_prob_opt_tol' ,XOFFSET=3 ,YOFFSET=33 ,TITLE='IDL'  $
;      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)
;
;  
;  lbl_prob_opt_tol = Widget_Label(base_prob_opt_tol,  $
;      UNAME='lbl_prob_opt_tol' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
;      ,VALUE='Probability Optimization Tol:')
;
;  
;  txt_prob_opt_tol = Widget_Text(base_prob_opt_tol,  $
;      UNAME='txt_prob_opt_tol' ,XOFFSET=134 ,YOFFSET=3  $
;      ,SENSITIVE=(*state).res_crossings $
;      ,/EDITABLE ,VALUE=string((*state).prob_opt_tol, format='(1F0.5)') ,XSIZE=8 ,YSIZE=1)
;
  
  base_prob_branch_thr = Widget_Base(WID_BASE_3,  $
      UNAME='base_prob_branch_thr' ,XOFFSET=3 ,YOFFSET=63  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_prob_branch_thr = Widget_Label(base_prob_branch_thr,  $
      UNAME='lbl_prob_branch_thr' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Probability Threshold:')

  
  txt_prob_branch_thr = Widget_Text(base_prob_branch_thr,  $
      UNAME='txt_prob_branch_thr' ,XOFFSET=146 ,YOFFSET=3 ,/EDITABLE  $
      ,SENSITIVE=(*state).res_crossings $
      ,VALUE=string((*state).prob_branch_tol, format='(1F0.4)') ,XSIZE=8 ,YSIZE=1)
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;; [ Branch Levels ];;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  WID_BASE_5 = Widget_Base(WID_BASE_4, UNAME='WID_BASE_5' ,FRAME=1  $
      ,XOFFSET=3 ,YOFFSET=99 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,ROW=1)

  
  base_maxl_branch_lvl = Widget_Base(WID_BASE_5,  $
      UNAME='base_maxl_branch_lvl' ,XOFFSET=3 ,YOFFSET=3 ,TITLE='IDL'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_maxlength_1 = Widget_Label(base_maxl_branch_lvl,  $
      UNAME='lbl_maxlength_1' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Max Branch Level:')

  
  txt_max_branch_lvl = Widget_Text(base_maxl_branch_lvl,  $
      UNAME='txt_max_branch_lvl' ,XOFFSET=121 ,YOFFSET=3 ,/EDITABLE  $
      ,SENSITIVE=(*state).res_crossings $
      ,VALUE=string((*state).max_branch_lvl, format='(3I0)') ,XSIZE=4 ,YSIZE=1)

  
  base_max_branch_voxel = Widget_Base(WID_BASE_5,  $
      UNAME='base_max_branch_voxel' ,XOFFSET=3 ,YOFFSET=33  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,ROW=1)

  
  lbl_max_branch_voxel = Widget_Label(base_max_branch_voxel,  $
      UNAME='lbl_max_branch_voxel' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Max Branches/Voxel:')

  
  txt_max_branch_voxel = Widget_Text(base_max_branch_voxel,  $
      UNAME='txt_max_branch_voxel' ,XOFFSET=134 ,YOFFSET=3 ,/EDITABLE  $
      ,SENSITIVE=(*state).res_crossings $
      ,VALUE=string((*state).max_branch_vxl, format='(3I0)') ,XSIZE=4 ,YSIZE=1)

;;;;;;;;;;;;;;;;;;;;; [ Directories ] ;;;;;;;;;;;;;;;;;;;;;;;;;;

  WID_BASE_7 = Widget_Base(WID_BASE_4, UNAME='WID_BASE_7' ,FRAME=1  $
      ,XOFFSET=3 ,YOFFSET=165 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,COLUMN=1)

  
  lbl_output_directory = Widget_Label(WID_BASE_7,  $
      UNAME='lbl_output_directory' ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT  $
      ,VALUE='Output Directory:')

  
  txt_output_directory = Widget_Text(WID_BASE_7,  $
      UNAME='txt_output_directory' ,XOFFSET=3 ,YOFFSET=21 ,/EDITABLE  $
      ,XSIZE=30 ,YSIZE=1, value=(*state).output_directory)

  junk = widget_base(WID_BASE_7, /row)

  lbl_output_trackfile = Widget_Label(junk, $;WID_BASE_7,  $
      UNAME='lbl_output_trackfile' ,XOFFSET=3 ,YOFFSET=45  $
      ,/ALIGN_LEFT ,VALUE='Output File Name:')

  
  txt_output_trs_filename = Widget_Text(junk, $;WID_BASE_7,  $
      UNAME='txt_output_trs_filename' ,XOFFSET=3 ,YOFFSET=63  $
      ,/EDITABLE ,XSIZE=20 ,YSIZE=1, value=(*state).output_filename)

  lbl_ext = widget_label(junk, value=".trs,.ROI,etc.")

  junk = widget_base(WID_BASE_7, /row)
  
  lbl_type = widget_label(junk, value="Output Types:")
  
  junk = widget_base(junk, /row, /nonexclusive)
  
  btn_surface = widget_button(junk, value="Surface", uname="btn_surface")
  widget_control, btn_surface, set_button=(*state).output_surface
  btn_stream  = widget_button(junk, value="Streamlines", uname="btn_stream")
  widget_control, btn_stream, set_button=(*state).output_stream
  btn_voxeldata  = widget_button(junk, value="VoxelData", uname="btn_voxeldata")
  widget_control, btn_voxeldata, set_button=(*state).output_voxeldata

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;[ ROI Tab ];;;;;;;;;;;;;;;;;;;;;;;

  b = widget_base(roi_tab, /column)
  for roi = 1, 5 do begin

      base_roi_1 = widget_base(b, /row)
      i = string(roi, format='(I0)')

      non_ex_base = widget_base(base_roi_1, /nonexclusive)
      btn = widget_button(non_ex_base, value='3D ROI #'+i+': ', $
                            uname='btn_roi_enable_'+i, sensitive=1, $
                            event_pro='mas_probtrack_gui_select_roi_event')
      
      txt = widget_text(base_roi_1, xsize=60, uname='txt_roi_path_'+i)

      btn = widget_button(base_roi_1, sensitive=0, $
                            uname='btn_roi_choose_'+i, value='Choose...',$
                            event_pro='mas_probtrack_gui_select_roi_event')

      dl = widget_droplist(base_roi_1, value=['Seeding', 'Trapping', 'Terminator'], uname='dl_roi_type_'+i, $
                             event_pro='mas_probtrack_gui_select_roi_event', sensitive=0)

  endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;; [ Buttons ] ;;;;;;;;;;;;;;;;;;;;;;
  
  WID_BASE_8 = Widget_Base(main, UNAME='WID_BASE_8'  $
      ,XOFFSET=450 ,YOFFSET=3 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,COLUMN=1)

  
  btn_specify_roi = Widget_Button(WID_BASE_8, UNAME='btn_specify_roi'  $
      ,XOFFSET=12 ,YOFFSET=3 ,SCR_XSIZE=115 ,SCR_YSIZE=32  $
      ,/ALIGN_CENTER ,VALUE='Begin Tracking')

 
  btn_repeat_track = Widget_Button(WID_BASE_8,  $
      UNAME='btn_repeat_track' ,XOFFSET=7 ,YOFFSET=28 ,SCR_XSIZE=115  $
      ,SCR_YSIZE=32 ,/ALIGN_CENTER ,VALUE='Repeat Track')

  btn_select_dir = Widget_Button(WID_BASE_8, UNAME='btn_select_dir'  $
      ,XOFFSET=3 ,YOFFSET=53 ,SCR_XSIZE=115 ,SCR_YSIZE=32  $
      ,/ALIGN_CENTER ,VALUE='Select Directory')

  
  btn_revisualize = Widget_Button(WID_BASE_8, UNAME='btn_revisualize'  $
      ,XOFFSET=13 ,YOFFSET=78 ,SCR_XSIZE=115 ,SCR_YSIZE=32  $
      ,/ALIGN_CENTER ,VALUE='Revisualize')

  base_show_vxl = Widget_Base(WID_BASE_8,  $
      UNAME='base_show_vxl' ,XOFFSET=3 ,YOFFSET=33, /align_center  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,/column, /frame)

  lbl_show_vxl = widget_label(base_show_vxl, value="Examine Voxel", /align_center)
  
  jb = widget_base(base_show_vxl, /row, /align_center)
  junk = widget_label(jb, value="X:")
  txt_show_vxl_x = Widget_Text(jb,  $
      UNAME='txt_show_vxl_x' ,XOFFSET=134 ,YOFFSET=3 ,/EDITABLE  $
      ,VALUE='' ,XSIZE=4 ,YSIZE=1)

  jb = widget_base(base_show_vxl, /row, /align_center)
  junk = widget_label(jb, value="Y:")
  txt_show_vxl_y = Widget_Text(jb,  $
      UNAME='txt_show_vxl_y' ,XOFFSET=134 ,YOFFSET=3 ,/EDITABLE  $
      ,VALUE='' ,XSIZE=4 ,YSIZE=1)

  jb = widget_base(base_show_vxl, /row, /align_center)
  junk = widget_label(jb, value="Z:")
  txt_show_vxl_z = Widget_Text(jb,  $
      UNAME='txt_show_vxl_z' ,XOFFSET=134 ,YOFFSET=3 ,/EDITABLE  $
      ,VALUE='' ,XSIZE=4 ,YSIZE=1)

  btn_examine_vxl = widget_button(base_show_vxl, value='Examine', $
                                  /align_center, uname='btn_examine_vxl')

  wid_state = ptr_new({$
                        dl_int_method: dl_int_method   , $
                        txt_stepsize: txt_stepsize     , $
                        txt_density: txt_density       , $
                        txt_maxlength: txt_maxlength   , $
                        txt_minlength: txt_minlength   , $
                        btn_res_crossings: btn_res_crossings , $
                        btn_show_crossings: btn_show_crossings , $
                        dl_anis_measure: dl_anis_measure, $
                        txt_anis_halt_thr: txt_anis_halt_thr, $
                        txt_anis_branch_thr: txt_anis_branch_thr, $
                        txt_thru_ang_lmt: txt_thru_ang_lmt, $
                        txt_branch_ang_lmt: txt_branch_ang_lmt, $
                        dl_prob_method: dl_prob_method, $
                        btn_config_prob: btn_config_prob, $
                        txt_prob_surf_tess: txt_prob_surf_tess, $
                        txt_prob_interp_tol: txt_prob_interp_tol, $
                        $;txt_prob_opt_tol: txt_prob_opt_tol, $
                        txt_prob_branch_thr: txt_prob_branch_thr, $
                        txt_show_vxl_x: txt_show_vxl_x, $
                        txt_show_vxl_y: txt_show_vxl_y, $
                        txt_show_vxl_z: txt_show_vxl_z, $
                        txt_max_branch_lvl: txt_max_branch_lvl, $
                        txt_max_branch_voxel: txt_max_branch_voxel, $
                        txt_output_directory: txt_output_directory , $
                        txt_output_trs_filename: txt_output_trs_filename, $
                        $ ;txt_output_roi_filename: txt_output_roi_filename, $
                        btn_specify_roi: btn_specify_roi, $
                        btn_repeat_track: btn_repeat_track, $
                        btn_select_dir: btn_select_dir, $
                        btn_revisualize: btn_revisualize, $
                        btn_stream: btn_stream, $
                        btn_surface: btn_surface, $
                        btn_voxeldata: btn_voxeldata $
                      }, /no_copy)

  widget_control, main, set_uvalue=wid_state
  widget_control, main, /realize

  xmanager, 'mas_probtrack_gui', main, cleanup='mas_probtrack_gui_cleanup', /NO_BLOCK  

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_probtrack
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; dummy procedure to call mas_probtrack_GUI
;;
;; Editing Information:
;;

pro mas_probtrack

    mas_probtrack_GUI

end
