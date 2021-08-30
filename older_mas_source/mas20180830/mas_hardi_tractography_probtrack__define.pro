;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                       ;;
;;    HIGHLY EXPERIMENTAL. Not in final form. No documentation exists.   ;;
;;                                                                       ;;
;;        Routines to perform probabilistic fiber tracking.              ;;
;;                                                                       ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: 
;; Created by: BT
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; 
;;
;; Editing Information:

function mas_hardi_tractography_probtrack::init, ndirections, theta_ptrs, phi_ptrs

    common scan_data, project
    
    ci = project.ci
    
    self.step_size = 0.25
    self.anis_halt_thr = 0.1
    self.anis_branch_thr = 0.3
    self.turn_angle_thr = cos(70.0 * !DTOR)
    self.interpolation_method = 0
    
    interp_factors = [ project.procpramarray[ci].freq_interp, $
                       project.procpramarray[ci].phase_interp, $
                       project.procpramarray[ci].slice_interp ] 

    voxel_dims = [ project.imndarray[ci].f_voxsz, $
                   project.imndarray[ci].p_voxsz, $
                   project.imndarray[ci].s_voxsz ]

    self.resolution_xyz = voxel_dims/interp_factors    
    self.resolution_xyz = self.resolution_xyz/min(self.resolution_xyz) * self.step_size

    self.theta_1_samp = theta_ptrs[0]
    self.theta_2_samp = theta_ptrs[1]
    self.phi_1_samp = phi_ptrs[0]
    self.phi_2_samp = phi_ptrs[1]
    self.diff_data_ptr = project.dataarray[project.ci].state1
    self.anis_data_ptr = project.dataarray[project.ci].frac_ani
    
    self.dist_size = (size(*theta_ptrs[0], /dimensions))[3]
    for d = 1, n_elements(theta_ptrs)-1 do begin
        tmp = (size(*theta_ptrs[d], /dimensions))[3]
        if (tmp lt self.dist_size) then self.dist_size = tmp
    endfor
    for d = 0, n_elements(phi_ptrs)-1 do begin
        tmp = (size(*phi_ptrs[d], /dimensions))[3]    
        if (tmp lt self.dist_size) then self.dist_size = tmp
    endfor
    
    return, 1
    
end

function mas_hardi_tractography_probtrack::initializeStreamlines, seedpt

    seed = self.random_seed
    streams = objarr(5000)
    nsamp = (size(*self.theta_1_samp, /dimensions))[3]-1
    
    for n = 0, 5000-1 do begin
    
        samp_ind = round(randomu(seed)*nsamp)
        if (randomu(seed) gt 0.5) then begin
           theta = !DPI/2.0-(*self.theta_1_samp)[seedpt[0], seedpt[1], seedpt[2], samp_ind]
           phi   = (*self.phi_1_samp)[seedpt[0], seedpt[1], seedpt[2], samp_ind]
        endif else begin
           theta = !DPI/2.0-(*self.theta_2_samp)[seedpt[0], seedpt[1], seedpt[2], samp_ind]
           phi   = (*self.phi_2_samp)[seedpt[0], seedpt[1], seedpt[2], samp_ind]
        endelse

        data = cv_coord(from_sphere=[ phi, theta, 1.0 ], /to_rect)
        
        dp = total([0.0,0.0,1.0] * data)
        
        max_dir = (dp lt 0) ? -data : data

        streams[n] = obj_new('mas_hardi_streamline', seedpt, max_dir)
    
    endfor
    
    self.random_seed = seed[0]
    
    return, streams
    
end

pro mas_hardi_tractography_probtrack::setExclusionMask, pmask

    if (ptr_valid(pmask)) then begin
        
        if (ptr_valid(self.excl_masks)) then begin
            ptr_free, self.excl_masks
            (*self.excl_masks) = ptr_new(*pmask)
        endif else begin
            self.excl_masks = ptr_new(*pmask)
        endelse

    endif
    
end

pro mas_hardi_tractography_probtrack::setTerminationMask, pmask

    if (ptr_valid(pmask)) then begin
        
        if (ptr_valid(self.term_masks)) then begin
            ptr_free, self.term_masks
            (*self.term_masks) = ptr_new(*pmask)
        endif else begin
            self.term_masks = ptr_new(*pmask)
        endelse

    endif
    
end

pro mas_hardi_tractography_probtrack::addWaypointMask, pmask

    if (ptr_valid(pmask)) then begin
        
        if (ptr_valid(self.wayp_masks)) then begin
            (*self.wayp_masks) = [ (*self.wayp_masks), ptr_new(*pmask) ]
        endif else begin
            self.wayp_masks = ptr_new([ ptr_new(*pmask) ])
        endelse

        ptr_free, self.wayp_check
        self.wayp_check = ptr_new(bytarr(n_elements(*self.wayp_masks)))

    endif
    
end

;function mas_hardi_tractography_probtrack::getData, coord
;    
;    random_seed = self.random_seed
;    samp_ind = round(randomu(random_seed, 2)*(200-1))
;    self.random_seed = random_seed[0]
;
;    case self.interpolation_method of
;    
;        2: begin ;; analogous to trilinear
;            data = mas_4d_data_average(coord, data_ptr=self.diff_data_ptr)
;        end
;        
;        else: begin
;            no_int = (round(coord) < (size(*self.diff_data_ptr, /dimensions)-1)) > [0,0,0]
;            theta = (*self.theta_1_samp)[no_int[0], no_int[1], no_int[2], samp_ind[0]]
;            phi   = (*self.phi_1_samp)[no_int[0], no_int[1], no_int[2], samp_ind[0]]
;        end
;        
;     endcase
;            
;end

function mas_hardi_tractography_probtrack::jump, stream

    curr_pt = stream->getCurrentPoint(direction=curr_dir)
    
    if (self->checkBoundsData(curr_pt)) then begin
        return, stream->switchDirection()
    endif
    
    no_int = (round(curr_pt) < (size(*self.diff_data_ptr, /dimensions)-1)) > [0,0,0]
    
    random_seed = self.random_seed
    samp_ind = round(randomu(random_seed, 2)*(self.dist_size-1))
    self.random_seed = random_seed[0]
    
    if (1) then begin
        ;; sample theta
        theta = [ !DPI/2.0-(*self.theta_1_samp)[no_int[0], no_int[1], no_int[2], samp_ind[0]], $
                  !DPI/2.0-(*self.theta_2_samp)[no_int[0], no_int[1], no_int[2], samp_ind[1]] ]
        ;; sample_phi
        phi = [ (*self.phi_1_samp)[no_int[0], no_int[1], no_int[2], samp_ind[0]], $
                (*self.phi_2_samp)[no_int[0], no_int[1], no_int[2], samp_ind[1]] ]
            
        data = cv_coord(from_sphere=[ transpose([ phi ]), transpose([ theta ]), transpose([replicate(1.0, 2)]) ], /to_rect)
        
        dps = fltarr(2)
        for d = 0, 1 do begin
            dps[d] = total(curr_dir * data[*,d]) / sqrt(total(curr_dir * curr_dir))
        endfor

        ord = reverse(sort(abs(dps)))
        rnd_dir_ind = 0;(randomu(random_seed) gt 0.5)? 0 : 1
        dp = dps[ord[rnd_dir_ind]]
        
        if (abs(dp) lt self.turn_angle_thr) then begin
            ;return, stream->switchDirection()
        endif
        
        max_dir = reform(data[*, ord[rnd_dir_ind]])
        
    endif else begin
    
        theta = !DPI/2.0 - (*self.theta_1_samp)[no_int[0], no_int[1], no_int[2], samp_ind[0]]
        phi   = (*self.phi_1_samp)[no_int[0], no_int[1], no_int[2], samp_ind[0]]
            
        data = cv_coord(from_sphere=[ phi , theta, 1.0 ], /to_rect)
        
        dp = total(curr_dir * data) / sqrt(total(curr_dir * curr_dir))
        
        if (abs(dp) lt self.turn_angle_thr) then begin
            return, stream->switchDirection()
        endif
        
        max_dir = data
        
    endelse
    
    if (dp lt 0) then max_dir = -max_dir
    next_pos = curr_pt + max_dir * self.resolution_xyz
    
    if (self->checkBoundsROI(next_pos)) then begin
        return, stream->switchDirection()
    endif
    
    if (not self->checkAnisotropy(curr_pt)) then begin
        return, stream->switchDirection()
    endif
    
    if (self->checkTerminationHit(curr_pt)) then  begin
        return, stream->switchDirection()
    endif
    
    if (self->checkExclusionHit(curr_pt)) then begin
        return, 0 ;;stream->switchDirection()
    endif
    
    self->checkWaypointHit, curr_pt
    
    stream->addPoint, curr_pt
    stream->setPosition, next_pos
    stream->setDirection, max_dir
    
    return, 1

end

pro mas_hardi_tractography_probtrack::checkWaypointHit, coord
    
    if (not ptr_valid(self.wayp_masks)) then return

    for p = 0, n_elements(*self.wayp_masks)-1 do begin
        if (ptr_valid((*self.wayp_masks)[p]) && (*((*self.wayp_masks))[p])[coord[0], coord[1], coord[2]] ne 0) then begin
            (*self.wayp_check)[p] = 1
            ;print, "Hit!"
        endif
    endfor

    
end

function mas_hardi_tractography_probtrack::checkTerminationHit, coord

    if (not ptr_valid(self.term_masks)) then begin
        return, 0
    endif else if ( (*self.term_masks)[coord[0], coord[1], coord[2]] ne 0) then begin
        return, 1
    endif else begin
        return, 0
    endelse
end

function mas_hardi_tractography_probtrack::checkExclusionHit, coord

    if (not ptr_valid(self.excl_masks)) then begin
        return, 0
    endif else if ( (*self.excl_masks)[coord[0], coord[1], coord[2]] ne 0) then begin
        self.excl_check = 1
        return, 1
    endif else begin
        return, 0
    endelse
end

pro mas_hardi_tractography_probtrack::resetCounters, exclude=excluse, madewaypoints=madewaypoints

    ;;self.excl_check = 0
    if (ptr_valid(self.wayp_check)) then begin
        (*self.wayp_check)[*] = 0
    endif
    
end

pro mas_hardi_tractography_probtrack::getHitCounts, exclude=excluse, madewaypoints=madewaypoints

    ;;exclude = self.excl_check
    if (ptr_valid(self.wayp_check)) then begin
        madewaypoints = (total((*self.wayp_check) eq 0)) eq 0 ? 1 : 0
    endif else begin
        madewaypoints = 1
    endelse
    
end

pro mas_hardi_tractography_probtrack__define

    struct = { MAS_HARDI_TRACTOGRAPHY_PROBTRACK, $
               inherits MAS_HARDI_TRACTOGRAPHY, $
               dist_size: 0L, $
               theta_1_samp: ptr_new(), $
               theta_2_samp: ptr_new(), $
               phi_1_samp: ptr_new(), $
               phi_2_samp: ptr_new(), $
               wayp_masks: ptr_new(), $
               term_masks: ptr_new(), $
               excl_masks: ptr_new(), $
               wayp_check: ptr_new(), $
               excl_check: 0B  }

end

