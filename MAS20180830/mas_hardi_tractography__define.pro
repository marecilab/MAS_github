function mas_hardi_tractography::init, recon_obj,  $
                                       data_ptr=data_ptr, $
                                       anis_ptr=anis_ptr, $
                                       step_size=step_size, $
                                       prob_thr=prob_thr, $
                                       tess_level=tess_level, $
                                       turn_angle_thr=turn_angle_thr, $
                                       anis_halt_thr=anis_halt_thr

    common scan_data, project
    
    ci = project.ci
    
    if (keyword_set(data_ptr) && ptr_valid(data_ptr)) then begin
        self.diff_data_ptr = data_ptr
    endif else begin
        self.diff_data_ptr = project.dataarray[ci].state1
    endelse
    
    if (keyword_set(anis_ptr) && ptr_valid(anis_ptr)) then begin
        self.anis_ptr = anis_ptr
    endif else begin
        self.anis_data_ptr = project.dataarray[ci].frac_ani
    endelse
    
    self.data_is_vector_field = 0
    
    self.step_size       = keyword_set(step_size) ? float(step_size) : 0.25
    self.anis_halt_thr   = keyword_set(anis_halt_thr) ? float(anis_halt_thr) : 0.08
    self.anis_branch_thr = 0.3
    self.turn_angle_thr  = keyword_set(turn_angle_thr) ? cos(float(turn_angle_thr) * !DTOR) : cos(70.0 * !DTOR)
    self.prob_thr        = keyword_set(prob_thr) ? float(prob_thr) : 0.5
    self.tess_level      = keyword_set(tess_level) ? fix(tess_level) : 3
    if (obj_valid(recon_obj)) then  self.recon_obj = recon_obj
    self.interpolation_method = 0
    
    interp_factors = [ project.procpramarray[ci].freq_interp, $
                       project.procpramarray[ci].phase_interp, $
                       project.procpramarray[ci].slice_interp ] 

    voxel_dims = [ project.imndarray[ci].f_voxsz, $
                   project.imndarray[ci].p_voxsz, $
                   project.imndarray[ci].s_voxsz ]

    self.resolution_xyz = voxel_dims/interp_factors
    self.resolution_xyz = (self.resolution_xyz/max(self.resolution_xyz)) * self.step_size
    ;self.resolution_xyz = self.step_size
    
    if (obj_valid(recon_obj)) then begin

        mas_tessellator_make_icos, level=tess_level, vertlist=vertlist, polylist=polylist
        self.recon_vertices = ptr_new(vertlist, /no_copy)
        self.recon_polygons = ptr_new(polylist, /no_copy)
    
        (self.recon_obj)->setReconstructionDirections, *self.recon_vertices
    
        mas_tessellator_make_icos, level=0, vertlist=vertlist

        self.restart_vertices = ptr_new(vertlist)
        
    endif
    
    ;; not used
    self.tx_field_to_rawdata   = diag_matrix(double([1.0/interp_factors, 1.0]))
    self.tx_streams_to_imagery = diag_matrix(double([interp_factors, 1.0]))

    return, 1
    
end

function mas_hardi_tractography::getData, coord

    p = (floor(coord) > [0,0,0]) < (size(*self.diff_data_ptr, /DIMENSIONS)-1)
    return, (*self.diff_data_ptr)[p[0], p1[1], p2[2]]

    return, mas_tricubic_4d(coord)
    
end

function mas_hardi_tractography::initializeStreamlines, seedpt

    forward_function mas_tricubic_4d
    forward_function mas_4d_data_average
    
    if (self.data_is_vector_field) then begin
    
        nvs = 0
        vs = fltarr(3, n_elements(*self.diff_data_ptr))
        for v = 0, n_elements(*self.diff_data_ptr)-1 do begin
             ;;tmp = mas_4d_data_average(seedpt, data_ptr=(*self.diff_data_ptr)[v])
             tmp = (*(*self.diff_data_ptr)[v])[floor(seedpt[0]), floor(seedpt[1]), floor(seedpt[2]),*]
             if (total(abs(tmp)) eq 0.0) then continue
             vs[*,nvs] = tmp
             nvs++
        endfor
        if (nvs eq 0) then return, [obj_new()]
    endif else begin
        data = mas_tricubic_4d(seedpt)
        Pr = self.recon_obj->getDisplacementProbability(data=data)
        max_Pr = max(Pr, min=min_Pr)
        Pr = (Pr - min_Pr)/(max_pr - min_pr)
    
        if (total(finite(pr, /nan)) ne 0) then begin
            return, [obj_new()]
        endif
        
        self->getMaxima, Pr, vs=vs, nvs=nvs
    endelse
    
    streams = objarr(nvs)
    for s = 0L, nvs-1 do begin
        streams[s] = obj_new('mas_hardi_streamline', seedpt, vs[*,s])
    endfor
    
    return, streams
    
end

pro mas_hardi_tractography::printSummary

    print, "========= Tractography Setting Summary ========="
    print, "interpolation_method: "+string(self.interpolation_method, format='(I0)')
    print, "        vector_field: "+string(self.data_is_vector_field, format='(I0)')
    print, "      resolution_xyz: "+string(self.resolution_xyz, format='(G0)')
    print, "           step_size: "+string(self.step_size, format='(G0)')
    print, "        max_branches: "+string(self.max_branching, format='(I0)')
    print, "       anis_halt_thr: "+string(self.anis_halt_thr, format='(G0)')
    print, "      turn_angle_thr: "+string(self.turn_angle_thr, format='(G0)')+' ('+string(1.0/!DTOR*acos(self.turn_angle_thr), format='(G0)')+' degrees)'
    print, "            prob_thr: "+string(self.prob_thr, format='(G0)')
    print, "  tessellation level: "+string(self.tess_level, format='(I0)')
    print, "================================================="

end

pro mas_hardi_tractography::setProperty, $
    diff_data_ptr=diff_data_ptr, $
    anis_data_ptr=anis_data_ptr, $
    step_size=step_size, $
    anis_halt_thr=anis_halt_thr, $
    anis_branch_thr=anis_branch_thr, $
    turn_angle_thr=turn_angle_thr, $
    prob_thr=prob_thr, $
    recon_obj=recon_obj, $
    recon_vertices=recon_vertices, $
    data_is_vector_field=data_is_vector_field, $
    interpolation_method=interpolation_method, $
    max_branching=max_branching, $
    resolution_xyz=resolution_xyz

    if (keyword_set(diff_data_ptr) && ptr_valid(diff_data_ptr)) then begin
        self.diff_data_ptr = diff_data_ptr
    endif
    
    if (keyword_set(data_is_vector_field)) then begin
        self.data_is_vector_field = 1
    endif
    
    if (keyword_set(anis_data_ptr) && ptr_valid(anis_data_ptr)) then begin
        self.anis_data_ptr = anis_data_ptr
    endif

    if (keyword_set(step_size)) then begin
        tmp = self.resolution_xyz / self.step_size
        self.step_size = float(step_size)
        self.resolution_xyz = tmp * self.step_size
    endif
    
    if (keyword_set(anis_halt_thr)) then begin
        self.anis_halt_thr = float(anis_halt_thr)
    endif
    
    if (keyword_set(anis_branch_thr)) then begin
        self.anis_branch_thr = float(anis_branch_thr)
    endif
    
    if (keyword_set(turn_angle_thr)) then begin
        self.turn_angle_thr = cos(turn_angle_thr * !DTOR)
    endif
    
    if (keyword_set(prob_thr)) then begin
        self.prob_thr = prob_thr
    endif
    
    if (keyword_set(recon_obj) && obj_valid(recon_obj)) then begin
        self.recon_obj = recon_obj
    endif

    if (keyword_set(recon_vertices) && ptr_valid(recon_vertices) ne 0) then begin
        self.recon_vertices = ptr_new(*recon_vertices)
    endif
    
    if (keyword_set(interpolation_method) or arg_present(interpolation_method)) then begin
        self.interpolation_method = interpolation_method
    endif
    
    if (keyword_set(resolution_xyz)) then begin
        self.resolution_xyz = resolution_xyz
    end

    if (arg_present(max_branching)) then begin
        self.max_branching = max_branching
    endif
end

pro mas_hardi_tractography::setBoundsROI, roi

    if (not obj_valid(roi)) then return
    
    obj_destroy, self.bound_roimask

    self.bound_roimask = roi

end

pro mas_hardi_tractography::setExclusionROI, roi

    if (not obj_valid(roi)) then return
    
    obj_destroy, self.exclusion_roimask

    self.exclusion_roimask = roi

end

pro mas_hardi_tractography::addROI, mask, roitype=roitype

    ;; seed roi
    sdpts = where(*mask ne 0, n_sdpts)
    if (n_sdpts gt 0) then begin
        sdpts = array_indices(*mask, sdpts)
        self.active_seedpts = ptr_new(sdpts, /no_copy)
        self.active_roimask = ptr_new(*mask)
    endif

end

pro mas_hardi_tractography::track

    forward_function mas_tricubic_4d
    ;; for each seed point
    nseeds = n_elements(*self.active_seedpts)/3

    for sd = 0, nseeds-1 do begin

        sdpt = (*self.active_seedpts)[*,sd]
        data = mas_tricubic_4d(sdpt)
        Pr = self.recon_obj->getDisplacementProbability(data=data)
        max_pr = (where(Pr eq max(Pr)))[0]
        self->newStream, sdpt, (*self.recon_vertices)[*,max_pr]

        status = 1
        while status ne 0 do begin
            status = self->jump()
        endwhile

        if (not ptr_valid(self.completed_streams)) then begin
            self.completed_streams = ptr_new([ self.active_streamline ])
        endif else begin
            (*self.completed_streams) = [ *self.completed_streams, self.active_streamline ]
        endelse

    endfor
    
    help, *self.completed_streams
    
end

pro mas_hardi_tractography::newStream, seed_point, seed_direction

    print, "mas_hardi_tractography::newStream", seed_point
    ;; create a new streamline object from current seed point
    self.active_streamline = obj_new('mas_hardi_streamline', seed_point, seed_direction)
    
end

pro mas_hardi_tractography::getMaxima, Pr, vs=vs, ps=ps, nvs=nvs, maxima_mask=maxima_mask

    nPr = n_elements(Pr)
    n_rsverts = n_elements(*self.restart_vertices)/3
    pr_interp_tol = 0.25
    
    vs = fltarr(3, n_rsverts)
    ps = fltarr(n_rsverts)
    ;;maxima_mask = bytarr(n_elements(*self.recon_vertices)/3)
    
    for rv = 0, n_rsverts-1 do begin
       
       ;; Select a vertex from the "restart" vertices; call it "test"
       test = reform((*self.restart_vertices)[*,rv])
       
       ;; find neighbors of "test" (within pr_interp_tol)
       ;; call these the "potential coordinates"
       pot_coords = where(abs((*self.recon_vertices)[0,*] - test[0]) lt pr_interp_tol and $
                          abs((*self.recon_vertices)[1,*] - test[1]) lt pr_interp_tol and $
                          abs((*self.recon_vertices)[2,*] - test[2]) lt pr_interp_tol, count)
       
       ;; proceed in the direction of greatest increase in Pr
       ;; assumes that the probability function is reasonably smooth
       ;; (which is it in this case)
       curr_pr = 0.0
       iter = 0L
       ;;best_coord_index = -1L
       while (1) do begin
          if (count eq 0) then begin
             print, "No vertices returned"
             break
          endif else begin
             
             ;; get the probabilites associated with the pot_coords
             ;; and choose the max out of all of them
             possible_pr = Pr[pot_coords]
             max_pr = max(possible_pr)
             possible_coord = (pot_coords[where(possible_pr eq max_pr)])[0]
             
             ;; if this probabilty is greater than current max, save
             ;; it + the coords
             if (max_pr gt curr_pr) then begin
                test = (*self.recon_vertices)[*, possible_coord]
                curr_pr = max_pr
                ;;best_coord_index = possible_coord
                iter++
             endif else begin
                ;; break out of the loop when the probabilites begin
                ;; to decrease
                iter++
                break
             endelse
          endelse

          ;; now find the neighbors for the next iteration
          pot_coords = where(abs((*self.recon_vertices)[0,*] - test[0]) lt pr_interp_tol and $
                              abs((*self.recon_vertices)[1,*] - test[1]) lt pr_interp_tol and $
                              abs((*self.recon_vertices)[2,*] - test[2]) lt pr_interp_tol, count)
          
       endwhile

       ;; at the end of the while loop, "test" will have the vector
       ;; associated witht he maximum probability 
       ;;maxima_mask[best_coord_index] = 1B
       vs[*,rv] = test ;;total(test * [0,0,1]) lt 0 ? -test : test
       ps[rv] = max_pr
       
    endfor
    
    eligible = where(ps gt self.prob_thr)
    
    vs = vs[*,eligible]
    ps = ps[eligible]

    ;print, vs ## transpose(vs)
    ;print, ''
    
    dps = ((vs ## transpose(vs)))[0,*]
    negs = where(dps lt 0, count)
    if (count gt 0) then begin
        vs[*, negs] *= (-1)
        dps[negs] *= (-1)
    endif
    
    dps_sorted = sort(dps)    
    dps = dps[dps_sorted]
    vs  = vs[*,dps_sorted]
    ps  = ps[dps_sorted]
    
    dps_uniq = uniq(dps)
    vs = vs[*,dps_uniq]
    ps = ps[dps_uniq]

    nvs = n_elements(dps_uniq)

end

function mas_hardi_tractography::checkAnisotropy, point

    if (not ptr_valid(self.anis_data_ptr)) then return, 1
    
    result = ((*self.anis_data_ptr)[point[0], point[1], point[2]] gt self.anis_halt_thr)
    return, result

end

function mas_hardi_tractography::checkBoundsData, pt_in

    if (self.data_is_vector_field) then begin
        sz = size((*(*self.diff_data_ptr)[0]), /dimensions)
        ;sz = (size(*self.diff_data_ptr, /dimensions))[1:3]
    endif else begin
        sz = size(*self.diff_data_ptr, /dimensions)
    endelse
    pt = round(pt_in)
    if (pt[0] le 0 or pt[1] le 0 or pt[2] le 0) or $
        (pt[0] ge sz[0] or pt[1] ge sz[1] or pt[2] ge sz[2]) then return, 1
    
    return, 0

end

function mas_hardi_tractography::checkBoundsROI, curr_pt

    if (not obj_valid(self.bound_roimask)) then return, 0
    return, (self.bound_roimask)->containsPoint(curr_pt) 

end

function mas_hardi_tractography::checkExclusionROI, curr_pt

    if (not obj_valid(self.exclusion_roimask)) then return, 0
    return, (self.exclusion_roimask)->containsPoint(curr_pt) 

end

function mas_hardi_tractography::checkTurnAngle, curr_pt, next_pt

    print, "mas_hardi_tractography::checkTurnAngle"
    return, total(curr_pt * next_pt) 

end

function mas_hardi_tractography::jumpVectorField, stream, glyph_out=glyph_out, make_glyph=make_glyph, split=split

    forward_function mas_tricubic_4d
    forward_function mas_4d_data_average
    
    ;;if (self.data_is_vector_field) then return, self->jumpVectorField()
    
    curr_pt = stream->getCurrentPoint(direction=curr_dir)
    ;print, "mas_hardi_tractography::jump", curr_pt

    if (self->checkBoundsData(curr_pt)) then begin
        return, stream->switchDirection()
    endif
    
    if (not self->checkAnisotropy(floor(curr_pt))) then begin
        return, stream->switchDirection()
    endif
    vs = fltarr(3, n_elements(*self.diff_Data_ptr))
    ps = fltarr(n_elements(*self.diff_Data_ptr))
    nvs = 0
    case self.interpolation_method of 
    
        0: begin ;; nearest neighbor
            no_int = (floor(curr_pt) < (size((*(*self.diff_data_ptr)[0]), /dimensions)-1)) > [0,0,0]
            for v = 0, n_elements(*self.diff_data_ptr)-1 do begin
                tmp = (*((*self.diff_data_ptr)[v]))[no_int[0],no_int[1],no_int[2], *]
                if (total(abs(tmp)) eq 0.0) then continue
                vs[*,nvs] = tmp
                ps[nvs] = 1.0
                nvs++
            endfor
        end
        
        1: begin ;; tricubic -- BROKEN when using vector field
            for v = 0, n_elements(*self.diff_data_ptr)-1 do begin
                tmp = mas_tricubic_4d(curr_pt, data_ptr=(*self.diff_data_ptr)[v])
                if (total(abs(tmp)) eq 0.0) then continue
                vs[*,nvs] = tmp
                ps[nvs] = 1.0
                nvs++
            endfor
        end
        
        2: begin ;; analogous to trilinear -- BROKEN when using vector field
            for v = 0, n_elements(*self.diff_data_ptr)-1 do begin
                tmp = mas_4d_data_average(curr_pt, data_ptr=(*self.diff_data_ptr)[v])
                if (total(abs(tmp)) eq 0.0) then continue
                vs[*,nvs] = tmp
                ps[nvs] = 1.0
                nvs++
            endfor
            
        end
        
    endcase
    
    if (1) then begin
        
        if (nvs gt 1) then begin
            dp = (vs ## [ [replicate(curr_dir[0], nvs)], [replicate(curr_dir[1], nvs)], [replicate(curr_dir[2], nvs)] ])[0,*]
            neg_dps = where(dp lt 0, n_neg_dps)
            if (n_neg_dps ne 0) then begin
                vs[*,neg_dps] *= (-1.0)
            endif
            abs_dp = abs(dp)
            max_dp = max(abs_dp, min=min_dp)
            sel = (where(abs_dp eq max_dp))[0]
            alt = (where(abs_dp ne max_dp, n_alt))
            max_dir = vs[*, sel]
            if (n_alt ne 0) then begin
                alt_dir = vs[*, alt]
            endif
            dps = dp
        endif else begin
            max_dir = vs[*,0]
        endelse
        
    endif else begin
    
        max_ind = (where(Pr eq max(Pr)))[0]
        max_dir = (*self.recon_vertices)[*,max_ind]
        
    endelse
    
    dp = total(max_dir * curr_dir) / sqrt(total(curr_dir * curr_dir))
    if (dp lt 0) then max_dir = -max_dir
    
    if (abs(dp) lt self.turn_angle_thr) then begin
        return, stream->switchDirection()
    endif
    
    next_pos = curr_pt + max_dir * self.resolution_xyz
    ;if (stream->compareAngleWithPivot(next_pos, self.turn_angle_thr) eq 0) then begin
    ;    return, stream->switchDirection()
    ;endif

    if (self->checkBoundsROI(next_pos)) then begin
        return, stream->switchDirection()
    endif else if (self->checkExclusionROI(next_pos)) then begin
        stream->setExcluded, 1
        junk = stream->switchDirection()
        return, stream->switchDirection()
    endif
    
    stream->addPoint, curr_pt
    stream->setPosition, next_pos
    stream->setDirection, max_dir
    
    if (self.max_branching ne 0) then begin
        if (n_elements(alt) ne 0 && n_alt ne 0) then begin
            nsplit = 0
            for a = 0, n_elements(alt)-1 do begin
                if (dps[alt[a]] le self.turn_angle_thr) then continue
                sp = stream->duplicate()
                sp->addPoint, curr_pt
                sp->setPosition, curr_pt + alt_dir[*,a] * self.resolution_xyz
                sp->setDirection, alt_dir[*,a]
                tmp_split = (n_elements(tmp_split) eq 0) ? sp : [ tmp_split, sp ]
                nsplit++
            endfor
            if (nsplit ne 0) then split = temporary(tmp_split)
        endif
    endif
    
    return, 1

end

function mas_hardi_tractography::jump, stream, glyph_out=glyph_out, make_glyph=make_glyph, split=split

    forward_function mas_tricubic_4d
    forward_function mas_4d_data_average
    
    if (self.data_is_vector_field) then return, self->jumpVectorField(stream, split=split)
    ;;if (self.data_is_vector_field) then return, self->jumpVectorField(stream, split=(n_elements(split) ? 1 : 0))
    
    curr_pt = stream->getCurrentPoint(direction=curr_dir)
    ;print, "mas_hardi_tractography::jump", curr_pt

    if (self->checkBoundsData(curr_pt)) then begin
        return, stream->switchDirection()
    endif
    
    if (not self->checkAnisotropy(curr_pt)) then begin
        return, stream->switchDirection()
    endif
    
    case self.interpolation_method of 
    
        0: begin ;; nearest neighbor
            no_int = (round(curr_pt) < (size(*self.diff_data_ptr, /dimensions)-1)) > [0,0,0]
            data = reform((*self.diff_data_ptr)[no_int[0], no_int[1], no_int[2], *])
        end
        
        1: begin ;; tricubic
            data = mas_tricubic_4d(curr_pt, data_ptr=self.diff_data_ptr)
        end
        
        2: begin ;; analogous to trilinear
            data = mas_4d_data_average(curr_pt, data_ptr=self.diff_data_ptr)
        end
        
    endcase
    
    Pr = self.recon_obj->getDisplacementProbability(data=data)
    max_Pr = max(Pr, min=min_Pr)
    Pr = (Pr - min_Pr)/(max_pr - min_pr)
    
    if (total(finite(Pr, /nan)) ne 0) then begin
        print, "mas_hardi_tractography::jump(): NAN in probability reconstruction at: "+$
            "["+strjoin(string(curr_pt, format='(G0)'), ',')+"]."
        return, stream->switchDirection()
    endif else if (total(abs(Pr)) eq 0) then begin
        ;; all zeros means something went wrong in computation
        print, "mas_hardi_tractography::jump(): probability reconstruction failed at: "+$
            "["+strjoin(string(curr_pt, format='(G0)'), ',')+"]."
        return, stream->switchDirection()
    endif

;    if (keyword_set(make_glyph)) then begin
;        vertlist = *self.recon_vertices
;        vertlist = (vertlist * transpose([[Pr*0.5], [Pr*0.5], [Pr*0.5]])) + $
;                    transpose([ [replicate(curr_pt[0],642)], [replicate(curr_pt[1],642)], [replicate(curr_pt[2],642)] ])
;                    
;        glyph_out = obj_new('idlgrpolygon', vertlist, polygons=*self.recon_polygons, color=[0,0,255], alpha_channel=0.5)
;    endif
    
    if (1) then begin

        self->getMaxima, Pr, vs=vs, nvs=nvs, ps=ps
        
        if (nvs gt 1) then begin
            dp = (vs ## [ [replicate(curr_dir[0], nvs)], [replicate(curr_dir[1], nvs)], [replicate(curr_dir[2], nvs)] ])[0,*]
            neg_dps = where(dp lt 0, n_neg_dps)
            if (n_neg_dps ne 0) then begin
                vs[*,neg_dps] *= (-1.0)
            endif
            abs_dp = abs(dp)
            max_dp = max(abs_dp, min=min_dp)
            sel = (where(abs_dp eq max_dp))[0]
            alt = (where(abs_dp ne max_dp, n_alt))
            max_dir = vs[*, sel]
            if (n_alt ne 0) then begin
                alt_dir = vs[*, alt]
            endif
            dps = dp
        endif else begin
            max_dir = vs
        endelse
        
    endif else begin
    
        max_ind = (where(Pr eq max(Pr)))[0]
        max_dir = (*self.recon_vertices)[*,max_ind]
        
    endelse
    ;; note that max_dir is already unit vector 
    dp = total(max_dir * curr_dir) / (sqrt(total(curr_dir * curr_dir)));; * sqrt(total(max_dir * max_dir)))
    if (dp lt 0) then max_dir = -max_dir

    if (abs(dp) lt self.turn_angle_thr) then begin
        return, stream->switchDirection()
    endif
;    if (finite(acos(abs(dp) < 1.0) * 1.0/!DTOR, /nan)) then begin
;        print, "NAN!"
;    endif
        
    next_pos = curr_pt + max_dir * self.resolution_xyz
    
    if (self->checkBoundsROI(next_pos)) then begin
        return, stream->switchDirection()
    endif else if (self->checkExclusionROI(next_pos)) then begin
        stream->setExcluded, 1
        junk = stream->switchDirection()
        return, stream->switchDirection()
    endif
    
    stream->addPoint, curr_pt
    stream->setPosition, next_pos
    stream->setDirection, max_dir
    
    if (self.max_branching ne 0 && n_elements(alt) ne 0 && n_alt ne 0) then begin
        nsplit = 0
        for a = 0, n_elements(alt)-1 do begin
            if (dps[alt[a]] le cos(70.0 * !DTOR)) then continue
            sp = stream->duplicate()
            sp->addPoint, curr_pt
            sp->setPosition, curr_pt + alt_dir[*,a] * self.resolution_xyz
            sp->setDirection, alt_dir[*,a]
            tmp_split = (n_elements(tmp_split) eq 0) ? sp : [ tmp_split, sp ]
            nsplit++
        endfor
        if (nsplit ne 0) then split = temporary(tmp_split)
    endif
    
    return, 1
        
end

pro mas_hardi_tractography::trackseed, sdpoint, streams=streams

    if (n_elements(sdpoint) eq 0) then return
    
    streams = self->initializeStreamlines(sdpoint)
    valid = where(obj_valid(streams), n_valid)
    if (n_valid eq 0) then return
    stream = streams[valid]
    
    ;; For each eligible stream in the seed point
    for s = 0L, n_elements(stream)-1 do begin
    
        status = -1
        iter   = 0L
        while (status ne 0) do begin
            status = self->jump(stream[s])
            iter++
        endwhile
        
        print, string(iter, format='(I0)')+ ' jumps'
                
    endfor ;; stream
    
end

pro mas_hardi_tractography::cleanup

    ;obj_destroy, self.recon_obj
    obj_destroy, self.active_streamline
    
    ptr_free, self.recon_vertices
    ptr_free, self.recon_polygons
    ptr_free, self.restart_vertices
    ptr_free, self.active_roimask
    ptr_free, self.active_seedpts
    ptr_free, self.roi_list
    ptr_free, self.roi_types
    
    if (ptr_valid(self.completed_streams)) then begin
        obj_destroy, *self.completed_streams
    endif

end

pro mas_hardi_tractography__define 

    struct = { MAS_HARDI_TRACTOGRAPHY, $
               random_seed: 0L, $
               diff_data_ptr: ptr_new(), $
               anis_data_ptr: ptr_new(), $
               interpolation_method: 0B, $
               data_is_vector_field: 0B, $
               resolution_xyz: fltarr(3), $
               step_size: fltarr(3), $
               tx_field_to_rawdata: dblarr(4,4), $
               tx_streams_to_imagery: dblarr(4,4), $
               max_branching: 0B, $
               anis_halt_thr: float(0), $
               anis_branch_thr: float(0), $
               turn_angle_thr: float(0), $
               prob_thr: float(0), $
               tess_level: 0, $
               recon_obj: obj_new(), $
               recon_vertices: ptr_new(), $
               recon_polygons: ptr_new(), $
               restart_vertices: ptr_new(), $
               bound_roimask: obj_new(), $
               exclusion_roimask: obj_new(), $
               active_roimask: ptr_new(), $
               active_seedpts: ptr_new(), $
               active_streamline: obj_new(), $
               roi_list: ptr_new(), $
               roi_types: ptr_new(), $
               completed_streams: ptr_new() }
               
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mas_hardi_streamline::init, seed_point, seed_direction

    self.propagation_dir = 0
    
    self.seed_point = seed_point
    self.seed_direction = seed_direction

    self.curr_point = seed_point
    self.curr_direction = seed_direction
    
    self.points_forward = ptr_new(seed_point)
    self.points_backward = ptr_new(seed_point)

    return, 1
end

function mas_hardi_streamline::switchDirection

    if (self.propagation_dir eq 0) then begin
        self.curr_point = self.seed_point
        self.curr_direction = -self.seed_direction
        self.propagation_dir = 1
        return, 1
    endif else begin
        return, 0
    endelse

end

pro mas_hardi_streamline::addPoint, point

    if (self.propagation_dir eq 0) then begin
        *(self.points_forward) = [ [*(self.points_forward)], [point] ]
    endif else begin
        *(self.points_backward) = [ [*(self.points_backward)], [point] ]
    endelse
    
end

function mas_hardi_streamline::getCurrentPoint, direction=direction

    direction = self.curr_direction
    return,  self.curr_point
    
end

function mas_hardi_streamline::getCurrentDirection

    return, self.curr_direction
    
end

pro mas_hardi_streamline::setPosition, position

    self.curr_point = position
    
end

pro mas_hardi_streamline::setDirection, direction

    self.curr_direction = direction
    
end

pro mas_hardi_streamline::setExcluded, bool

    self.is_excluded = (n_elements(bool) ne 0 && bool ne 0) ? 1 : 0
    
end

function mas_hardi_streamline::isExcluded

    return, self.is_excluded

end

function mas_hardi_streamline::compareAngleWithPivot, point, cos_angle

    vec_pivot2pt  = point - self.pivot_point
    vec_pivot2pt  = vec_pivot2pt/(sqrt(total(vec_pivot2pt^2)))
    dp = total(self.pivot_direction * vec_pivot2pt)
   
    dist_pivot2pt = sqrt( total( (point - self.pivot_point)^2 ) )
    ;print, dp, acos(dp)*1.0/!DTOR, dist_pivot2pt
    
    if (dist_pivot2pt gt 1) then begin
        self.pivot_direction = self.curr_direction
        self.pivot_point = self.curr_point
    endif

    if (dp lt cos_angle) then return, 0
    
    return, 1
end


function mas_hardi_streamline::getTerminalPoints, tol

    if (n_elements(tol) eq 0) then tol = 2
    
    n_fwd = n_elements(*self.points_forward)
    n_bwd = n_elements(*self.points_backward)
    
    end1 = (*self.points_forward)[*,n_fwd-tol:n_fwd-1]
    end2 = (*self.points_forward)[*,n_bwd-tol:n_bwd-1]

end

function mas_hardi_streamline::getStreamAsPoints, count=count

    bwd = *self.points_backward
    if (n_elements(*self.points_forward)/3 gt 1) then begin
        fwd = reverse(*self.points_forward,2)
    endif else begin
        fwd = *self.points_forward
    endelse
    
    nfwd = n_elements(fwd)/3
    nbwd = n_elements(bwd)/3
    if (nfwd gt 2) then begin
        fwd = fwd[*,0:nfwd-2]
    endif
    if (nbwd gt 2) then begin
        bwd = bwd[*,2:nbwd-1]
    endif
    
    points = [ [fwd], [bwd] ]
    count = n_elements(points)/3
    
    return, points
    
end

pro mas_hardi_streamline::getProperty,       $
            propagation_dir=propagation_dir, $
            seedpt=seedpt, seeddir=seeddir, $
            midpt_idx=midpt_idx, $
            currpt=currpt, currdir=currdir

    propagation_dir = self.propagation_dir
    seedpt = self.seed_point
    seeddir = self.seed_direction
    currpt = self.curr_point
    currdir = self.curr_direction
    midpt_idx = n_elements(*self.points_forward)/3

end

pro mas_hardi_streamline::cleanup

    ptr_free, self.points_forward
    ptr_free, self.points_backward

end

function mas_hardi_streamline::duplicate
    
     dupe = obj_new('mas_hardi_streamline', self.seed_point, self.seed_direction)
     dupe.propagation_dir = self.propagation_dir
     dupe.seed_point = self.seed_point
     dupe.curr_point = self.curr_point
     dupe.seed_direction = self.seed_direction
     dupe.curr_direction = self.curr_direction
     if (ptr_valid(self.points_forward)) then begin
        dupe.points_forward = ptr_new(*self.points_forward)
     endif
     
     if (ptr_valid(self.points_backward)) then begin
        dupe.points_backward = ptr_new(*self.points_backward)
     endif
     
     return, dupe
     
end

pro mas_hardi_streamline__define

    struct = { MAS_HARDI_STREAMLINE, $
               is_excluded: 0B, $
               propagation_dir: 0B, $
               seed_point: fltarr(3), $
               curr_point: fltarr(3), $
               seed_direction: fltarr(3), $
               curr_direction: fltarr(3), $
               points_forward: ptr_new(), $
               points_backward: ptr_new() }

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mas_hardi_ROI::init, mask, name=name, type=type

    self.mask = ptr_new(*mask)
    self.mask_dim = size(*mask, /dimensions)
    
    self.name = (keyword_set(name)) ? string(name) : 'Untitled'
    self.type = (keyword_set(type)) ? self->getType(type) : 0
    
    ind = reverse(( (*mask)[ sort((*mask)[*]) ] ))
    ind = ind[uniq(ind)]
    if (n_elements(ind) gt 1) then begin
        ind = reverse(ind[0:n_elements(ind)-2])
    endif
    self.subregions = ptr_new(ind, /no_copy)

    return, 1
end

function mas_hardi_ROI::getSubRegions

    if (ptr_valid(self.subregions)) then begin
        return, *self.subregions
    endif else begin
        return, -1
    endelse

end

function mas_hardi_ROI::makeSeeds, roi_num, density=density, exponent=exponent, randomize=randomize

    sz_mask = self.mask_dim
    
    if (keyword_set(randomize)) then begin
    
        if (n_elements(roi_num) ne 0) then begin
            active = where(*self.mask eq roi_num, n_active)
        endif else begin
            active = where(*self.mask ne 0, n_active)
        endelse
        
        if (n_active eq 0) then return, -1
        
        active_ind = array_indices(*self.mask, active)
        seed_rnd = systime(1)
        seeds = fltarr(3, n_active*randomize)
        ct = 0L
        for ind = 0L, n_active-1 do begin
            for r = 0L, randomize-1 do begin
                seeds[*,ct++] = active_ind[*,ind] + randomu(seed_rnd, 3)
            endfor
        endfor
        
        return, seeds
        
    endif else begin
    
        density  = (keyword_set(density))  ? density  : 1
        exponent = (keyword_set(exponent)) ? exponent : 3
        
        if (n_elements(roi_num) ne 0) then begin
            active = where(*self.mask eq roi_num, n_active)
        endif else begin
            active = where(*self.mask ne 0, n_active)
        endelse
        
        if (n_active eq 0) then return, -1
        
        active_ind = array_indices(*self.mask, active)
        
        if (density eq 1) then return, active_ind + 0.5
        
        work = bytarr(density,density,density)
        
        sds_per_vxl = density^exponent
        
        seeds = fltarr(3, n_elements(active_ind)/3 * sds_per_vxl)
        
        oldPT = !P.t
        
        t3d, /reset, scale=double([1,1,1])/density
        sd = 0L
        vox = array_indices(work, indgen(sds_per_vxl))
        for n = 0L, n_elements(active_ind)/3 - 1 do begin
                    
            tmp = vert_t3d(vox)
            
            tmp[0,*] += active_ind[0,n]
            tmp[1,*] += active_ind[1,n]
            tmp[2,*] += active_ind[2,n]
            
            seeds[*, sd:sd+sds_per_vxl - 1] = tmp
            sd += sds_per_vxl
            ;seeds = (n eq 0) ? vox : [ [seeds], [vox] ]
            
        endfor
        
        !P.t = oldPT
        
        seeds += (1.0/density/2.0)
        return, seeds

    endelse
    

end

function mas_hardi_ROI::getType, strcode

    case strcode of 
        'SEEDING' : return, 0
        'CATCHING' : return, 1
        'BOUNDING' : return, 2
        else: return, 0
    endcase
    
end

pro mas_hardi_ROI::computeMesh, vertlist, polylist

    if (ptr_valid(self.mask)) then begin
        shade_volume, *self.mask, 0, vertlist, polylist, /low
    endif

end

function mas_hardi_ROI::containsPoint, point

    pt = round(point)
    pt[0] = ((round(point[0])) > 0) < (self.mask_dim[0]-1)
    pt[1] = ((round(point[1])) > 0) < (self.mask_dim[1]-1)
    pt[2] = ((round(point[2])) > 0) < (self.mask_dim[2]-1)
    return, ((*self.mask)[pt[0], pt[1], pt[2]] ne 0)

end

function mas_hardi_ROI::intersectsStream, stream, subregion=subregion

    pts = stream->getStreamAsPoints()
    
    for i = 0, n_elements(pts)/3-1 do begin
    
        pt = round(pts[*,i])
        
        if (keyword_set(subregion)) then begin
            if ((*self.mask)[pt[0],pt[1],pt[2]] eq subregion) then begin
                return, 1
            endif
        endif else begin
            if ((*self.mask)[pt[0],pt[1],pt[2]] ne 0) then begin
                return, 1
            endif
        endelse
        
    endfor

    return, 0

end

pro mas_hardi_ROI::cleanup

    ptr_free, self.mask
    ptr_free, self.seeds
    ptr_free, self.subregions
    
end

pro mas_hardi_ROI__define

    struct = { MAS_HARDI_ROI, $
               subregions: ptr_new(), $
               mask_dim: lonarr(3), $
               mask: ptr_new(), $
               seeds: ptr_new(), $
               name: '', $
               type: byte(0) }
               
end
