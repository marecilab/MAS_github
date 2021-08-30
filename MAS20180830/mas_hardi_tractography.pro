;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_seedROIs
;; Created by:
;; Calling Information:   
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;
;; Editing Information:

function mas_hardi_tractography_seedROIs, roiptrs, density=density, randomize=randomize

    n_rois = n_elements(roiptrs)
    if (n_rois eq 0) then return, -1
    
    for roi = 0, n_rois-1 do begin
    
        if (not ptr_valid(roiptrs[roi])) then continue
        
        oroi = obj_new('mas_hardi_ROI', roiptrs[roi], name='ROI_'+string(roi, format="(I0)"))
        
;        subs =  orois[roi]->getSubRegions()
        
;        for i = 0, n_elements(subs)-1 do begin
;            tmp_seeds = orois[roi]->makeSeeds(subs[i], density=1)
;            if (tmp_seeds[0] eq -1) then continue
;            roi_arr = (nrois++ eq 0) ? ptr_new(tmp_seeds) : [ roi_arr, ptr_new(tmp_seeds) ]
;        endfor
        
        if (keyword_set(randomize)) then begin
            tmp = oroi->makeSeeds(randomize=density)
        endif else begin
            tmp = oroi->makeSeeds(density=density, exponent=3)        
        endelse

        obj_destroy, oroi
        
        seeds = ( n_elements(seeds) eq 0 ) ? temporary(tmp) : [ [seeds], [temporary(tmp)] ]

    endfor
        
    seedpoints = ptr_new(seeds, /no_copy)
     
    return, seedpoints
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_preprocess
;; Created by: BT, 2008-07
;; Calling Information:
;; 
;; method_spec   IN: the reconstructor
;; threshold     IN: data threshold in percent
;; even/odd      IN: set one of these two process only even or odd slices
;; ndirections   IN: the number fo directions to resolve. Results in 
;;                   ndirections number of output files
;; destdir       IN: the destination directory in which to save the output files
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; Processes an entire data set for HARDI directionality. 
;;
;; Editing Information:
;;
pro mas_hardi_tractography_preprocess, method_spec, $
                                prob_thr=prob_thr, $
                                data_ptr=data_ptr, $
                                threshold=threshold, $
                                even=even, odd=odd, $
                                ndirections=ndirections, $
                                tess_level=tess_level, $
                                destdir=destdir

    common scan_data, project
    common probtrack_objective, Pr, vlist, pr_interp_tol

    if (not obj_valid(method_spec)) then return

    if (n_elements(pr_interp_tol) eq 0) then pr_interp_tol = 0.350

    if (not keyword_set(normalize_probabilities)) then begin
        normalize_probabilities = 0
    endif

    if (not keyword_set(ndirections)) then begin
        ndirections = 4
    endif
    
    if (not keyword_set(prob_thr)) then begin
        prob_thr = 0.5
    endif 
    
    if (not keyword_set(tess_level)) then begin
        tess_level = 3
    endif
    
    threshold  = n_elements(threshold) ne 0 ? float(threshold) : 1.5

    if (keyword_set(data_ptr) and ptr_valid(data_ptr)) then begin
        data = data_ptr
    endif else begin
        ci = project.ci
        data = project.dataarray[ci].state1
    endelse

    ht = obj_new('mas_hardi_tractography', method_spec, data_ptr=data, prob_thr=prob_thr, tess_level=tess_level)
    
    data_sz = size((*data), /dimensions)

    fdim = data_sz[0]
    pdim = data_sz[1]
    sdim = data_sz[2]
    adim = data_sz[3]
    
    if (keyword_set(destdir)) then begin
        if (not file_test(destdir, /directory, /write)) then begin
            void = dialog_message(['Cannot write to:', destdir], /error, /center)
            return
        endif
    endif else begin
        destdir = dialog_pickfile(path=project.current_path, /directory, /must_exist, /write, $
                                title='Please choose a destination directory.')
        if (destdir eq '') then return
    endelse

    pbar = obj_new('progressbar', $
                   text='Processing Directionality', $
                   title='Processing', /fast_loop)
    pbar->start

    ;;mask = mas_roi_get_current_mask([fdim, pdim], crop=crop_dims, /no_transform)
    
    pbar->setProperty, text='Processing Threshold...'
    max_signal = max(*data)
    pct_threshold = max_signal * threshold/100.0

    pbar->setProperty, text='Processing Diffusion Data'
    
    dir_ptrs = ptrarr(ndirections)
    for v = 0,ndirections-1 do begin
        dir_ptrs[v] = ptr_new(fltarr(fdim,pdim,sdim,3))
    endfor
    
    iter = 0.
    ;;n  = n_elements(vlist)/3 

;;;; The following is an alternative way to keep track of maxima.
;;;; Given a fixed set of n possible reconstrction directions, we may
;;;; keep track of the maxima of Pr using an array of n bytes, setting
;;;; byte[i] = 1 if direction[i] is a maxima, byte[i] = 0 otherwise. 
;;;; this method requires a great deal more space, and is not used at this time.
;    maxima_data = bytarr(n, fdim, pdim, sdim)
;    maxima_hdr = create_struct(name='MAS_NIFTI_HEADER')
;    maxima_hdr.dim = fix([4, n, fdim, pdim, sdim, 1, 1, 1])
;    maxima_hdr.pixdim = float([1.0, 1,1,1,1,1,1,1])
;    maxima_hdr.datatype = 2
;    maxima_hdr.bitpix = 8
;    maxima_hdr.sizeof_hdr = 348
;    maxima_hdr.magic = byte('n+1')
;    maxima_hdr.vox_offset = 352
;    maxima_hdr.scl_slope = 1.0
;    maxima_hdr.scl_inter = 0
;    maxima_hdr.xyzt_units = 10
;    maxima_hdr.qform_code = 0
;    maxima_hdr.sform_code = 1
;    maxima_hdr.srow_x = [1.0,0.0,0.0,0.0]
;    maxima_hdr.srow_y = [0.0,1.0,0.0,0.0]
;    maxima_hdr.srow_z = [0.0,0.0,1.0,0.0]
;    openw, lun, destdir+'hardi_directional_maxima.nii', /get_lun
;    writeu, lun, maxima_hdr
;    point_lun, lun, maxima_hdr.vox_offset
    
    st = systime(1)
    
    ;; process even or odd slices -- cheap multithreading by 
    ;; loading the data set in two instances of MAS and having one
    ;; process even, one odd slices, then add the resulting data files
    if (keyword_set(even)) then begin
        start = 0
        skip = 2
        f_append = 'even'
    endif else if (keyword_set(odd)) then begin
        start = 1
        skip = 2
        f_append = 'odd'
    endif else begin
        start = 0
        skip = 1
        f_append = 'all'
    endelse
    
    for ss = start, sdim-1, skip do begin
       print, "Starting Slice", ss
        for pp = 0, pdim-1 do begin
            iter++
            for rr = 0, fdim-1 do begin

                voxel_data = reform((*data)[rr,pp,ss,*])
                ;maxima_mask = bytarr(n)
                
                if (voxel_data[0] lt pct_threshold) then begin
                    goto, WRITE_MAXIMA
                endif
                if (ptr_valid(mask)) then begin
                    if (*mask)[rr,pp] eq 0 then begin
                        goto, WRITE_MAXIMA
                    endif
                endif

                Pr = method_spec->getDisplacementProbability(data=voxel_data)

                nans = finite(Pr, /nan)
                if (total(nans) ne 0) then begin
                    goto, WRITE_MAXIMA
                endif
                
                max_pr = max(Pr, min=min_pr)
                if (max_pr eq min_pr) then begin
                    goto, WRITE_MAXIMA
                endif
                
                Pr = (Pr - min_pr)/(max_pr - min_pr)
;;;;;;;;;;;;;;;;;;;;
                ht->getMaxima, Pr, vs=vs, ps=ps, nvs=nvs
                ;mas_hardi_tractography_get_maxima, Pr, vlist_restart, 0.65, $
                ;    vs=vs, ps=ps, nvs=nvs, maxima_mask=maxima_mask
                                
                for v = 0, ((nvs-1)<(ndirections-1)) do begin
                    dir = reform(vs[*,v])
                    if (dir[2] lt 0) then dir *= (-1.0)
                    (*dir_ptrs[v])[rr,pp,ss,*] = dir/(sqrt(total(dir*dir)))
                endfor

                WRITE_MAXIMA:
                ;; not used
                ;writeu, lun, maxima_mask
                ;maxima_data[*,rr,pp,ss] = maxima_mask

;;;;;;;;;;;;;;;;;;;;

            endfor
            
            if rr mod 2 eq 0 then begin
                pbar->update, iter/(fdim*sdim) * 100.0
                if (pbar->checkCancel()) then begin
                    pbar->destroy
                    anis = 0
                    obj_destroy, ht
                    return
                endif
            endif

        endfor
        print, "Done Slice", ss
    endfor

    pbar->destroy
    ;; not used
    ;close, lun
    ;free_lun, lun
    
    print, "Time: ", systime(1)-st

    for i = 0, n_elements(dir_ptrs)-1 do begin
        mas_export_nifti, data_ptr=dir_ptrs[i], file_name=destdir+'V_'+string(i, format='(I02)')+'_'+f_append+'.nii'
    endfor
    
    ;; not used
    ;openw, lun, destdir+'HARDI_RECON_DIRS.txt', /get_lun
    ;printf, lun, vlist
    ;close, lun
    ;free_lun, lun
    
    obj_destroy, ht
    ptr_free, dir_ptrs

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_track
;; Created by:
;; Calling Information:   
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;
;; Editing Information:

pro mas_hardi_tractography_track, roiptrs, $
        tess_level=tess_level, $
        exclusion_roi=exclusion_roi, $
        bound_roi=bound_roi, $
        interp_method=interp_method, $
        handle=handle, $
        randomize=randomize, $
        density=density, $
        odf_obj=odf_obj, $
        step_size=step_size, $
        prob_thr=prob_thr, $
        anis_thr=anis_thr, $
        angle_thr=angle_thr, $
        project_dataset=project_dataset, $
        transformation_matrix=transformation_matrix, $
        maximum_split=maximum_split, $
        working_dir=working_dir, $
        preprocessed_directions=preprocessed_directions, $
        track_vector_field=track_vector_field
        
    forward_function mas_hardi_tractography_make_trackvis_header
        
    common scan_data
    
    if (not arg_present(interp_method)) then begin
        interp_method = 0
    endif
    
    if (not arg_present(step_size)) then begin
        step_size = 0.25
    endif
    
    if (not arg_present(prob_thr)) then begin
        prob_thr = 0.5
    endif

    if (not arg_present(tess_level)) then begin
        tess_level = 3
    endif
    
    if (not arg_present(anis_thr)) then begin
        anis_thr = 0.1
    endif
        
    if (n_elements(working_dir) eq 0) then begin
        working_dir = project.current_path
    endif
    
    if (not arg_present(angle_thr)) then begin
        angle_thr_deg = 70.0 ;;cos(70.0*!DTOR)
    endif else begin
        angle_thr_deg = float(angle_thr) ;;cos(angle_thr * !DTOR)
    endelse
    
    if (not keyword_set(maximum_split)) then maximum_split=0

    if (keyword_set(track_vector_field)) then begin
        
        check_1 = total(ptr_valid(*preprocessed_directions))
        check_2 = n_elements(*preprocessed_directions)
        if ((check_1 ne 0 and check_2 ne 0) and (check_1 ne check_2)) then begin
            void = dialog_message('Invalid data detected in preprocessed directions.', /error, /center)
            return
        endif 
        
        data_ptr = preprocessed_directions
        
        tracker = obj_new('MAS_HARDI_TRACTOGRAPHY', odf_obj)
        tracker->setProperty, interpolation_method=interp_method, $
            step_size=step_size, $
            turn_angle_thr=angle_thr_deg, $
            prob_thr=prob_thr, max_branching=maximum_split, $
            anis_halt_thr=anis_thr, diff_data_ptr=data_ptr, data_is_vector_field=1
        odf_obj = obj_new()
        
    endif else begin
        if (not obj_valid(odf_obj)) then begin
            ;;mas_dot_init_obj, radius=15.0, bval_thr=300, lmax=4, /bulk, dot_obj=odf_obj
            mas_odf_init_obj, b_thr=300, lmax=4, /bulk, odf_obj=odf_obj
            ;;mas_mow_init_obj, comp_meth=1, deco_level=2, radius=15.0, /bulk, mow_obj=odf_obj
        endif
        
        tracker = obj_new('MAS_HARDI_TRACTOGRAPHY', odf_obj, tess_level=tess_level)
        tracker->setProperty, interpolation_method=interp_method, $
            step_size=step_size, $
            turn_angle_thr=angle_thr_deg, $
            prob_thr=prob_thr, max_branching=maximum_split, $
            anis_halt_thr=anis_thr
    endelse
              
    nstreams = 0L
    fiber_arr = ptr_new()
    seed_dens = keyword_set(density) ? density : 1
    n_rois = n_elements(roiptrs)

    if (keyword_set(bound_roi) && ptr_valid(bound_roi)) then begin
        bound_roi = obj_new('mas_hardi_ROI', bound_roi, name='ROI_bound')
        tracker->setBoundsROI, bound_roi
    endif else begin
        bound_roi = obj_new()
    endelse

    if (keyword_set(exclusion_roi) && ptr_valid(exclusion_roi)) then begin
        exclusion_roi = obj_new('mas_hardi_ROI', exclusion_roi, name='ROI_exclusion')
        tracker->setExclusionROI, exclusion_roi
    endif else begin
        exclusion_roi = obj_new()
    endelse

    file = dialog_pickfile(title='Save tract file to?', $
                             path=working_dir, filter=['*.trk'], default_extension='trk', $
                             /write, /OVERWRITE_PROMPT)

    if (file eq '') then return

    progressbar = obj_new('progressbar', color='Red', $
                           title='Tractography in Progress...', $
                           text='HARDI Tractography Session')

    progressbar->Start
    print, "Seeding...."
    seedpoints = mas_hardi_tractography_seedROIs(roiptrs, density=seed_dens, randomize=randomize)
    n_seeds = n_elements(*seedpoints)/3
    print, "...Done. Number of Seeds: ", n_seeds

    time_start = systime(1)
    om = obj_new('idlgrmodel')
    
    tracker->printSummary
    
    trackvis_hdr = mas_hardi_tractography_make_trackvis_header()
    openw, lun, file, /get_lun
    writeu, lun, trackvis_hdr
    mat = diag_matrix([[trackvis_hdr.voxel_size],1.0])
    ;; For each seed point in the ROI
    for i = 0L, n_seeds - 1 do begin
    
        sdpoint = reform((*seedpoints)[*,i])
        if (tracker->checkBoundsData(sdpoint)) then continue
        
        stream  = tracker->initializeStreamlines(sdpoint)
        valid = where(obj_valid(stream), n_valid)
        if (n_valid eq 0) then continue
        stream = stream[valid]
        
        ;; For each eligible stream in the seed point
        curr_strm = 0L
        splt_strm = 0
        nstrms = n_elements(stream)
        while (curr_strm lt nstrms) do begin

            status = -1
            iter   =  0
            while (status ne 0) do begin
                status = tracker->jump(stream[curr_strm], split=split)
                if (maximum_split && total(obj_valid(split)) ne 0 && splt_strm lt maximum_split) then begin
                    stream = [stream, split]
                    nstrms++
                    splt_strm++
                    split = 0
                endif
                gl = 0
                iter++
                if (iter gt 10000L) then begin
                    print, "Runaway streamline!"
                    status = 0
                endif
            endwhile
            ;print, "Streams: "+string(n_elements(stream), format='(I0)')
            
            if (not stream[curr_strm]->isExcluded()) then begin
                pts = stream[curr_strm]->getStreamAsPoints()
                ;stream[curr_strm]->getProperty, midpt_idx=midpt, seedpt=seedpt
                nstreams++
                p = n_elements(pts)/3
                pts = vert_t3d(pts, matrix=mat, /no_divide)
                writeu, lun, p
                writeu, lun, pts
            endif
            
            obj_destroy, stream[curr_strm]
            curr_strm++
            
        endwhile ;for ;; stream
                
        obj_destroy, stream[*]
        
        pc_complete = float(i)/float(n_seeds) * 100.0
        if (i mod 10 eq 0) then begin
            progressbar->Update, float(i)/n_seeds  * 100.0
            if (progressbar->CheckCancel()) then begin
                progressbar->Destroy
                goto, CLEANUP
            endif
        endif
        
        if (i mod 50 eq 0) then begin
           heap_gc
        endif
        
    endfor ;; seed point
    
    print, "Total Time: "+string(systime(1)-time_start, format='(F0)')
    
    progressbar->Destroy
    
    CLEANUP:
    ;; rewrite header since now we know what n_count is.
    trackvis_hdr.n_count = nstreams
    point_lun, lun, 0
    writeu, lun, trackvis_hdr
    close, lun
    free_lun, lun

    obj_destroy, tracker
    obj_destroy, odf_obj
    
    return
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_make_trackvis_header
;; Created by: BT
;; Calling Information:   
;;
;;    No arguments
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Creates a trackvis file (.trk) header using the parameters
;;   from the currently selected scan in the scan info window.
;;
;;
;; Editing Information:

function mas_hardi_tractography_make_trackvis_header

    common scan_data

    ci = project.ci
    interp_f = project.procpramarray[ci].freq_interp
    interp_p = project.procpramarray[ci].phase_interp
    interp_s = project.procpramarray[ci].slice_interp
    interp = [ interp_f, interp_p, interp_s ]
    
    fov_f = project.imndarray[ci].f_fov
    fov_p = project.imndarray[ci].p_fov
    fov_s = project.imndarray[ci].s_fov

    voxdims = [project.imndarray[ci].f_voxsz/interp_f,$
               project.imndarray[ci].p_voxsz/interp_p,$
               project.imndarray[ci].s_voxsz/interp_s] * 10.0
    
    hdr = create_struct(name='MAS_HARDI_TRACTOGRAPHY_TRACKVIS_HEADER')
    hdr.id_string[0:4] = byte('TRACK')
    hdr.dim = (size(*project.dataarray[project.ci].state1, /dimensions))[0:2]
    hdr.voxel_size = voxdims
    hdr.origin = [0.0,0,0]
    if (project.procpramarray[ci].have_orientation) then begin
        print, "Using computed or user-defined orientation" ;; Use oriention first.
        sf = float(mas_orient_get_matrix()) ;* (-1)
        print, sf
        vo = byte(mas_orient_from_matrix(sf, /lettercode))
    endif else if (ptr_valid(project.imndarray[project.ci].sform_matrix)) then begin
        print, "Using acquired SFORM matrix."
        sf = float(round(*project.imndarray[project.ci].sform_matrix)) ;; sform if available
        sf = diag_matrix(-1.0/voxdims/interp) # (sf)[0:2,0:2]          ;; scale according to voxel dimensions. (-1) factor due to trackvis's weirdness)
        sf = float(round(sf))                                          ;; round to nearest orthogonal axis
        print, sf
    endif else begin ;; last stand
        print, "No orientation -- using LAS"
        vo = byte('LAS')
    endelse
    print, "Voxel Order: "+strjoin(vo, /single)
    hdr.voxel_order[0:2] = vo

    ;hdr.vox_to_ras = fltarr(4,4)
    ;hdr.vox_to_ras[0:2,0:2] = sf & hdr.vox_to_ras[3,3] = 1.0
    ;print, hdr.vox_to_ras

    hdr.pad2[0:2] = vo
    hdr.image_orientation_patient = [1,0,0,0,1,0]
    hdr.pad1 = byte([32,32])

    hdr.invert_x = 0
    hdr.invert_y = 0
    hdr.invert_z = 0

    hdr.swap_xy = 0
    hdr.swap_yz = 0
    hdr.swap_xz = 0

    hdr.version = 1L
    hdr.hdr_size = 1000L
    return, hdr
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_convert_sav2trackvis
;; Created by:
;; Calling Information:   
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;
;; Editing Information:

pro mas_hardi_tractography_convert_sav2trackvis, infile=infile, outfile=outfile, voxel_dims=voxel_dims

    common scan_data
    
    catch, error_state
    if (error_state) then begin
        catch, /cancel
        void = dialog_message('Operation Failed.', /center, /error)
        return
    endif
    
    if (not keyword_set(infile)) then begin
        infile = dialog_pickfile(title='Select a .sav file to convert.', $
                                 path=project.current_path, /read, filter='*.sav', /must_exist)
        if (infile eq '') then return
    endif
    
    destpath = file_dirname(infile)
    destfile = file_basename(infile, '.sav')+'.trk'
    if (not keyword_set(outfile)) then begin
        outfile = dialog_pickfile(title='Select a destination for the TrackVis file.', $
                                 path=destpath, filter='*.trk', default_extension='trk', /write, $
                                 file=destfile)
        if (outfile eq '') then return
    endif

    restore, infile
    
    hdr = create_struct(name='MAS_HARDI_TRACTOGRAPHY_TRACKVIS_HEADER')
    hdr.id_string[0:4] = byte('TRACK')
    hdr.dim = (*fiber_data).matrix
    hdr.voxel_size = voxel_dims
    hdr.origin = [0.0,0,0]
;    hdr.vox_to_ras = mat
    hdr.voxel_order[0:2] = byte('LPS')
    hdr.pad2[0:2] = byte('LPS')
    hdr.image_orientation_patient = [1,0,0,0,1,0]
    hdr.pad1 = byte([32,32])
    hdr.invert_x = 32
    hdr.invert_y = 32
    hdr.invert_z = 32
    hdr.swap_xy = 32
    hdr.swap_yz = 32
    hdr.swap_xz = 32
    hdr.version = 1L
    hdr.hdr_size = 1000L
    hdr.n_count = (*fiber_data).nfibers
    
    openw, lun, outfile, /get_lun
    writeu, lun, hdr
    mat = float(diag_matrix([[voxel_dims], 1]))
    
    for f = 0L, (*fiber_data).nfibers-1 do begin
        tmp = (*fiber_data).fibers[f]
        fsize = n_elements(*tmp)/3
        writeu, lun, fsize
        *tmp = vert_t3d(*tmp, matrix=mat, /no_divide)
        writeu, lun, *tmp
        ptr_free, tmp
    endfor
    
    close, lun & free_lun, lun
    ptr_free, fiber_data
    heap_gc
    catch, /cancel
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_convert_trackvis2sav
;; Created by:
;; Calling Information:   
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;
;; Editing Information:

pro mas_hardi_tractography_convert_trackvis2sav, infile=infile, outfile=outfile

    common scan_data
    
    ;catch, error_state
    error_state = 0
    if (error_state) then begin
        catch, /cancel
        void = dialog_message('Operation Failed.', /center, /error)
        help, /traceback
        return
    endif
    
    if (not keyword_set(infile)) then begin
        infile = dialog_pickfile(title='Select a .trk file to convert.', $
                                 path=project.current_path, /read, filter='*.trk', /must_exist)
        if (infile eq '') then return
    endif
    
;    destpath = file_dirname(infile)
;    destfile = file_basename(infile, '.trk')+'.sav'
;    if (not keyword_set(outfile)) then begin
;        outfile = dialog_pickfile(title='Select a destination for the .sav file.', $
;                                 path=destpath, filter='*.sav', default_extension='sav', /write, $
;                                 file=destfile)
;        if (outfile eq '') then return
;    endif

    openr, lun, infile, /get_lun    
    hdr = create_struct(name='MAS_HARDI_TRACTOGRAPHY_TRACKVIS_HEADER')
    readu, lun, hdr
    nf = 0L
    fibtmp = fltarr(3)
    fiber_arr = ptrarr(hdr.n_count)
    mid_points = fltarr(3, hdr.n_count)
    mat = float(diag_matrix([[1.0/hdr.voxel_size], 1]))
    
    for f = 0L, hdr.n_count-1 do begin
        readu, lun, nf
        fiber = fltarr(3, nf)
        for i = 0, nf-1 do begin
            readu, lun, fibtmp
            fiber[*,i] = fibtmp
        endfor
        fiber = vert_t3d(fiber, matrix=mat, /no_divide)
        mid_points[*,f] = fiber[*,0]
        fiber_arr[f] = ptr_new(fiber, /no_copy)
    endfor
    
    fdata = { fibers: temporary(fiber_arr), $
        seed_points: temporary(mid_points), $
        transform: diag_matrix(replicate(1D, 4)), $
        branch_levels: [0], $
        matrix: hdr.dim, $
        nfibers: hdr.n_count, $
        filedir: '', $
        tractfile: '', $
        roifile: '' }
    
    fiber_data = ptr_new(fdata, /no_copy)
    view_fibers_hardi, fiber_data, /show_axis, alpha=.5, thick=1, /in_mas, handle=handle, working_dir=project.current_path;, supplemental_omodel=om

    close, lun & free_lun, lun
    heap_gc
    catch, /cancel
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_trackvis_header__define
;; Created by:
;; Calling Information:   
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;
;; Editing Information:

pro mas_hardi_tractography_trackvis_header__define

    struct = { MAS_HARDI_TRACTOGRAPHY_TRACKVIS_HEADER, $
               id_string: bytarr(6), $
               dim: intarr(3), $
               voxel_size: fltarr(3), $
               origin: fltarr(3), $
               n_scalars: 0, $
               scalar_name: bytarr(10,20), $
               n_properties: 0, $
               property_name: bytarr(10,20), $
               vox_to_ras:fltarr(4,4), $
               reserved: bytarr(444), $
               voxel_order: bytarr(4), $
               pad2: bytarr(4), $
               image_orientation_patient: fltarr(6), $
               pad1: bytarr(2), $
               invert_x: byte(0), $
               invert_y: byte(0), $
               invert_z: byte(0), $
               swap_xy: byte(0), $
               swap_yz: byte(0), $
               swap_xz: byte(0), $
               n_count: 0L, $
               version: 1L, $
               hdr_size: 1000L }

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_probtrack
;; Created by: BT
;; Calling Information:
;;
;; Not called by MAS explicitly. This is an experimental processing method
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; 
;;  Given some data files that contain a distribution of directions, this
;;  procedure performs probabilistic fiber tracking. 
;;
;;  It is currently in development and its usage has not been finalized.
;;
;;
;; Editing Information:

pro mas_hardi_tractography_probtrack, $
        roiptrs, $                        ;; array of pointers containing seed rois
        excl_roi=excl_roi, $              ;; array of pointers containing exclusion rois
        wayp_roi=wayp_roi, $              ;; array of pointers containing waypoint rois
        bound_roi=bound_roi, $            ;;   "    "    "        "       boundary rois
        handle=handle, $                  ;;
        randomize=randomize, $            ;; set this keyword to randomize seed point placement within a voxel
        density=density, $                ;; set this to the number of seed points per voxel to create
        step_size=step_size, $            ;; step size in fractions of a voxel per tracking step
        anis_thr=anis_thr, $              ;; fa threshold at which to stop tracking
        angle_thr=angle_thr, $
        project_dataset=project_dataset, $
        interp_method=interp_method, $    ;; interpolation method
        theta_ptrs=theta_ptrs, $          ;; array of pointers containing the theta distribution(s)
        phi_ptrs=phi_ptrs                 ;; array of pointers contining the phi distribution(s)

    common scan_data
    
    if (not arg_present(interp_method)) then begin
        interp_method = 0
    endif
    
    if (not arg_present(step_size)) then begin
        step_size = 0.25
    endif

    if (not arg_present(angle_thr)) then begin
        angle_thr_rad = cos(70.0*!DTOR)
    endif else begin
        angle_thr_rad = cos(angle_thr * !DTOR)
    endelse
    
    tracker = obj_new('MAS_HARDI_TRACTOGRAPHY_PROBTRACK', 2, $
               theta_ptrs, phi_ptrs)
    tracker->setProperty, interpolation_method=interp_method, $
                          step_size=step_size, $
                          turn_angle_thr=angle_thr_rad
                                                        
    nstreams = 0L
    fiber_arr = ptr_new()
    seed_dens = keyword_set(density) ? density : 1
    n_rois = n_elements(roiptrs)

    track_mask = bytarr((size(*theta_ptrs[0], /dimensions))[0:2])
    master_track_mask = lonarr(size(track_mask, /dimensions))

    if (keyword_set(bound_roi) && ptr_valid(bound_roi)) then begin
        bound_roi = obj_new('mas_hardi_ROI', bound_roi, name='ROI_bound')
        tracker->setBoundsROI, bound_roi
    endif else begin
        bound_roi = obj_new()
    endelse
    
    if (keyword_set(excl_roi) && ptr_valid(excl_roi)) then begin
        tracker->setExclusionRoi, excl_roi[0]
    endif
    
    if (keyword_set(wayp_roi) && n_elements(wayp_roi) ne 0) then begin
        valid = ptr_Valid(wayp_roi)
        for v = 0, n_elements(valid)-1 do begin
            tracker->addWaypointMask, wayp_roi[v]
        endfor
    endif
    
    seedpoints = mas_hardi_tractography_seedROIs(roiptrs, density=seed_dens, randomize=randomize)
    n_seeds = n_elements(*seedpoints)/3
    print, "# Seeds: ", n_seeds;; & wait, 2

    file = dialog_pickfile(title='Save tract file to?', $
                           path=project.current_path, $
                           /write, /OVERWRITE_PROMPT)

    progressbar = obj_new('progressbar', color='Red', $
                           title='Tractography in Progress...', $
                           text='HARDI Tractography Session')
    progressbar->Start
    
    waytotal = 0L
    ;; For each seed point in the ROI
    for i = 0L, n_seeds - 1 do begin
    
        sdpoint = reform((*seedpoints)[*,i])
        stream  = tracker->initializeStreamlines(sdpoint)
        valid = where(obj_valid(stream), n_valid)
        if (n_valid eq 0) then continue
        stream = stream[valid]
        
        ;; For each eligible stream in the seed point
        for s = 0L, n_elements(stream)-1 do begin
        
            status = -1
            iter   =  0
            tracker->resetCounters
            while (status ne 0) do begin
                status = tracker->jump(stream[s])
                iter++
            endwhile
            
            ;print, string(iter, format='(I0)')+ ' jumps'
            
            tracker->getHitCounts, madewaypoints=mw
            if (mw eq 0) then continue
            
            ;; need for view_fibers
            track_mask[*] = 0
            pts = stream[s]->getStreamAsPoints()
            for p = 0L, n_elements(pts)/3 -1 do begin
                track_mask[ pts[0,p], pts[1,p], pts[2,p] ] = 1
            endfor
            
            master_track_mask += track_mask
            waytotal++
            
        endfor ;; stream
        
        obj_destroy, stream

        ;pc_complete = float(i)/float(n_seeds) * 100.0
        progressbar->Update, float(i)/n_seeds  * 100.0
        if (progressbar->CheckCancel()) then begin
            progressbar->Destroy
            goto, CLEANUP
        endif
        
    endfor ;; seed point
    
    progressbar->Destroy

    print, "WAYTOTAL: ", waytotal
    if (file ne '') then begin
        mtrack = ptr_new(master_track_mask, /no_copy)
        mas_export_nifti, data_ptr=mtrack, file_name=file
        ptr_free, mtrack
    endif

    CLEANUP:
    obj_destroy, tracker
    obj_destroy, stream
    ;;obj_destroy, streams

    return
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_hardi_tractography_addROI
;; Created by: BT
;; Calling Information:   
;;
;;  Called by the view_fibers_hardi procedure.
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;  Used to add an ROI into the viewer. 
;;
;; Editing Information:

pro mas_hardi_tractography_addROI, state, roimask_ptrs

    nrois = 0
    
    for r = 0, n_elements(roimask_ptrs)-1 do begin
        oroi = obj_new('mas_hardi_ROI', roimask_ptrs[r], name="ROI_N")
    
        subs =  oroi->getSubRegions()
    
        for i = 0, n_elements(subs)-1 do begin
            tmp_seeds = oroi->makeSeeds(subs[i], density=1)
            if (tmp_seeds[0] eq -1) then continue
            roi_arr = (nrois++ eq 0) ? ptr_new(tmp_seeds) : [ roi_arr, ptr_new(tmp_seeds) ]
            endfor
    
        obj_destroy, oroi
    
    endfor

    VFH_add_roi, state, roi_arr
    
end

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
    
    ht = obj_new('mas_hardi_tractography', method_spec, data_ptr=data, $
                 prob_thr=0.5, tess_level=4)

    data_sz = size((*data), /dimensions)
    
    fdim = data_sz[0]
    pdim = data_sz[1]
    sdim = data_sz[2]
    adim = data_sz[3]
       
    ;mas_tessellator_make_icos, level=5, vertlist=tmp;;, polylist=polylist
    ;vlist = tmp[*,where(tmp[2,*] ge 0)]
    
    ;method_spec->setReconstructionDirections, vlist
    ;tx = *project.imndarray[project.ci].acq_matrix
    ;if (obj_class(method_spec) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
    ;    vlist = vert_t3d(vlist, matrix=tx, /double)
    ;endif
    
    ;; prepare the "restart vertices" for finding max Pr
    ;mas_tessellator_make_icos, level=0, vertlist=tmp
    ;vlist_restart = tmp[*,where(tmp[2,*] ge 0)]
    
    ;if (obj_class(method_spec) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
    ;    vlist_restart = vert_t3d(vlist_restart, matrix=tx, /double)
    ;endif
    
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
                if (total(abs(Pr)) eq 0) then continue
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
                        ht->getMaxima, Pr, vs=vs, ps=ps, nvs=nvs
                        ;mas_hardi_tractography_get_maxima, Pr, vlist_restart, 0.65, $
                        ;    vs=vs, ps=ps, nvs=nvs
                            
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

pro mas_hardi_tractography_get_maxima, Pr_in, vlist_restart, prob_thr, $
                                 vs=vs, ps=ps, nvs=nvs, maxima_mask=maxima_mask

    common probtrack_objective, Pr, vlist, pr_interp_tol

    nPr = n_elements(Pr_in)
    n_rsverts = n_elements(vlist_restart)/3
    pr_interp_tol = 0.25
    
    vs = fltarr(3, n_rsverts)
    ps = fltarr(n_rsverts)
    maxima_mask = bytarr(n_elements(vlist)/3)
    
    for rv = 0, n_rsverts-1 do begin
       
       ;; Select a vertex from the "restart" vertices; call it "test"
       test = reform((vlist_restart)[*,rv])
       
       ;; find neighbors of "test" (within pr_interp_tol)
       ;; call these the "potential coordinates"
       pot_coords = where(abs((vlist)[0,*] - test[0]) lt pr_interp_tol and $
                          abs((vlist)[1,*] - test[1]) lt pr_interp_tol and $
                          abs((vlist)[2,*] - test[2]) lt pr_interp_tol, count)
       
       ;; proceed in the direction of greatest increase in Pr
       ;; assumes that the probability function is reasonably smooth
       ;; (which is it in this case)
       curr_pr = 0.0
       iter = 0L
       best_coord_index = -1L
       while (1) do begin
          if (count eq 0) then begin
             print, "No vertices returned"
             break
          endif else begin
             
             ;; get the probabilites associated with the pot_coords
             ;; and choose the max out of all of them
             possible_pr = Pr_in[pot_coords]
             max_pr = max(possible_pr)
             possible_coord = (pot_coords[where(possible_pr eq max_pr)])[0]
             
             ;; if this probabilty is greater than current max, save
             ;; it + the coords
             if (max_pr gt curr_pr) then begin
                test = (vlist)[*, possible_coord]
                curr_pr = max_pr
                best_coord_index = possible_coord
                iter++
             endif else begin
                ;; break out of the while() loop when the probabilites begin
                ;; to decrease
                iter++
                break
             endelse
          endelse

          ;; now find the neighbors for the next iteration
          pot_coords = where(abs((vlist)[0,*] - test[0]) lt pr_interp_tol and $
                             abs((vlist)[1,*] - test[1]) lt pr_interp_tol and $
                             abs((vlist)[2,*] - test[2]) lt pr_interp_tol, count)
          
       endwhile

       ;; at the end of the while loop, "test" will have the vector
       ;; associated with the maximum probability 

       ;; correct for anipodal symmetry 
       vs[*,rv] = total(test * [0,0,1]) lt 0 ? -test : test
       ps[rv] = max_pr
       if (best_coord_index ne -1) then maxima_mask[best_coord_index] = 1B
       
    endfor
    
    ;; weed out small probabilities; it is assumes that there will
    ;; be at least one
    eligible = where(ps gt prob_thr)
    vs = vs[*,eligible]
    ps = ps[eligible]
    
    dps = ((vs ## transpose(vs)))[0,*]
    
    dps_sorted = sort(dps)    
    dps = dps[dps_sorted]
    vs  = vs[*,dps_sorted]
    ps  = ps[dps_sorted]
    
    dps_uniq = uniq(dps)
    vs = vs[*,dps_uniq]
    ps = ps[dps_uniq]

    ps_sorted = reverse(sort(ps))
    ps = ps[ps_sorted]
    vs = vs[*,ps_sorted]
    
    nvs = n_elements(dps_uniq)

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
pro mas_hardi_tractography_show_voxel, reco_obj, coord, $
        prob_thr=prob_thr, $
        thru_angle_lmt=thru_angle_lmt, $
        branch_angle_lmt=branch_angle_lmt, $
        pr_int_tol=pr_int_tol, $
        powell_tol=powell_tol, $
        recon_tess_level=recon_tess_level, $
        restart_tess_level=restart_tess_level, $
        nbootstrap=nbootstrap, $
        use_state_ptr=use_state_ptr, compare_adt=compare_adt
        
    common scan_data
    common probtrack_objective, Pr, vlist, pr_interp_tol
    
    ci = project.ci
    pr_interp_tol    = keyword_set(pr_int_tol)      ? float(pr_int_tol)     : 0.085
    thru_ang_lmt     = keyword_set(thru_ang_lmt)    ? float(thru_ang_lmt)   : 45.0
    branch_ang_lmt   = keyword_set(branch_ang_lmt)  ? float(branch_ang_lmt) : 45.0
    recon_tess_level = keyword_set(recon_tess_level) ? fix(recon_tess_level) : 4
    prob_thr         = keyword_set(prob_thr)        ? float(prob_thr)       : 0.01
    prob_branch_tol= 0.5
        
    mas_tessellator_make_icos, level=(1 > recon_tess_level)<5, $
        vertlist=vlist, $
        polylist=polylist
        
    reco_obj->setReconstructionDirections, vlist
    
    restart_tess = (keyword_set(restart_tess_level)) ? floor(restart_tess_level) : 0
    
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
    mas_tessellator_make_icos, level=restart_tess_level, vertlist=vlist_restart
    
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
        
    nbootstrap = keyword_set(nbootstrap) ? fix(nbootstrap) : 0
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if (nbootstrap gt 0) then begin
        seed = systime(1)
        can_handle_bootstrapping = reco_obj->makeBootstrap(data=data)
        if (can_handle_bootstrapping eq 0) then begin
            void = dialog_message(["This reconstruction method does not support", $
                                   "statistical bootstrapping. Please use another method."], $
                                   /center, /error)
            return
        endif
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
            if (ps[j] lt prob_thr) then continue
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
    
    if (keyword_set(compare_adt)) then begin
        tmp = vlist
        evec = reform((*project.dataarray[project.ci].eign_vec)[coord[0], coord[1], coord[2], * , *])
        eval = diag_matrix(reform((*project.dataarray[project.ci].eign_val)[coord[0], coord[1], coord[2], *]))
        eval = eval/max(eval) & print, eval
        mat = diag_matrix(replicate(1.0,4))
        mat[0:2,0:2] = evec # eval # transpose(evec)
        tmp = vert_t3d(tmp, matrix=mat, /no_copy) 
        adtsurface = obj_new('idlgrpolygon', tmp, $
                name=strcompress("Anis=", /remove_all),$
                polygons=polylist, alpha_channel=0.65, $
                color=[180,180,180], style=1, shading=0)
        oadtdirs = objarr(3)
        for j = 0, 3-1 do begin ;;;(nthru+ncross)-1 do begin
            d = reform(evec[*,j])
            if j eq 0 then col = [255,0,0] else col = [0,0,255]
            oadtdirs[j] = obj_new('idlgrpolyline', $
                [ [-(d*step_size*0.65)], [(d*step_size*0.65)] ], $
                color = col, $
                thick=1)
        endfor
                
        xobjview, [adtsurface, oadtdirs], background=[0,0,0], renderer=0
        
    endif
    
    HEAP_GC
    
end
