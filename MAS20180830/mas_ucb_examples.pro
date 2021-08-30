;+
; :Description:
;    Displays an inverted version of the in-plane image of the currently selected slice.
;
; :Params:
;    coord -- the coordinate of the crosshairs in the orthogonal viewer
;    state -- the state structure pertaining to the orthoginal viewer.
;
;
;
; :Author: btt
;-
pro UCB_ex_invert, coord, state

    common scan_data, project

    ;; Get the project index -- the scan id (in mas) of the data currently
    ;; being shown in the orthogonal viewer. This index is used to get
    ;; data and scan parameters from MAS internals.
    proj_index = (*state).proj_index

    ;; Check to see if a window already exists (we are using window id 5,
    ;; but any window id will do (up to 32).
    device, window_state=w    
    if (w[5] eq 0) then begin
        ;; Here we set the window size to the x/y dims of the current scan
        window, 5, xsize=project.imndarray[proj_index].fdim, $
                   ysize=project.imndarray[proj_index].pdim
    endif
    
    ;; extract the 2D slice referred to be coord[2] from the
    ;; data array in MAS.
    temp = (*project.dataarray[proj_index].state1)[*,*, coord[2], 0]
    
    ;; Find the max value over the 2D slice
    max_temp = max(temp)
  
    ;; Invert the image and display in our window. TVSCL automatically
    ;; scales the intensity to byte range (0-255)
    tvscl, (max_temp - temp)
    
end

;+
; :Description:
;    Example of how to compute DTI on-the-fly for a single voxel.
;
; :Params:
;    coord -- the coordinate of the crosshairs in the orthogonal viewer
;    state -- the state structure pertaining to the orthoginal viewer.
;
;
;
; :Author: btt
;-
pro UCB_ex_dti, coord, state

    common scan_data, project
    
    proj_index = (*state).proj_index
    
    ;; if there isn't a valid b_matrix, then we should probably not try to do DTI
    if (not ptr_valid(project.imndarray[proj_index].b_matrix)) then return

    do_adt_regress_standard_singlevoxel, coord, fa=fa, eval=eval, evec=evec, ad=ad, fit=fit
    
    print, 'FA.....: '+strtrim(fa,2)
    print, 'AD.....: '+strtrim(ad,2)
    print, 'S0.....: '+strtrim(fit[0],2)
    print, 'EVALS..: '+strjoin(strtrim(eval,2), ', ')
    print, 'EVEC1..: '+string(evec[*,0], format='(%"[ %+0.4f, %+0.4f, %+0.4f ]")')
    print, 'EVEC2..: '+string(evec[*,1], format='(%"[ %+0.4f, %+0.4f, %+0.4f ]")')
    print, 'EVEC3..: '+string(evec[*,2], format='(%"[ %+0.4f, %+0.4f, %+0.4f ]")')
    print, 'TENSOR.: '+strjoin(strtrim(fit[1:*], 2), ', ')

end

;+
; :Description:
;    Computes the maximum intensity projection for each slicing plane
;
; :Params:
;    coord -- the coordinate of the crosshairs in the orthogonal viewer
;    state -- the state structure pertaining to the orthoginal viewer.
;
; :Author: btt
;-
pro UCB_ex_mip, coord, state

    common scan_data, project
    common UCB_ex_mip_data, last_coord
    
    proj_index = (*state).proj_index
    
    ;; grab the pointer so we don't have to type all this. 
    d = project.dataarray[proj_index].state1

    dims = size(*d, /dimensions)
    
    ;; doesn't work for single slice data.
    if (dims[0] eq 1 or dims[1] eq 1 or dims[2] eq 1) then return
    
    device, window_state=w    
    if (w[5] eq 0) then begin
        window, 5, xsize=2*dims[0]+dims[1], $
                   ysize=max(dims[1:2])
    endif
    
    ;; Quick way to find out which plane to update so
    ;; all there don't have to be updated every time.
    if (n_elements(last_coord) eq 0) then begin
        diff = [1,1,1]
    endif else begin
        diff = abs(last_coord - coord)
    endelse

    last_coord = coord

    ;; compute the maximum intensity proj.
    if (diff[2] gt 0.5) then begin
        mipZ = max((*d)[*,*,0:(coord[2]-1)>1,0], dimension=3)
        tvscl, mipZ, 0,0
    endif
    
    if (diff[0] gt 0.5) then begin
        mipX = max((*d)[0:(coord[0]-1)>1,*,*,0], dimension=1)
        tvscl, mipX, dims[0], 0
    endif
    
    if (diff[1] gt 0.5) then begin
        mipY = max((*d)[*,0:(coord[1]-1)>1,*,0], dimension=2)
        tvscl, mipY, dims[0]+dims[1],0
    endif
    
end
