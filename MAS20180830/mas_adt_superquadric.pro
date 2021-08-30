;; $Id$
;;

; Subroutine name: sp_make_sphere
; Created by: BT, 2008-06
; Calling Information:
;
; vdim  - IN - the pixel radius
; vert  - OUT - the polygon vertices suitable for IDLgrPolygon or
;         POLYSHADE
; poly  - OUT - polygon adjency list for above vertices

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; creates a polygon representation of a sphere

; Editing Information:
function sp_make_sphere, vdim, vert=vert, poly=poly

    volume = dblarr(vdim,vdim,vdim)

    for x=0, vdim-1 do begin
        for y=0, vdim-1 do begin
            for z=0, vdim-1 do begin
                volume[x,y,z] = -total( (vdim/2.-.5-[x,y,z])^2 )
            endfor
        endfor
    endfor
    
    scale3, xrange=[0, vdim], yrange=[0, vdim], zrange=[0, vdim]    
    st=-(vdim/4.)^2*3.25
    shade_volume, volume, st, vert, poly, low=1

    return, (n_elements(vert) eq 0) ? 0 : 1

end

; Subroutine name: sp_superquadric
; Created by: BT, 2008-06
; Calling Information:
;  dim  - IN - volume of enclosing cube in pixels
;  evals - IN - the eigenvalues. used to compute cp and c;
;  gamma - IN - the paramater that controls the 'squareness'
;  vert  - OUT - the polygon vertices suitable for IDLgrPolygon or
;          POLYSHADE
;  poly  - OUT - polygon adjency list for above vertices

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Creates a polygon represtentation of a superquadric glyph scaled
; according to the eigenvalues of a rank-2 tensor
; (cf. http://www.cs.utah.edu/~gk/papers/vissym04/)
; all variable names are taken from the paper.

; Editing Information:

function sp_superquadric, dim, evals, gamma, vert=vert, poly=poly

    volume = dblarr(dim,dim,dim)

    vert = 0
    poly = 0
    gamma = float(gamma)

    coord_xyz = abs(findgen(dim) - 0.5*dim)/(0.4*dim)

    cl = (evals[0]-evals[1])/(total(evals))
    cp = (2*(evals[1]-evals[2]))/(total(evals))

    if (cl ge cp) then begin

        alpha = (1.0 - cp)^gamma
        beta  = (1.0 - cl)^gamma
        
    endif else begin

        alpha = (1.0 - cl)^gamma
        beta  = (1.0 - cp)^gamma

    endelse

    alpha_beta = alpha/beta

    c2a = coord_xyz^(2.0/alpha)
    c2b = coord_xyz^(2.0/beta)
    
    for x = 1, dim-2 do begin
        for y = 1, dim-2 do begin
            for z = 1, dim-2 do begin
                
                test = (cl ge cp) ? $
                       (c2a[y] + c2a[z])^(alpha_beta) + c2b[x] : $
                       (c2a[y] + c2a[x])^(alpha_beta) + c2b[z]
                
                if (test - 1.0 le 1e-2) then volume[x,y,z] = 1   

            endfor
        endfor
    endfor
                
    shade_volume, volume, 0, vert, poly
    
    return, (n_elements(vert) eq 0) ? 0 : 1
    
end

; Subroutine name: mas_adt_glyph_configure_event
; Created by: BT, 2008-06
; Calling Information:
;
; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; event handler for configure gui

; Editing Information:
pro mas_adt_sp_configure_event, event

    common scan_data
    ci = project.ci

    adt_sp_state = project.procpramarray[ci].adt_sp_state
    difftools_state = project.procpramarray[ci].difftools_state

    if (not ptr_valid(adt_sp_state)) then return

    uname = widget_info(event.id, /uname)

    widget_control, event.id, get_value=value
    
    case uname of 

        'glyph_edge': (*adt_sp_state).edge = value
        'fa_thr':  (*adt_sp_state).fa_thr = value
        'eval_divisor': (*adt_sp_state).eval_divisor = value
        'btn_create': begin

           if ((*difftools_state).disp_3d) then begin
                
                case project.procpramarray[ci].slice_axis of
                    0: slice = project.procpramarray[ci].sdim_start
                    1: slice = project.procpramarray[ci].pdim_start
                    2: slice = project.procpramarray[ci].fdim_start
                endcase
                
                mas_adt_3d, [ slice ]
                
            endif else begin

                mas_adt_superquadric, /display, /image_only

            endelse

        end
        'btn_batch': begin
            slices = mas_get_multi_slice_selection(event.top, count=count)
            if (count eq 0) then return

            if ((*difftools_state).disp_3d) then begin
                mas_adt_3d, [ slices ]
                
            endif else begin
 
                for sl = 0, n_elements(slices)-1 do begin
                    mas_adt_superquadric, $
                      slice=slices[sl], /display, /image_only
                endfor
                
            endelse
        end

        else: print, 'mas_adt_glyph_configure_event: unknown event: '+uname+', value:'+value

    endcase

end

; Subroutine name: mas_adt_glyph_configure_event
; Created by: BT, 2008-06
; Calling Information:
;
; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; makes a state pointer with default parameters

; Editing Information:

function mas_adt_sp_make_default_state

    return, ptr_new({ ADT_STATE, $
                      fa_thr: 0.01, $
                      eval_divisor: 1, $
                      edge: 0 $
                    }, /no_copy)
end

; Subroutine name: mas_adt_sp_redraw
; Created by: BT, 2008-06
; Calling Information:
;
; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; redraws the ADT tab within diffusion tools window

; Editing Information:

pro mas_adt_sp_redraw, adt_base

    common scan_data
    ci = project.ci
    
    adt_sp_state = project.procpramarray[ci].adt_sp_state
    ;; create if not exists
    if (not ptr_valid(adt_sp_state)) then begin
        
        if (project.ni eq 0) then begin

            sen = 0
            widget_control, adt_base, get_uvalue=widget_state
            widget_control, (*widget_state).sl_fa_thr           , sensitive=sen
            widget_control, (*widget_state).sl_glyph_edge       , sensitive=sen
            return

        endif else begin

            adt_sp_state = mas_adt_sp_make_default_state()
            project.procpramarray[ci].adt_sp_state = adt_sp_state

        endelse

    endif

    sen = 1
    widget_control, adt_base, get_uvalue=widget_state
    widget_control, (*widget_state).sl_fa_thr      , sensitive=sen, set_value=(*adt_sp_state).fa_thr
    widget_control, (*widget_state).sl_glyph_edge  , sensitive=sen, set_value=(*adt_sp_state).edge
    widget_control, (*widget_state).sl_eval_divisor, sensitive=sen, set_value=(*adt_sp_state).eval_divisor
    
end

; Subroutine name: mas_adt_sp_gui_make
; Created by: BT, 2008-06
; Calling Information:
;
; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; creates the gui holding options for this particular visualization

; Editing Information:

pro mas_adt_sp_gui_make, base, state

    common scan_data
    common common_widgets

    ci = project.ci

    if (ptr_valid(project.procpramarray[ci].adt_sp_state)) then begin
        adt_sp_state = project.procpramarray[ci].adt_sp_state
    endif else begin
        adt_sp_state = mas_adt_sp_make_default_state() 
        project.procpramarray[ci].adt_sp_state = adt_sp_state
    endelse
    
    param_base = widget_base(base, column=2, /grid, xpad=0)

    sl_glyph_edge  = widget_slider(param_base, $
                                minimum = 0, $
                                maximum = 50, $
                                scroll  = 2, $
                                value   = (*adt_sp_state).edge, $
                                uname   = 'glyph_edge', $
                                event_pro='mas_adt_sp_configure_event', $
                                title   = "Squareness")

    sl_fa_thr = cw_fslider(param_base, uname='fa_thr', $
                        minimum=0.01, $
                        maximum=1.0, $
                        scroll=0.1, $
                        value=(*adt_sp_state).fa_thr, $
                        /edit, $
                        title='FA Threshold')

    sl_eval_divisor  = widget_slider(param_base, $
                                     minimum = 1, $
                                     maximum = 10, $
                                     scroll  = 1, $
                                     value   = ((*adt_sp_state).eval_divisor) > 1, $
                                     uname   = 'eval_divisor', $
                                     event_pro='mas_adt_sp_configure_event', $
                                     title   = "E-value Divisor")

    widget_state = ptr_new({ ADT_WIDGET_STATE, $
                             sl_glyph_edge: sl_glyph_edge, $
                             sl_eval_divisor: sl_eval_divisor, $
                             sl_fa_thr: sl_fa_thr }, /no_copy)
                            
    widget_control, base, set_uvalue=widget_state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mas_adt_3d, slice_range, min_fa=min_fa, max_fa=max_fa

    common scan_data
    ci = project.ci

    mas_load_state_1

;    catch, error_status
;    if (error_status ne 0) then begin
;        print, "Error!"
;    endif

    if (not ptr_valid(project.imndArray[ci].bval_Array)) then begin
        void = dialog_message(['B values are not present for this scan.', $
                               'Please make sure that this is a diffusion-',$
                               'weighted scan.'], /error, /center)
        return
    endif

    if n_elements(vdim) eq 0  then vdim = 10
    if n_elements(vmult) eq 0 then vmult = 1.0
 
    if not keyword_set(min_fa) then min_fa = 0.0
    if not keyword_set(max_fa) then max_fa = 1.0

    if (ptr_valid(project.procpramarray[ci].adt_sp_state)) then begin
        adt_state = project.procpramarray[ci].adt_sp_state
    endif else begin
        adt_state = mas_adt_sp_make_default_state()
    endelse

    if (ptr_valid(project.procpramarray[ci].difftools_state)) then begin
        difftools_state = project.procpramarray[ci].difftools_state
    endif else begin
        difftools_state = mdt_make_default_state()
    endelse

    fa_thr       = (*adt_state).fa_thr
    gamma        = (*adt_state).edge/5.0
    crop_roi     = (*difftools_state).crop_roi
    mask_roi     = (*difftools_state).mask_roi
    vdim         = (*difftools_state).glyph_size
    vmult        = 1;(*difftools_state).glyph_size_mult
    l0_intensity = (*difftools_state).l0_intensity / 100.0
    l1_intensity = (*difftools_state).l1_intensity / 100.0
    l2_intensity = (*difftools_state).l2_intensity / 100.0
    image_opac   = (*difftools_state).image_opacity   / 100.0
    glyph_opacity= (*difftools_state).glyph_opacity   / 100.0
    glyph_shiny  = (*difftools_state).glyph_shiny
    use_adt_color= (*difftools_state).use_adt_color
    ;; incorporate freq, phase, slice viewing direction
    sz_img = size(*project.dataarray[ci].state1, /dimensions)
    adim = project.imndarray[ci].adim

    ;; Set up dimensions and grab the slice data
    slice_axis = project.procPramArray[project.ci].slice_axis
    case slice_axis of
        0: begin                ; freq_phase
            xdim = sz_img[0]
            ydim = sz_img[1]
            zdim = sz_img[2]
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].p_fov
        end
        
        1: begin                ; freq_slice
            xdim = sz_img[0]
            ydim = sz_img[2]
            zdim = sz_img[1]
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].s_fov
        end
        
        2: begin                ; phase_slice
            xdim = sz_img[1]
            ydim = sz_img[2]
            zdim = sz_img[0]
            xfov = project.imndarray[ci].p_fov
            yfov = project.imndarray[ci].s_fov
        end
    endcase
    
    xstart = 0
    ystart = 0

    glyph_color    = [50,50,100] 
    glyph_diffuse  = [90,90,110]
    glyph_specular = [150,150,170]

    vdim_effective = vmult*vdim
    mas_adt_sp_init, vdim_effective, adt_spec=adt_spec

    if gamma eq 0.0 then surf = sp_make_sphere(vdim_effective, vert=ell_verts, poly=ell_polys)
    osurfs = objarr(xdim, ydim, n_elements(slice_range))

    progressBar = obj_new('progressbar', title='Computing ADT', text='Computing ADT...', /fast_loop)
    progressBar->start
    iter = 0
    cols = indgen(255)
    junk = get_orient_rgb_axis(tx_matrix=sa_tx)

    syst = systime(1)

    mask = mas_roi_get_current_mask([xdim, ydim], crop=crop_dims, /no_transform)

    for sl = 0, n_elements(slice_range)-1 do begin
        
        glyph_color=abs([255-8*cols[2*sl], sl*20, 10*cols[3*((sl mod 2)+1)]])
        slice = slice_range[sl]
        ;;print, 'slice: ', slice

        for x = 0, xdim-1 do begin

            for y = 0, ydim-1 do begin

                nevermind = 0

                if (ptr_valid(mask)) then begin
                    if ((*mask)[x,y] eq 0) then continue
                endif

                ;; get the eigvenvalues
                case slice_axis of
                    0: begin
                        if (*project.dataarray[ci].state1)[x,y,slice,0] eq 0 then nevermind = 1
                        evals = reform((*project.dataarray[ci].eign_val)[x,y,slice,*])
                        fa    = reform((*project.dataarray[ci].frac_ani)[x,y,slice])
                        evecs = transpose(reform((*project.dataarray[ci].eign_vec)[x,y,slice,*,*]))
                    end
                    1: begin
                        if (*project.dataarray[ci].state1)[x,slice,y,0] eq 0 then nevermind = 1
                        evals = reform((*project.dataarray[ci].eign_val)[x,slice,y,*])
                        fa    = reform((*project.dataarray[ci].frac_ani)[x,slice,y])
                        evecs = transpose(reform((*project.dataarray[ci].eign_vec)[x,slice,y,*,*]))
                    end
                    2: begin
                        if (*project.dataarray[ci].state1)[slice,x,y,0] eq 0 then nevermind = 1
                        evals = reform((*project.dataarray[ci].eign_val)[slice,x,y,*])
                        fa    = reform((*project.dataarray[ci].frac_ani)[slice,x,y])
                        evecs = transpose(reform((*project.dataarray[ci].eign_vec)[slice,x,y,*,*]))
                    end
                endcase
                
                
                if (total(evals) eq 0 or total(finite(evals, /nan)) gt 0 or nevermind eq 1) then continue
                
                ;; normalize so that eval[0] = 1
                evals = abs(evals/max(evals))
                
                if (fa lt fa_thr) then begin 
                    evals = evals * 0.1 
                    continue 
                endif
                
                if (max(evals) ne evals[0]) then begin
                    ;print, "Unordered eigenvalue error at ("+strcompress(string(x)+','+string(y), /remove_all)+') ', evals
                    continue
                endif
                
                ;; Eqn (N) from the paper, evec matrix use as
                ;; transformation matrix
                evecs = [[evecs[*,0], 0], $
                         [evecs[*,1], 0], $
                         [evecs[*,2], 0], $
                         [0, 0, 0, 1] ]
                
                ;; get the surface
                if (gamma ne 0.0) then begin
                    surf = sp_superquadric(vdim, evals, gamma, poly=polylist, vert=vertlist)
                    if (surf eq 0) then continue
                endif else begin
                    polylist = ell_polys
                    vertlist = ell_verts
                endelse
                
                if (total(finite(vertList, /nan)) ne 0) then continue
                
                if iter mod 20 eq 0 then progressBar->Update, float(iter)/(xdim*ydim)*100.0
                
                ;; transformation matrix: evecs (rotation) * evals (scale)
                tx = evecs ## diag_matrix([evals, 1.0]) # sa_tx
                
                ;; before applying the trans. matrix, we'd like to have
                ;; the center of the surface at the origin. So we make a
                ;; shift array
                sz = (size(vertlist))[2]
                shift = transpose([ [replicate(-(vdim/2.), sz)], $
                                    [replicate(-(vdim/2.), sz)], $
                                    [replicate(-(vdim/2.), sz)] ])
                ;; shift, rotate, de-shift
                vertList += shift
                vertList = vert_t3d(vertList, /no_copy, matrix=tx)
                vertList -= shift
                
                ;; this will translate the vertex coordinates to their
                ;; position relative to the background image
                shift = transpose([ [replicate(x*vdim, sz)], $
                                    [replicate(y*vdim, sz)], $
                                    [replicate(slice*vdim, sz)] ])
                vertList += shift
                
                ;; compute the glyph color based on FA * principal evector
                if (keyword_set(use_adt_color) and project.procpramarray[ci].adt_proccess_flag) then begin
                    glyph_color = adt_get_orientation_color(x, y, slice) 
                endif
                
                osurf = obj_new('idlgrpolygon', vertList, polygons=polylist, $
                                name      = strcompress('Glyph: ('+string(x)+','+string(y)+')'), $
                                specular  = glyph_specular, $
                                diffuse   = glyph_diffuse, $ ;glyph_color, $
                                ambient   = glyph_color, $
                                shading   = glyph_shading, $
                                color     = glyph_color, $
                                shininess = glyph_shiny,$
                                alpha_channel = glyph_opacity, $
                                reject    = 0, $
                                texture_interp  = 1) 
                
                osurfs[x,y,sl] = osurf
                
            endfor
            
            progressbar->update, (100.0) * x/xdim
            if (progressBar->checkcancel()) then begin
                progressBar->destroy
                if (ptr_valid(mask)) then ptr_free, mask
                obj_destroy, osurfs
                return
            endif
            
        endfor
        
    endfor
    
    progressBar->destroy
    
    print, '... Done. (Time = '+string(systime(1)-syst, format='(F0.5)')+')'
    syst = systime(1)
    
    sf = obj_valid(osurfs)
    ptr_free, mask
    ptr_free, adt_spec
    
    xobjview, osurfs[where(sf eq 1)], /block, background=[0,0,0]
    obj_destroy, osurfs
    return
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mas_adt_sp_init, vdim, adt_spec=adt_spec

    common scan_data

    ci = project.ci
    vdim_effective = vdim

    junk = sp_make_sphere(vdim_effective, vert=vertlist, poly=polylist)
    result = mesh_decimate(vertlist, polylist, new_plist, vertices=new_vlist, percent_vertices=40)
    vertlist=new_vlist
    polylist=new_plist
    
    ;;for i = 0, 2 do vertlist[i,*] *= vdim_effective/2.0
    nverts = (size(vertList))[2]
    ;; shift to center
    for i = 0, 2 do begin
        vertlist[i,*] -= vdim_effective/2.0 ; + 0.5
    endfor
    
    for i = 0 ,nverts-1 do begin
        vertlist[*,i] /= norm(vertlist[*,i])
    endfor
    
    print, "Number of ADT vertices: ", nverts
    
    vertList0 = fltarr(3,nverts)
    phi_theta_r = cv_coord(from_rect=vertList, /to_sphere)
    
    b_thres = mean(*project.imndarray[ci].bval_array)
    hi_b = where(*project.imndarray[ci].bval_array gt b_thres)
    lo_b = where(*project.imndarray[ci].bval_array le b_thres)

    adt_spec = ptr_new({ $
                         hi_b: hi_b, $
                         lo_b: lo_b, $
                         phi_theta_r: phi_theta_r, $
                         vertlist:cv_coord(from_sphere=phi_theta_r, /to_rect), $
                         polylist:polylist, $
                         n_verts: nverts $
                       }, /no_copy)
    
end

; Subroutine name: mas_adt_superquadric
; Created by: BT, 2008-06
; Calling Information:

; vdim             - IN  - size of glyph
; gamma_in         - IN  - superquadric parameter that controls
;                          edgeiness of glyph
; shiny            - IN  - controls how shiny the glyph is
; colorize         - IN  - enable color
; use_roi_mask     - IN  - mask with selected roi. if this is set and
;                          there is no active roi, it is ignored
; crop_roi         - IN  - set this to crop the image to the roi boundary
; l0_intensity     - IN  - intensity of light
; l1_intensity             (same)
; l2_intensity             (same)
; bg_image_opacity - IN  - opacity of background image
; glyph_opacity    - IN  - opacity of glyph
; slice            - IN  - which slice to view (if not set, current
;                          slice will be used
; image_only
; display          - IN  - set this to have the resulting image
;                          displayed to the user
; image_out        - OUT - set this keyword to a named variable to
;                          and the resulting image will be returned
; xobjview_display - IN  - Show the result in xobjview

;
; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
pro mas_adt_superquadric, $
                          vdim             = vdim, $
                          vmult            = vmult, $
                          gamma_in         = gamma_in, $
                          shiny            = shiny, $
                          colorize         = colorize, $
                          use_roi_mask     = use_roi_mask, $
                          crop_roi         = crop_roi, $
                          fa_thr           = fa_thr, $
                          l0_intensity     = l0_intensity, $
                          l1_intensity     = l1_intensity, $
                          l2_intensity     = l2_intensity, $
                          bg_image_opacity = bg_image_opacity, $
                          glyph_opacity    = glyph_opacity, $
                          slice            = slice, $
                          image_only       = image_only, $
                          display          = display, $
                          image_out        = image_out, $
                          xobjview_display = xobjview_display

    common scan_data

;;;
;;; get some project characteristics
;;;

    ci = project.ci

    do_adt_regress_mono

    ;; ADT didn't get processed? - can't continue. ADT should have
    ;;                             notified the user if there was and
    ;;                             error. They also might have
    ;;                             cancelled it.
    if (project.procpramarray[ci].adt_proccess_flag eq 0) then return

    rint = project.procpramarray[ci].freq_interp
    pint = project.procpramarray[ci].phase_interp
    sint = project.procpramarray[ci].slice_interp

    fdim = project.imndarray[ci].fdim * rint
    pdim = project.imndarray[ci].pdim * pint
    sdim = project.imndarray[ci].sdim * sint

    slice_axis = project.procpramarray[ci].slice_axis

    if not keyword_set(vdim)   then vdim  = 10
    if not keyword_set(fa_thr) then fa_thr  = 1.0
    if not keyword_set(vmult)  then begin
       if (vdim ge 20) then vmult = 0.5 else vmult = 1.0
    endif

    if not keyword_set(gamma_in) then gamma = 0.0 else gamma = float(gamma_in)

    if n_elements(lmax) eq 0      then L=4 else L=lmax
    if n_elements(vdim) eq 0      then vdim = 10
    if n_elements(lambda) eq 0    then lambda = 0.1
    if n_elements(vmult) eq 0     then vmult = 1.0
    if n_elements(threshold) eq 0 then threshold = 0.0

    if (ptr_valid(project.procpramarray[ci].adt_sp_state)) then begin
        adt_sp_state = project.procpramarray[ci].adt_sp_state
    endif else begin
        adt_sp_state = mas_adt_sp_make_default_state()
    endelse

    if (ptr_valid(project.procpramarray[ci].difftools_state)) then begin
        difftools_state = project.procpramarray[ci].difftools_state
    endif else begin
        difftools_state = mdt_make_default_state()
    endelse

    gamma        = (*adt_sp_state).edge/5.0
    fa_thr       = (*adt_sp_state).fa_thr
    eval_divisor = (*adt_sp_state).eval_divisor > 1

    crop_roi     = (*difftools_state).crop_roi
    mask_roi     = (*difftools_state).mask_roi
    vdim         = (*difftools_state).glyph_size
    vmult        = (*difftools_state).glyph_size_mult
    l0_intensity = (*difftools_state).l0_intensity / 100.0
    l1_intensity = (*difftools_state).l1_intensity / 100.0
    l2_intensity = (*difftools_state).l2_intensity / 100.0
    image_opac   = (*difftools_state).image_opacity   / 100.0
    glyph_opacity= (*difftools_state).glyph_opacity   / 100.0
    glyph_shiny  = (*difftools_state).glyph_shiny
    use_adt_color = (*difftools_state).use_adt_color

    print, 'mas_adt_superquadric: starting with gamma='+strcompress(string(gamma), /remove_all)

;;;
;;; Determine the dimensions of the input and output 
;;; images. Also use the selected slice if none is specified 
;;; in the parameters
;;;
    voxel_dim = [ project.imndarray[ci].f_voxsz/project.procpramarray[ci].freq_interp , $
                  project.imndarray[ci].p_voxsz/project.procpramarray[ci].phase_interp, $
                  project.imndarray[ci].s_voxsz/project.procpramarray[ci].slice_interp ]
    voxel_dim /= min(voxel_dim)

    case slice_axis of 
        0: begin
            xdim = fdim
            ydim = pdim 
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].p_fov
            sl = keyword_set(slice) ? slice : project.procpramarray[ci].sdim_start
            xstart = 0
            ystart = 0
        end
        1: begin
            xdim = fdim; * voxel_dim[0] 
            ydim = sdim; * voxel_dim[2]
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].s_fov
            sl = keyword_set(slice) ? slice : project.procpramarray[ci].pdim_start
            xstart = 0
            ystart = 0
        end
        2: begin
            xdim = pdim
            ydim = sdim
            xfov = project.imndarray[ci].p_fov
            yfov = project.imndarray[ci].s_fov
            sl = keyword_set(slice) ? slice : project.procpramarray[ci].fdim_start
            xstart = 0
            ystart = 0
        end
    endcase

    ;; Obtain the background image from ADT routines
    tmp = ptr_new(adt_make_image(sl,0), /no_copy)
    
    case ((size(*tmp))[0]) of
        2: begin
            ;; adt_make_image returned a b/w image, use that
            bg_image = *tmp
            
            if (project.procpramarray[project.ci].adt_display_type eq 4) then begin
               ;; adt_make_image is inconsistent when it comes to the
               ;; diffusion-weighted image, so we un-de-rotate it. why
               ;; does this happen?
                if (project.procpramarray[project.ci].flip_direction gt 0) then begin
                    bg_image = reverse(bg_image, project.procpramarray[project.ci].flip_direction)
                endif
                
                ;; inverse rotation
                bg_image = rotate(bg_image, (4-(project.procpramarray[project.ci].rotate_direction)) mod 4)
                
            endif
        end
        3: begin
            bg_image = bytscl(*tmp)
        end
        else: begin
            junk = dialog_message(['Could not acquire ADT image.', $
                                   'Please check ADT display selection'], $
                                  /center, /error)
            return
        end
    endcase
    ptr_free, tmp

   ;; Apply ROI Masking and cropping, recompute dimensions
    if (keyword_set(mask_roi)) then begin
        if (keyword_set(crop_roi)) then begin
            mask = mas_roi_get_current_mask([xdim, ydim], crop=crop_dims, /no_transform)
        endif else begin
            mask = mas_roi_get_current_mask([xdim, ydim], /no_transform)
        endelse

        if (not ptr_valid(mask)) then begin
            mask_roi = 0
            if (n_elements(crop_dims) ne 0) then begin
                junk = temporary(crop_dims)
            endif
        endif
    endif

    if (n_elements(crop_dims) ne 0) then begin
        crop = temporary(bg_image[ crop_dims[0]:crop_dims[2]-1, $
                                   crop_dims[1]:crop_dims[3]-1, $
                                   * ])
        bg_image = temporary(crop)

        crop = (*mask)[ crop_dims[0]:crop_dims[2]-1, $
                        crop_dims[1]:crop_dims[3]-1 ]
        ptr_free, mask
        mask = ptr_new(crop, /no_copy)
        
        old_xdim = xdim
        old_ydim = ydim

        xdim = crop_dims[2]-crop_dims[0]
        ydim = crop_dims[3]-crop_dims[1]
        
        xfov *= float(xdim)/float(old_xdim)
        yfov *= float(ydim)/float(old_ydim)

        xstart = crop_dims[0]
        xend   = crop_dims[2]-1
        ystart = crop_dims[1]
        yend   = crop_dims[3]-1

    endif

   ;; recompute optimal vdim based on image dimensions
    vdim_opt = mas_glyphscene_renderer_getOptimalVdim(xdim, ydim, vdim)

    if (vdim_opt ne vdim) then begin
        vdim = vdim_opt
;        junk = dialog_message(['Glyph size too large for render buffer.', $
;                               'Reducing to '+string(vdim_opt, format='(I0.0)')+'.'], $
;                              /center)
    endif

    ;; at least 20 vertices are required for a decent shape
    if (vdim lt 24) then vmult = 1.0
    if (gamma ne 0.0) then begin
;        vmult += 0.75
    endif

    gs = obj_new('mas_glyphscene_renderer', bg_image=bg_image, $
                 vmult=vmult, $
                 vdim=vdim, $
                 light0_int=l0_intensity, $
                 light1_int=l1_intensity, $
                 light2_int=l2_intensity, $
                 bg_img_opac=image_opac, $
                 glyph_opac=glyph_opacity, $
                 /skip_bufsize_warnings)

    prendered_img = ptr_new(bytarr(xdim*vdim, ydim*vdim, 3), /no_copy)
    vdim = vmult*vdim

    ;; In this case, we don't need to create a superquadric --
    ;; just a sphere will do
    if gamma eq 0.0 then begin
       surf = sp_make_sphere(vdim, vert=ell_verts, poly=ell_polys)
    endif

    osurfs_xy      = objarr(xdim,ydim)
    syst = systime(1)
    curr_iter = 0

    glyph_color    = [50,50,100] ;fltarr(3)
    glyph_diffuse  = [90,90,110]
    glyph_specular = [150,150,170]

    if (n_elements(glyph_opacity) eq 0) then glyph_opacity = 1.0

    if (keyword_set(shiny)) then begin
        glyph_shiny = float(shiny)
    endif else begin
        glyph_shiny = (gamma eq 0.0) ? 5.0 : 7.0 ;; "shininess' property
    endelse

    glyph_shading = 1         ;; shading property
    
    junk = get_orient_rgb_axis(tx_matrix=sa_tx)

    progressBar = obj_new('progressbar', title='Progress Bar', text='Creating Glyph Image...')
    progressBar->Start

    for x = 0, xdim-1 do begin 
        for y =  0, ydim-1 do begin 
            
            x_ind = x + xstart
            y_ind = y + ystart

            nevermind = 0
            curr_iter++

            ;; check the mask.
            if (keyword_set(mask_roi)) then begin
                if (*mask)[x,y] eq 0 then continue
            endif

            ;; get the eigvenvalues
            case slice_axis of
                0: begin
                    if (*project.dataarray[ci].state1)[x_ind,y_ind,sl,0] eq 0 then nevermind = 1
                    evals = reform((*project.dataarray[ci].eign_val)[x_ind,y_ind,sl,*])
                    fa    = reform((*project.dataarray[ci].frac_ani)[x_ind,y_ind,sl])
                    evecs = transpose(reform((*project.dataarray[ci].eign_vec)[x_ind,y_ind,sl,*,*]))
                end
                1: begin
                    if (*project.dataarray[ci].state1)[x_ind,sl,y_ind,0] eq 0 then nevermind = 1
                    evals = reform((*project.dataarray[ci].eign_val)[x_ind,sl,y_ind,*])
                    fa    = reform((*project.dataarray[ci].frac_ani)[x_ind,sl,y_ind])
                    evecs = transpose(reform((*project.dataarray[ci].eign_vec)[x_ind,sl,y_ind,*,*]))
                end
                2: begin
                    if (*project.dataarray[ci].state1)[sl,x_ind,y_ind,0] eq 0 then nevermind = 1
                    evals = reform((*project.dataarray[ci].eign_val)[sl,x_ind,y_ind,*])
                    fa    = reform((*project.dataarray[ci].frac_ani)[sl,x_ind,y_ind])
                    evecs = transpose(reform((*project.dataarray[ci].eign_vec)[sl,x_ind,y_ind,*,*]))
                end
            endcase

            if (total(evals) eq 0 or total(finite(evals, /nan)) gt 0 or nevermind eq 1) then continue

            ;; normalize so that eval[0] = 1
            evals = abs(evals/max(evals))
            ;; eval_divisor used to diminish the non-primary evs so
            ;; the main direction is more prominent
            evals = [evals[0], evals[1]/eval_divisor, evals[2]/eval_divisor]

            if (fa lt fa_thr) then begin evals = evals * 0.1 & continue & endif

            if (max(evals) ne evals[0]) then begin
                print, "Unordered eigenvalue error at ("+strcompress(string(x_ind)+','+string(y_ind), /remove_all)+') ', evals
                continue
            endif
            
            ;;  evec matrix use as
            ;; transformation matrix
            evecs = [[evecs[*,0], 0], $
                     [evecs[*,1], 0], $
                     [evecs[*,2], 0], $
                     [0, 0, 0, 1] ]

            ;; get the surface
            if (gamma ne 0.0) then begin
                surf = sp_superquadric(vdim, evals, gamma, poly=polylist, vert=vertlist)
                if (surf eq 0) then continue
            endif else begin
                polylist = ell_polys
                vertlist = ell_verts
            endelse
            
            if (total(finite(vertList, /nan)) ne 0) then continue

            if curr_iter mod 20 eq 0 then progressBar->Update, float(curr_iter)/(xdim*ydim)*100.0
            
            ;; transformation matrix: evecs (rotation) * evals (scale)
            ;; * slice_axis transformation
            ;; NOTE NOTE NOTE: because Bruker's b-matrix and
            ;; its angle_{theta,phi} are reported differently,
            ;; we DON'T apply acq_matrix here. Also note that
            ;; for PAR/REC files, we adjust the b-matrix when the par
            ;; files is read.
            tx = evecs ## diag_matrix([evals, 1.0]) # sa_tx

            ;; before applying the trans. matrix, we'd like to have
            ;; the center of the surface at the origin. So we make a
            ;; shift array
            sz = (size(vertlist))[2]
            shift = transpose([ [replicate(-(vdim/2.), sz)], $
                                [replicate(-(vdim/2.), sz)], $
                                [replicate(-(vdim/2.), sz)] ])
            ;; shift, rotate, de-shift
            vertList += shift
            vertList = vert_t3d(vertList, /no_copy, matrix=tx)
            vertList -= shift

            ;; this will translate the vertex coordinates to their
            ;; position relative to the background image
            shift = transpose([ [replicate(0, sz)], $
                                [replicate(y*vdim, sz)], $
                                [replicate(0, sz)] ])
            vertList += shift

            ;; compute the glyph color based on FA * principal evector
            if (keyword_set(use_adt_color) ne 0) then begin
                glyph_color = adt_get_orientation_color(x_ind, y_ind, sl)
            endif
            
            osurf = obj_new('idlgrpolygon', vertList, polygons=polylist, $
                            name      = strcompress('Glyph: ('+string(x)+','+string(y)+')'), $
                            specular  = glyph_specular, $
                            diffuse   = glyph_diffuse, $ ;glyph_color, $
                            ambient   = glyph_color, $
                            shading   = glyph_shading, $
                            color     = glyph_color, $
                            shininess = glyph_shiny,$
                            alpha_channel = glyph_opacity, $
                            reject    = 0, $
                            texture_interp  = 1) 
            ;; accumulate
            gs->addGlyph, osurf, x, y, /replace
            osurfs_xy[x,y] = osurf
            
            if (progressBar->checkcancel()) then begin
                progressBar->destroy
                obj_destroy, gs
                if (ptr_valid(mask)) then ptr_free, mask
                obj_destroy, osurfs_xy
                return
            endif
            
        endfor

    endfor

    progressBar->destroy

    gs->render, rendered_img=image_out, use_this_ptr=prendered_img
    
    obj_destroy, gs
    obj_destroy, osurfs_xy
    if (ptr_valid(mask)) then ptr_free, mask

    if (not ptr_valid(prendered_img)) then return
    
    mas_rotate_flip, prendered_img; image_out
    image_out = temporary(*prendered_img) ;; gets "returned"
    ptr_free, prendered_img

    if (keyword_set(use_adt_color)) then begin
        orient_circle = circle(150,150,1)
        pcircle = ptr_new(orient_circle, /no_copy)
        mas_rotate_flip, pcircle
    endif

    mas_display_multi, image_out, fov_x=xfov, fov_y=yfov, fov_units=1, $
      tab_title='ADT Scan '+strcompress(string(project.ci+1), /remove_all)+ $
      '/Slice '+strcompress(string(sl), /remove_all), $
      extra_thumbnail=ptr_valid(pcircle) ? temporary(*pcircle) : 0

    return

end
