;; $Id$
;;

;; Subroutine name: mas_dot_gui_make_default_state
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Creates a pointer to a state structure containing the default
;;     DOT parametrs
;;
;; Editing Information:

function mas_dot_gui_make_default_state

    common scan_data

    if (not ptr_valid(project.imndarray[project.ci].big_delta)) then return, ptr_new()

    return, ptr_new({ DOT_STATE, $
                      dot_sav_flt:0,    $
                      dot_diff_time: (*project.imndarray[project.ci].big_delta)[0], $
                      dot_diff_sep:  (*project.imndarray[project.ci].small_delta)[0], $
                      dot_load_flt:0,   $
                      dot_signal_thr:5,  $
                      dot_anis_thr:0.001, $
                      dot_res_idx:0,    $
                      dot_rzero:10.0,   $
                      dot_lmax:4,       $
                      dot_bval:400,     $
                      dot_img_ct: 1,    $
                      dot_glyph_ct: 3,  $
                      dot_glyph_size:10,$
                      dot_glyph_size_mult:0.5, $
                      disp_ga_only:0,   $
                      use_obj_gfx:1,    $
                      bg_image:0,       $
                      ga_stat:0    $
                      
                      
                       ;; 0 = "GA", 1 = "Current ADT Selection"
                    }, /no_copy)

end

; Subroutine name: mas_dot_make_initial_surface
; Created by: BT, 2008-06
; Calling Information:
;
; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Creates the vertices that will be the "domain" for each probability
; distribution function. The output will be an array of vertices and
; the polygon list. The vertices will be altered based on the
; probabilities thus creating the probability surface

; Editing Information:

;; pro mas_dot_make_initial_surface, $
;;                                  dims, $
;;                                  lmax=lmax, $
;;                                  vdim=vdim, $
;;                                  surface_def=surface_def
                                 
;;     common scan_data
    
;;     xdim=dims[0]
;;     ydim=dims[1]
;;     pdim=dims[2]
    
;;     if (n_elements(lmax)       eq 0) then lmax = (2*(sqrt(pdim)-1) < 6)

;;     volume = fltarr(vdim,vdim,vdim)

;;     for x = 0,vdim-1 do begin
;;         for y = 0,vdim-1 do begin
;;             for z = 0,vdim-1 do begin
;;                 volume[x,y,z] = -total( (vdim/2.0 - 0.5 - [x,y,z])^2 )
;;             endfor
;;         endfor
;;     endfor
    
;;     SCALE3, XRANGE=[0, vdim], YRANGE=[0, vdim], ZRANGE=[0, vdim]

;;     st = -(vdim / 4.0)^2 * 3.0
;;     shade_volume, volume, st, vertList, polyList, low=1
;;     nverts = (size(vertList))[2]

;;     for i = 0, 2 do begin
;;         vertlist[i,*] -= vdim/2.0 + 0.5
;;     endfor

;;     vertList0 = fltarr(3,nverts)
;;     phi_theta_r = cv_coord(from_rect=vertList, /to_sphere)
    
;;     ;; normalize
;;     for i = 0L, nverts-1 do begin
;;         vertList0[*,i] = vertList[*,i]/phi_theta_r[2,i]
;;     endfor
    
;;     verts_theta0 = -phi_theta_r[1,*] + !pi/2.0
;;     verts_phi0   =  phi_theta_r[0,*] 
    
;;     l_arr = [0,2,2,2,4,4,4,4,4,REPLICATE(6,7),REPLICATE(8,9),REPLICATE(10,11)]
;;     m_arr = [0,0,1,2,0,1,2,3,4,INDGEN(7),INDGEN(9),INDGEN(11)]
        
;;     ;; this is a weighting factor that doubles the weighting of non-zero m values, since -m values
;;     ;; are not calculated because of equation 11 in Ozarslan 2006
;;     w_arr = 2.D - DOUBLE(m_arr eq 0)
;;     N_lm  = (lmax/2+1)^2

;;     spharm = dcomplexarr(nverts, N_lm)
;;     for ll = 0, N_lm-1 do begin
;;         spharm[*, ll] = SPHER_HARM(verts_theta0, verts_phi0, l_arr[ll], m_arr[ll])
;;     endfor

;;     surface_def = ptr_new({ vertList:vertList0, $
;;                             polyList:polyList,  $
;;                             vdim: vdim, $
;;                             spharm: spharm, $
;;                             w_arr:w_arr, $
;;                             N_lm:N_lm }, /no_copy)

;; end

;; Subroutine name: mas_dot_color_mapper
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;     IMAGE:  The 2D image data
;;
;;     {X,Y}DIM: the dimensions of the image
;;
;;     TABLE: IDL color table to use 
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Colorize a b/w image using an IDL color table
;;
;; Editing Information:

function mas_dot_color_mapper, image, xdim, ydim, table

    color_im = bytarr(xdim, ydim, 3)

    LOADCT, table
    TVLCT, r, g, b, /get

    pix_idx = image[indgen(xdim*ydim, /LONG)]

    color_im[*,*,0] = reform(temporary(r[pix_idx]), xdim, ydim)
    color_im[*,*,1] = reform(temporary(g[pix_idx]), xdim, ydim)
    color_im[*,*,2] = reform(temporary(b[pix_idx]), xdim, ydim)
    
    return, color_im

end


;; Subroutine name: mas_dot_compute_ga
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;     DOT_OBJ: DOT reconstructor object
;;   
;;     VIS_SLICE: slice to use. If not set, current MAS selection
;;                according to slice axis is used
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; 
;;     Computes Generalized Anisotropy from DOT weights
;;
;; Editing Information:

function mas_dot_compute_ga, dot_obj, vis_slice=vis_slice
    
    common scan_data
    
    ci = project.ci
    n_slices = n_elements(vis_slice)
    dot_obj->getProperty, l_max=lmax

    if (n_slices eq 0) then begin
       case project.procpramarray[ci].slice_axis of
          0: slice = project.procpramarray[ci].sdim_start
          1: slice = project.procpramarray[ci].pdim_start
          2: slice = project.procpramarray[ci].fdim_start
       endcase
    endif else begin
       slice = vis_slice[0]
    endelse
    data_sz = size(*project.dataarray[ci].state1, /dimensions)

    case project.procpramarray[ci].slice_axis of
       0: image = reform((*project.dataarray[ci].state1)[*,*,slice,*])
       1: image = reform((*project.dataarray[ci].state1)[slice,*,*,*])
       2: image = reform((*project.dataarray[ci].state1)[*,slice,*,*])
    endcase

    xdim = (size(image, /dimensions))[0]
    ydim = (size(image, /dimensions))[1]
    N_lm  = (lmax/2+1)^2
    sum = FLTARR(xdim,ydim)
    ga  = FLTARR(xdim,ydim)
    plm = dcomplexarr(xdim,ydim, N_lm)

    for x = 0, xdim-1 do begin
       for y = 0, ydim-1 do begin

          tmp = dot_obj->computeWeights(reform(image[x,y,*]))
          if total(finite(tmp, /nan)) ne 0 then continue
          plm[x,y,*] = tmp
          
       endfor
    endfor

    ww  = WHERE(plm[*,*,0] NE 0.)
    
    l_arr = [0,2,2,2,4,4,4,4,4,REPLICATE(6,7),REPLICATE(8,9),REPLICATE(10,11)]
    m_arr = [0,0,1,2,0,1,2,3,4,INDGEN(7),INDGEN(9),INDGEN(11)]
        
    ;; this is a weighting factor that doubles the weighting of non-zero m values, since -m values
    ;; are not calculated because of equation 11 in Ozarslan 2006
    w_arr = 2.D - DOUBLE(m_arr eq 0)

    IF ww[0] EQ -1 THEN  RETURN, -1
    
    FOR ll = 0, N_lm-1 DO BEGIN
        sum[ww] = sum[ww]+w_arr[ll]*ABS((plm[*,*,ll])[ww])^2
    ENDFOR
    
    ga[ww] = (sum[ww]/((plm[*,*,0])[ww]^2)-1.)/9.
    ga[ww] = 1. - 1./(1.+(250.*ga[ww])^(1.+1./(1.+5000.*ga[ww]))) > 1./255.
    
    return, ga

end
 
;; Subroutine name: getTesAngles
;; Created by: CD 7/03/07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; returns vectors of tesselation vertices in polar coordinates (dimensions 1 and 2)
;; based on tesselation level requested
;; not sure what third dimension is
;; data from files from Evren Ozarslan
;;
;; third dimension contains integration weights for that particular
;; gradient direction. They were supposedly calculated by computing the
;; areas of the voronoi polygons determined by each gradient direction. However,
;; there is no readily available code to confirm how they were
;; computed. Also, setting them all to 1 in the computation has a very
;; negligible effect. (BT 2008-06-23)
;;
;; Editing Information:

;; function getTesAngles, level

;;     case level of 
        
;;         0: begin
;;             tes_angles = [[0.000000,63.4350,63.4350,63.4350,63.4350,63.4350,116.565, $
;;                            116.565,116.565,116.565,116.565,180.000],$
;;                           [90.0000,0.000000,72.0000,144.000, $
;;                            216.000,288.000,36.0000,108.000,180.000,252.000,324.000,180.000],$
;;                           [1.04720, $
;;                            1.04720,1.04720,1.04720,1.04720,1.04720,1.04720,1.04720,1.04720,1.04720, $
;;                            1.04720,1.04720]]
;;         end

;;         1: begin
;;             tes_angles = [[0.000000,63.4350,63.4350,63.4350,63.4350,63.4350,116.565, $
;;                            116.565,116.565,116.565,116.565,180.000,31.7175,31.7175,31.7175,31.7175, $
;;                            31.7175,58.2826,58.2826,90.0000,90.0000,58.2826,90.0000,90.0000,58.2826, $
;;                            90.0000,90.0000,58.2826,90.0000,90.0000,90.0000,90.0000,121.717,121.717, $
;;                            148.283,121.717,148.283,121.717,148.283,121.717,148.283,148.283], $
;;                           [90.0000, $
;;                            0.000000,72.0000,144.000,216.000,288.000,36.0000,108.000,180.000,252.000, $
;;                            324.000,180.000,0.000000,72.0000,144.000,216.000,288.000,36.0000,324.000, $
;;                            18.0000,342.000,108.000,54.0000,90.0000,180.000,126.000,162.000,252.000, $
;;                            198.000,234.000,-90.0000,306.000,72.0000,0.000000,36.0000,144.000,108.000, $
;;                            216.000,180.000,288.000,252.000,324.000], $
;;                           [0.273844,0.273844,0.273844, $
;;                            0.273844,0.273844,0.273844,0.273846,0.273846,0.273846,0.273846,0.273846, $
;;                            0.273836,0.309342,0.309342,0.309342,0.309342,0.309342,0.309341,0.309341, $
;;                            0.309340,0.309340,0.309341,0.309340,0.309340,0.309341,0.309340,0.309340, $
;;                            0.309341,0.309340,0.309340,0.309340,0.309340,0.309343,0.309343,0.309343, $
;;                            0.309343,0.309343,0.309343,0.309343,0.309343,0.309343,0.309343]]
;;         end

;;         2: begin
;;             tes_angles = [ $
;;                            [0.000000,63.4350,63.4350,63.4350,63.4350,63.4350,116.565, $
;;                             116.565,116.565,116.565,116.565,180.000,21.1450,42.2900,37.3774,37.3774, $
;;                             21.1450,42.2900,37.3774,21.1450,42.2900,37.3774,21.1450,42.2900,37.3774, $
;;                             21.1450,42.2900,58.8818,58.8818,79.1877,58.8818,58.8818,79.1877,81.0208, $
;;                             98.9792,100.812,81.0208,98.9792,58.8818,58.8818,79.1877,81.0208,98.9792, $
;;                             100.812,81.0208,98.9792,58.8818,58.8818,79.1877,81.0208,98.9792,100.812, $
;;                             81.0208,98.9792,58.8818,58.8818,79.1877,81.0208,98.9792,100.812,81.0208, $
;;                             98.9792,81.0208,98.9792,100.812,81.0208,98.9792,121.118,121.118,142.623, $
;;                             121.118,121.118,142.623,137.710,158.855,121.118,121.118,142.623,137.710, $
;;                             158.855,121.118,121.118,142.623,137.710,158.855,121.118,121.118,142.623, $
;;                             137.710,158.855,137.710,158.855] , $
;;                            [90.0000,0.000000,72.0000,144.000,216.000, $
;;                             288.000,36.0000,108.000,180.000,252.000,324.000,180.000,0.000000,0.000000, $
;;                             36.0000,324.000,72.0000,72.0000,108.000,144.000,144.000,180.000,216.000, $
;;                             216.000,252.000,288.000,288.000,23.6244,48.3756,36.0000,336.376,311.624, $
;;                             324.000,12.3957,23.6043,0.000000,347.604,336.396,95.6244,120.376,108.000, $
;;                             59.6043,48.3957,72.0000,84.3957,95.6043,167.624,192.376,180.000,131.604, $
;;                             120.396,144.000,156.396,167.604,239.624,264.376,252.000,203.604,192.396, $
;;                             216.000,228.396,239.604,275.604,264.396,288.000,300.396,311.604,59.6244, $
;;                             84.3756,72.0000,12.3756,347.624,0.000000,36.0000,36.0000,131.624,156.376, $
;;                             144.000,108.000,108.000,203.624,228.376,216.000,180.000,180.000,275.624, $
;;                             300.376,288.000,252.000,252.000,324.000,324.000] , $
;;                            [0.122797,0.122796, $
;;                             0.122796,0.122796,0.122796,0.122796,0.122796,0.122796,0.122796,0.122795, $
;;                             0.122795,0.122797,0.137284,0.137283,0.142790,0.142790,0.137284,0.137284, $
;;                             0.142790,0.137284,0.137283,0.142790,0.137284,0.137283,0.142790,0.137284, $
;;                             0.137283,0.137284,0.137284,0.142790,0.137284,0.137284,0.142790,0.137283, $
;;                             0.137283,0.142789,0.137282,0.137283,0.137285,0.137284,0.142790,0.137283, $
;;                             0.137283,0.142789,0.137283,0.137283,0.137284,0.137284,0.142790,0.137282, $
;;                             0.137282,0.142789,0.137283,0.137283,0.137284,0.137284,0.142790,0.137283, $
;;                             0.137283,0.142789,0.137282,0.137283,0.137283,0.137283,0.142789,0.137283, $
;;                             0.137283,0.137285,0.137285,0.142791,0.137286,0.137285,0.142791,0.137284, $
;;                             0.137282,0.137286,0.137286,0.142791,0.137284,0.137282,0.137285,0.137286, $
;;                             0.142791,0.137284,0.137282,0.137286,0.137286,0.142791,0.137284,0.137282, $
;;                             0.137284,0.137282]]
;;         end
        
;;         3: begin
;;             tes_angles = [[ $
;;                             0.000000,63.4350,63.4350,63.4350,63.4350,63.4350,116.565,116.565,116.565,116.565,116.565,180.000,31.7175,31.7175, $
;;                             31.7175,31.7175,31.7175,58.2826,58.2826,90.0000,90.0000,58.2826,90.0000,90.0000,58.2826,90.0000,90.0000,58.2826, $
;;                             90.0000,90.0000,90.0000,90.0000,121.717,121.717,148.283,121.717,148.283,121.717,148.283,121.717,148.283,148.283, $
;;                             15.8587,15.8587,15.8587,15.8587,15.8587,47.5762,59.6208,59.6208,76.5584,76.5584,47.5762,59.6208,59.6208,76.5584, $
;;                             76.5584,47.5762,59.6208,59.6208,76.5584,76.5584,47.5762,59.6208,59.6208,76.5584,76.5584,47.5762,59.6208,59.6208, $
;;                             76.5584,76.5584,103.442,103.442,120.379,120.379,132.424,103.442,103.442,120.379,120.379,132.424,103.442,103.442, $
;;                             120.379,120.379,132.424,103.442,103.442,120.379,120.379,132.424,103.442,103.442,120.379,120.379,132.424,164.141, $
;;                             164.141,164.141,164.141,164.141,26.5651,26.5651,43.6470,43.6470,26.5651,43.6470,43.6470,26.5651,43.6470,43.6470, $
;;                             26.5651,43.6470,43.6470,43.6470,43.6470,73.9550,73.9550,73.9550,73.9550,90.0000,90.0000,106.045,90.0000,106.045, $
;;                             73.9550,73.9550,90.0000,106.045,90.0000,106.045,73.9550,73.9550,90.0000,106.045,90.0000,106.045,73.9550,73.9550, $
;;                             90.0000,106.045,90.0000,106.045,90.0000,106.045,106.045,136.353,136.353,136.353,136.353,153.435,153.435,136.353, $
;;                             136.353,153.435,136.353,136.353,153.435,136.353,136.353,153.435], $
;;                           [ $
;;                             90.0000,0.000000,72.0000,144.000,216.000,288.000,36.0000,108.000,180.000,252.000,324.000,180.000,0.000000,72.0000, $
;;                             144.000,216.000,288.000,36.0000,324.000,18.0000,342.000,108.000,54.0000,90.0000,180.000,126.000,162.000,252.000, $
;;                             198.000,234.000,-90.0000,306.000,72.0000,0.000000,36.0000,144.000,108.000,216.000,180.000,288.000,252.000,324.000, $
;;                             0.000000,72.0000,144.000,216.000,288.000,0.000000,17.5330,342.467,9.50570,350.494,72.0000,54.4670,89.5330,62.4943, $
;;                             81.5057,144.000,126.467,161.533,134.494,153.506,216.000,198.467,233.533,206.494,225.506,288.000,305.533,270.467, $
;;                             278.494,297.506,26.4943,45.5057,53.5330,18.4670,36.0000,98.4943,117.506,90.4670,125.533,108.000,170.494,189.506, $
;;                             162.467,197.533,180.000,242.494,261.506,234.467,269.533,252.000,333.506,314.494,341.533,306.467,324.000,36.0000, $
;;                             108.000,180.000,252.000,324.000,36.0000,324.000,22.3862,337.614,108.000,49.6138,94.3862,180.000,121.614,166.386, $
;;                             252.000,193.614,238.386,310.386,265.614,26.2677,45.7323,333.732,314.268,6.28397e-006,36.0000,9.73230,324.000,350.268, $
;;                             98.2677,117.732,72.0000,62.2677,108.000,81.7323,170.268,189.732,144.000,134.268,180.000,153.732,242.268,261.732, $
;;                             216.000,206.268,252.000,225.732,288.000,278.268,297.732,58.3862,85.6138,13.6138,346.386,72.0000,1.20440e-005,130.386, $
;;                             157.614,144.000,202.386,229.614,216.000,274.386,301.614,288.000], $
;;                           [ $
;;                             0.0692920,0.0692920,0.0692920,0.0692920,0.0692920,0.0692920,0.0692910,0.0692910,0.0692910,0.0692910,0.0692910,0.0692940,0.0782400,0.0782400, $
;;                             0.0782390,0.0782390,0.0782390,0.0782400,0.0782390,0.0782400,0.0782400,0.0782390,0.0782400,0.0782400,0.0782390,0.0782400,0.0782400,0.0782390, $
;;                             0.0782390,0.0782400,0.0782400,0.0782400,0.0782400,0.0782390,0.0782390,0.0782390,0.0782390,0.0782390,0.0782390,0.0782390,0.0782390,0.0782390, $
;;                             0.0741900,0.0741900,0.0741900,0.0741900,0.0741900,0.0741900,0.0741900,0.0741900,0.0741900,0.0741910,0.0741900,0.0741900,0.0741900,0.0741900, $
;;                             0.0741900,0.0741900,0.0741900,0.0741900,0.0741910,0.0741910,0.0741900,0.0741900,0.0741900,0.0741910,0.0741910,0.0741900,0.0741900,0.0741900, $
;;                             0.0741910,0.0741910,0.0741900,0.0741900,0.0741900,0.0741900,0.0741920,0.0741900,0.0741910,0.0741900,0.0741900,0.0741910,0.0741910,0.0741910, $
;;                             0.0741900,0.0741900,0.0741910,0.0741910,0.0741910,0.0741900,0.0741900,0.0741910,0.0741910,0.0741910,0.0741900,0.0741900,0.0741910,0.0741900, $
;;                             0.0741900,0.0741900,0.0741900,0.0741900,0.0822710,0.0822720,0.0822710,0.0822710,0.0822710,0.0822710,0.0822710,0.0822710,0.0822710,0.0822710, $
;;                             0.0822720,0.0822710,0.0822710,0.0822710,0.0822710,0.0822710,0.0822710,0.0822700,0.0822700,0.0822710,0.0822720,0.0822700,0.0822720,0.0822700, $
;;                             0.0822710,0.0822710,0.0822710,0.0822700,0.0822720,0.0822700,0.0822700,0.0822700,0.0822710,0.0822690,0.0822720,0.0822690,0.0822700,0.0822700, $
;;                             0.0822710,0.0822690,0.0822720,0.0822690,0.0822710,0.0822690,0.0822690,0.0822720,0.0822720,0.0822720,0.0822720,0.0822710,0.0822710,0.0822720, $
;;                             0.0822720,0.0822710,0.0822720,0.0822720,0.0822710,0.0822720,0.0822730,0.0822710]]
;;         end
        
;;         else: tes_angles = [0,0,0]
        
;;     endcase
    
;;     return, tes_angles
;; end

;; Subroutine name: match_angles
;; Created by: E. Ozarslan
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; given arrays of theta and phi angles, returns the locations of
;; the angles in the tabulated files
;;
;; Editing Information:

;; function match_angles, tes_level, theta_angle, phi_angle, a

;;     on_error, 2

;;     precision = 10.

;;     bt_listed = ROUND(REFORM(a[*,0])* precision)
;;     bp_listed = ROUND(REFORM(a[*,1])* precision)
    
;;     bt = ROUND(theta_angle*precision)
;;     bp = ROUND(phi_angle*precision)
    
;;     n_angles = n_elements(bt)
;;     locations = intarr(n_angles)
    
;;     for i = 0, n_angles-1 do begin

;;         checks = BYTARR(N_ELEMENTS(bt_listed))
;;         wt = WHERE(ABS(bt_listed - bt[i]) LE 1, wt_ct)
;;         wp = WHERE(ABS(bp_listed - bp[i]) LE 1, wp_ct)

;;         checks[wt] = checks[wt]+1
;;         checks[wp] = checks[wp]+1
        
;;         locations[i] = WHERE(checks EQ 2, lll)
        
;;         IF lll NE 1 THEN BEGIN
            
;;             ;PRINT, 'error, number of matches=', lll, '   at i=', i
;;             ;PRINT, 'theta=', bt[i]
;;             ;PRINT, 'phi=', bp[i]
            
;;             IF bt[i] EQ 0 THEN locations[i]=0
            
;;             ;;PRINT, i, '       theta=', bt[i], bt_listed[locations[i]], $
;;             ;;          '      phi=',    bp[i], bp_listed[locations[i]], locations[i]

;;         ENDIF

;;     ENDFOR
    
;;     RETURN, locations

;; END

;; Subroutine name: sht_weights_icosahedral
;; Created by: E. Ozarslan on 02/09/2005.
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; This function evaluates the weights (s_i) associated with the
;; spherical harmonic transform of a function whose values at
;; the points theta_angle and phi_angle are provided. These angles
;; are assumed to be those produced by the tes.pro routine and the
;; integration weights are calculated by angles.pro algorithm.
;; The returned weights are :
;; 	s_i^lm = w_i Y_lm(theta_i, phi_i)^*
;;
;; These weights are stored in the double precision complex
;; array s_arr, whose dimensions are
;;   10(tes_level+1)^2 +2 by (lmax/2+1)^2
;;
;; Only weights corresponding to l-even and m>=0 are calculated.
;;
;; Editing Information:

;; FUNCTION sht_weights_icosahedral, tes_level, theta_angle, phi_angle, lmax=lmax

;;     on_error, 2
    
;;     a = getTesAngles(tes_level)
    
;;     tha = theta_angle*!dtor
;;     pha = phi_angle*!dtor
    
;; ;;    locations = match_angles(tes_level, theta_angle, phi_angle, a)
    
;;     n_of_lms = (lmax/2+1)^2
;;     n_of_is  = 5*(tes_level+1)^2+1 ; hemisphere

;; ;;    n_of_is  = n_elements(tha)

;;     l_arr=[0,2,2,2,4,4,4,4,4,REPLICATE(6,7),REPLICATE(8,9),REPLICATE(10,11)]
;;     m_arr=[0,0,1,2,0,1,2,3,4,INDGEN(7),INDGEN(9),INDGEN(11)]
;;     s_arr=DCOMPLEXARR(n_of_is, n_of_lms)
    
;;     FOR lm=0, n_of_lms-1 DO BEGIN
;;         s_arr[*,lm]=  CONJ(SPHER_HARM(tha,pha,l_arr[lm],m_arr[lm]))
;; ;;        s_arr[*,lm]=a[locations, 2] * CONJ(SPHER_HARM(tha,pha,l_arr[lm],m_arr[lm]))
;;     ENDFOR
    
;;     RETURN, s_arr
;; END

;; Subroutine name: sph_har_transform_weights
;; Created by: E. Ozarslan and CD 6/15/07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; passthrough function, if not using diffusion directions corresponding
;; to tesselated ocosahedron, then need to modify to calculate integration
;; weights according to modified diffusion gradient direction scheme
;;
;; Editing Information:

;; FUNCTION sph_har_transform_weights, theta_angle, phi_angle, tes_level=tes_level, lmax=lmax

;;     on_error, 2

;;     IF N_ELEMENTS(lmax) EQ 0 THEN lmax=([2,4,8,10])[tes_level]
    
;;     IF N_ELEMENTS(tes_level) NE 0 THEN BEGIN
;;         RETURN, sht_weights_icosahedral(tes_level, theta_angle, phi_angle, lmax=lmax)
;;     ENDIF
    
;;     void = DIALOG_MESSAGE('DOT of this dataset is not supported.',/ERROR)
;;     RETURN, 0

;; END

;; All paper references in mas_DOT refer to following:
;; 1: Ozarslan et al. NeuroImage 31 (2006) 1086-1103

;; From Ozarslan, Parametric Reconstruction Step 2
;; use intermediate values to calculate Ilu
;; from Appendix A, Equation 28

;; reconstructed by CD, left out of original code
;; used by mas_DOT
;; updated by BT 20080605 - eliminate common block
;; function fast_fibers_prad, l,  para

;;     CASE l OF
;;         0:	BEGIN
;;             Al = 1
;;             Bl = 0
;;         END
;;         2: 	BEGIN
;;             Al = -(1 + 6 * (*para).bi2)
;;             Bl = 3
;;         END
;;         4:	BEGIN
;;             Al = 1 + 20 * (*para).bi2 + 210 * (*para).bi4
;;             Bl = 15/2 * ( 1 - 14 * (*para).bi2)
;;         END
;;         6: 	BEGIN
;;             Al = -(1 + 42 * (*para).bi2 + 1575/2 * (*para).bi4 + 10395 * (*para).bi6)
;;             Bl = 105/8 * (1 - 36 * (*para).bi2 + 396 * (*para).bi4)
;;         END
;;         8:	BEGIN
;;             Al = 1 + 72 * (*para).bi2 + 10395/4 * (*para).bi4 + 45045 * (*para).bi6 + 675675*(*para).bi8
;;             Bl = 315/16 * (1 - 66 * (*para).bi2 + 1716 * (*para).bi4 - 17160 * (*para).bi6)
;;         END
;;     ENDCASE
    
;;     RETURN, (Al * (*para).expo + Bl * (*para).errf)

;; end


;; Subroutine name: mas_dot_gui_event
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; event handler for configure gui
;;
;; Editing Information:

pro mas_dot_gui_event, event

    common scan_data
    
    ci = project.ci

    if ptr_valid(project.procpramarray[ci].difftools_state) then begin
        difftools_state = project.procpramarray[ci].difftools_state
    endif

    dot_state = project.procpramarray[ci].dot_state

    if (not ptr_valid(dot_state)) then return

    uname = widget_info(event.id, /uname)
    widget_control, event.id, get_value=val

    case uname of
        'sl_glyph_size'  : (*dot_state).dot_glyph_size  = long(val)
        'sl_r0'          : (*dot_state).dot_rzero       = float(val)
        'sl_lmax'        : (*dot_state).dot_lmax        = long(val)
        'sl_b_thr'       : (*dot_state).dot_bval        = long(val)
        'sl_thr'         : (*dot_state).dot_signal_thr  = long(val)
        'sl_fa_thr'      : (*dot_state).dot_anis_thr    = float(val)
        'btn_use_ga'     : (*dot_state).bg_image        = 1-event.select
        ;'com_ga_stats'   : (*dot_state).ga_stat         = 1-event.select
        'txt_diff_sep'   : begin
           temp = float(val)
           (*dot_state).dot_diff_sep = temp
        end
        'txt_diff_time'   : begin
           temp = float(val)
           (*dot_state).dot_diff_time = temp
        end
        'btn_create': begin
            
           mas_dot_init_obj, dot_obj=reco_obj, $
                             lmax=(*dot_state).dot_lmax, $
                             radius=(*dot_state).dot_rzero, $
                             bval_thr=(*dot_state).dot_bval, $
                             diff_time=(*dot_state).dot_diff_time, $
                             diff_sep=(*dot_state).dot_diff_sep, $
                             /bulk

            if ((*difftools_state).disp_3d) then begin

                case project.procpramarray[ci].slice_axis of
                    0: slice = project.procpramarray[ci].sdim_start
                    1: slice = project.procpramarray[ci].pdim_start
                    2: slice = project.procpramarray[ci].fdim_start
                endcase

                mdt_make_glyph3d, reco_obj, [ slice ], threshold=(*dot_state).dot_signal_thr

            endif else begin

               if ((*dot_state).bg_image eq 0) then begin
                  ga = mas_dot_compute_ga(reco_obj)
               endif

               mdt_make_glyphimage, reco_obj, $
                                    image_in=ga, $
                                    anis_threshold=(*dot_state).dot_anis_thr, $
                                    threshold=(*dot_state).dot_signal_thr, $
                                    /display
              endelse  
         end

        'btn_batch': begin

            slices = mas_get_multi_slice_selection(event.top, count=count)
            if (count eq 0) then return

            mas_dot_init_obj, dot_obj=reco_obj, $
                              lmax=(*dot_state).dot_lmax, $
                              radius=(*dot_state).dot_rzero, $
                              bval_thr=(*dot_state).dot_bval, $
                              diff_time=(*dot_state).dot_diff_time, $
                              diff_sep=(*dot_state).dot_diff_sep, $
                              /bulk

            if ((*difftools_state).disp_3d) then begin
               mdt_make_glyph3d, reco_obj, [ slices ], threshold=(*dot_state).dot_signal_thr
                
            endif else begin

                for sl = 0, n_elements(slices)-1 do begin
                   slice=slices[sl]
                   if ((*dot_state).bg_image eq 0) then begin
                      ga = mas_dot_compute_ga(reco_obj, vis_slice=slice)
                   endif
                   
                   mdt_make_glyphimage, reco_obj, $
                                        image_in=ga, $
                                        vis_slice=slice, $
                                        anis_threshold=(*dot_state).dot_anis_thr, $
                                        threshold=(*dot_state).dot_signal_thr, $
                                        /display
                endfor
            endelse
        end
        'get_ga_vol': begin 
        
        common scan_data
        
        sdim = project.iMNDArray[CI].sdim 
        pdim = project.IMNDArray[ci].pdim
        fdim = project.IMNDArray[ci].fdim
        
        lmax=(*project.procpramarray[ci].dot_state).dot_lmax
        rzero=(*project.procpramarray[ci].dot_state).dot_rzero
        bval=(*project.procpramarray[ci].dot_state).dot_bval
        diff_time=(*project.procpramarray[ci].dot_state).dot_diff_time
        diff_sep=(*project.procpramarray[ci].dot_state).dot_diff_sep
        
        mas_dot_init_obj, dot_obj=dot, $
                          lmax=lmax, $
                          radius=rzero, $
                          bval_thr=bval, $
                          diff_time=diff_time, $
                          diff_sep=diff_sep, $
                          /bulk
        ga_vol = ptr_new(fltarr(fdim,pdim,sdim),/no_copy)  
                            
        (*ga_vol)[*]=0.0
        
        progressbar = obj_new('progressbar', color='Red', $
                           title='Calculating GA...', $
                           text='GA map')
        progressbar->Start
        
        ;; Save the user-selected slice axis as we are going
        ;; to set it to the default to compute the volume
        slice_axis_save = project.procpramarray[ci].slice_axis
        
        catch, error_state
        if (error_state ne 0) then begin
          catch, /cancel
          project.procpramarray[ci].slice_axis = slice_axis_save
          void = dialog_message(['GA computation failed:', $
                                 !ERROR_STATE.msg], /error, /center)
          return
        endif
        
        project.procpramarray[ci].slice_axis = 0
        for i=0, sdim-1 do begin
            (*ga_vol)[*,*,i] = mas_dot_compute_ga(dot, vis_slice=i)
             if (i mod 10 eq 0) then begin
            progressbar->Update, float(i)/sdim  * 100.0
            if (progressbar->CheckCancel()) then begin
                progressbar->Destroy
                
            endif
        endif
        endfor
        progressbar->Destroy
        ;; restore user-selected slice axis
        project.procpramarray[ci].slice_axis = slice_axis_save
        
     if (not keyword_set(file_name)) then begin   
        dest = dialog_pickfile(path=project.current_path, $
                               title='Where would you like to save the GA map?', $
                               filter='*.nii', default_extension='nii', /overwrite_prompt, $
                               /write)
        if (dest eq '') then return
    endif else begin
             CD, CURRENT=c & path_name=c 
             file_name= string(path_name+'/ga.nii')
    endelse
         
        dest_dir = file_dirname(dest)
                  
        if (file_test(dest)) then begin
        if (file_test(dest, /write)) then begin
            void = dialog_message(['File exists:', $
                dest, $
                'Do you want to overwrite?'], /question, /center)
            if (void eq 'No') then return
        endif else begin
            void = dialog_message('Unable to save file, permission denied.', /error, /center)
            return
        endelse
    endif else if (not file_test(dest_dir, /write)) then begin
        void = dialog_message('Unable to save file, permission denied.', /error, /center)
        return
    endif  
        mas_export_nifti, data_ptr=ga_vol, file_name= dest
        end    
                
        'com_ga_stats': begin 
        
        sdim = project.iMNDArray[CI].sdim 
        pdim = project.IMNDArray[ci].pdim
        fdim = project.IMNDArray[ci].fdim
        sdim_Start = project.procPramArray[CI].sdim_Start
        pdim_start = project.procPramArray[ci].pdim_start
        fdim_start = project.procPramArray[ci].fdim_start
      
        Multi_check = project.procPramArray[ci].single_multi_flag
        if Multi_check eq 1 then begin 
          void =  dialog_message('Please change to single slice and make sure the slice is correct',/error)
                    return
                    end
        
        mas_dot_init_obj, dot_obj=reco_obj, $
                              lmax=(*dot_state).dot_lmax, $
                              radius=(*dot_state).dot_rzero, $
                              bval_thr=(*dot_state).dot_bval, $
                              diff_time=(*dot_state).dot_diff_time, $
                              diff_sep=(*dot_state).dot_diff_sep, $
                              /bulk
                         
        if not ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
                    void =  dialog_message('Please either make an ROI or change selection',/error)
                    return
                    end
                    regions = (*project.roi.pROIs[project.roi.ci])
                 
                ;check to make sure the object regions is valid
                if not(total(obj_valid(regions))) then begin
                  update_status_bar, 'No region returned from ROI tool'
                  return
                  end

                ;print, 'size of regions', size(regions)

                  if     (size(regions))[0] eq 0 then numRoi = 1 $
                  else numRoi =  (size(regions))[1]

                  ;make an array to hold the name of the roi's
                  roi_name_array = strarr(numRoi)
    
                  ; make a pointer array to hold all the masks
                  mask = ptrarr(numRoi)
                  dim = [fdim, pdim]
             for temp=0 , numRoi-1 do begin
                if (not obj_valid(regions[temp])) then continue

                mask[temp] = ptr_new(regions[temp]->ComputeMask(dimensions=dim, $
                                                        MASK_RULE=2))

                regions[temp] -> GETPROPERTY, name= name
                roi_name_array[temp] = name
            end      
                
                  ga_roi = ptr_new(fltarr(fdim, pdim))
                  
                
                 (*ga_roi)[*,*] = mas_dot_compute_ga(reco_obj, vis_slice=sdim_Start) 
            
                 stat_array = fltarr(8, numRoi)
                  
                  
             for maskCounter = 0, numRoi-1 do begin

              maskPixels = where( *(mask[maskCounter]) gt 0, count)
              if count eq 0 then continue
             ; mask = mas_roi_get_current_mask()
                  IMAGE_STATISTICS, *ga_roi, mask=*(mask[maskCounter])   $
                               , MEAN=image_mean $
                               , STDDEV=image_stddev $
                               , DATA_SUM=image_sum $
                               , SUM_OF_SQUARES=image_sum_of_squares $
                               , MINIMUM=image_min $
                               , MAXIMUM=image_max $
                               , VARIANCE=image_var 

                stat_array[*, maskCounter] = [ image_mean, image_stddev, image_sum, image_sum_of_squares, image_min, image_max, image_var, count]
                endfor 
              print_format = '(E0, E0, E0, E0, E0, E0, E0, E0, E0)'
              htab         = string(09B)
              crlf         = string(13B)+string(10B)
              header       = 'Slice'+htab+'Region Name'+htab+ $
                    'Mean'+htab+'Std. Dev.'+htab+'Sum'+htab+'Sum of Squares'+htab+ $
                    'Min'+htab+'Max'+htab+'Variance'+htab+'Count'
              output_lines = ['File: '+project.imndarray[ci].file_path, ''] 
              line         = 0
               
              output_lines = [output_lines, header]
         
         for roi = 0, n_elements(roi_name_array)-1 do begin
            
                pref = string(sdim_Start, format='(I0)')+htab+roi_name_array[roi]
                stat = strjoin(string(stat_array[*, roi], format='(E0)'), htab)
                output_lines = [ output_lines, pref+htab+stat ]
                output_lines = [output_lines, '']
            
            output_lines = [output_lines, '']
         endfor 
         
         display_stats, output_lines, 'GA Statistics'         
                  
        end

              else: begin
            print, 'mas_dot_gui_event: recv''d event from unknown widget: '+uname
        end
       
      
    endcase

end

;; Subroutine name: mas_dot_gui_make
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; creates the GUI of DOT-specific options.
;;
;; Editing Information:

pro mas_dot_gui_make, base, state

    common scan_data
    ci = project.ci

    if (ptr_valid(project.procpramarray[ci].dot_state)) then begin
        dot_state = project.procpramarray[ci].dot_state
    endif else begin
        dot_state = mas_dot_gui_make_default_state()
        project.procpramarray[ci].dot_state = dot_state
    endelse
    
    param_base = widget_base(base, column=2, grid=0, xpad=0)
    
    pbase_c1 = widget_base(param_base, /column, /frame)
    sl_r0 = cw_fslider(pbase_c1, title='Radius (um)', $
                       value=(*dot_state).dot_rzero,  $
                       minimum=0, $
                       maximum=20,$
                       scroll=1,  $
                       /edit,     $
                       uname='sl_r0')

    sl_lmax = widget_slider(pbase_c1, title='Harmonic Order (L)', $
                            value=(*dot_state).dot_lmax,  $
                            minimum=2, $
                            maximum=8,$
                            scroll=1,  $
                            event_pro='mas_dot_gui_event', $
                            uname='sl_lmax')

    sl_b_thr = widget_slider(pbase_c1, title='Low B-value Threshold', $
                             minimum=100,  $
                             maximum=1000, $
                             value=(*dot_state).dot_bval,    $
                             scroll=100,   $
                             event_pro='mas_dot_gui_event', $
                             uname='sl_b_thr')
    
    btn_base = widget_base(pbase_c1, /column)
    com_ga_stats = widget_button(btn_base, $
                               value='Compute GA stats', $
                               event_pro='mas_dot_gui_event', $
                               uname='com_ga_stats')
                               btn_base = widget_base(pbase_c1, /column)
   
    get_ga_vol = widget_button(btn_base, $
                               value='Save GA map', $
                               event_pro='mas_dot_gui_event', $
                               uname='get_ga_vol')
                               
    non_ex_base = widget_base(pbase_c1, /nonexclusive, /align_center)
    
    btn_use_ga = widget_button(non_ex_base, value="Use GA Background", uname='btn_use_ga')
    
    widget_control, btn_use_ga, set_button=1-project.procpramarray[ci].adt_proccess_flag
    (*dot_state).bg_image = project.procpramarray[ci].adt_proccess_flag

    pbase_c2 = widget_base(param_base, /column, /frame)
    sl_dot_thr = widget_slider(pbase_c2, MINIMUM=0, MAXIMUM=100, $
                               event_pro='mas_dot_gui_event', $
                               UNAME='sl_thr', TITLE='% Signal Threshold', $
                               VALUE=(*dot_state).dot_signal_thr)

    sl_fa_thr = cw_fslider(pbase_c2, title='FA Threshold', $
                           value=(*dot_state).dot_anis_thr,  $
                           minimum=0.001, $
                           maximum=1.0,$
                           scroll=0.1,  $
                           /edit,     $
                           uname='sl_fa_thr')
    widget_control, sl_fa_thr, sensitive=project.procpramarray[ci].adt_proccess_flag

        
    dt_base = widget_base(pbase_c2, /row)
    lbl_dt = widget_label(dt_base, value="Diffusion Time:")
    txt_diff_time = widget_text(dt_base, xsize=7, $
                                uname='txt_diff_time', $
                                event_pro='mas_dot_gui_event', $
                                /editable, /all_events, $
                                value=strcompress((*dot_state).dot_diff_time, /remove_all))

    ds_base = widget_base(pbase_c2, /row)
    lbl_ds = widget_label(ds_base, value="Diff Pulse Length:")
    txt_diff_sep = widget_text(ds_base, xsize=7, $
                               uname='txt_diff_sep', $
                               event_pro='mas_dot_gui_event', $
                               /editable, /all_events, $
                               value=strcompress((*dot_state).dot_diff_sep, /remove_all))
    
    widget_state = ptr_new({ DOT_WIDGET_STATE, $
                             btn_use_ga: btn_use_ga, $
                             txt_diff_sep: txt_diff_sep, $
                             txt_diff_time: txt_diff_time, $
                             sl_r0: sl_r0, $
                             sl_lmax: sl_lmax, $
                             sl_b_thr: sl_b_thr, $
                             sl_fa_thr: sl_fa_thr, $
                             com_ga_stats: com_ga_stats, $
                             sl_dot_thr: sl_dot_thr }, /no_copy)

    widget_control, base, set_uvalue=widget_state

end

;; Subroutine name: mas_dot_redraw
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Redraws the GUI, triggered by mas control
;;
;; Editing Information:

pro mas_dot_redraw, dot_base

    common scan_data, project

    ci = project.ci

    dot_state = project.procpramarray[ci].dot_state
    
    if (not ptr_valid(dot_state)) then begin
        
        if (project.ni eq 0) then begin

            sen = 0
            widget_control, dot_base, get_uvalue=widget_state
            
            widget_control, (*widget_state).sl_r0     , sensitive=sen
            widget_control, (*widget_state).sl_lmax   , sensitive=sen
            widget_control, (*widget_state).sl_b_thr  , sensitive=sen
            widget_control, (*widget_state).sl_dot_thr, sensitive=sen
            widget_control, (*widget_state).btn_use_ga, sensitive=sen
            widget_control, (*widget_state).sl_fa_thr , sensitive=sen
            
            return

        endif else begin

            dot_state = mas_dot_gui_make_default_state()
            project.procpramarray[ci].dot_state = dot_state

        endelse

    endif

    sen = 1
    widget_control, dot_base, get_uvalue=widget_state

    have_adt = project.procpramarray[ci].adt_proccess_flag
    if not have_adt then (*dot_state).bg_image = 0

    widget_control, (*widget_state).sl_r0     , sensitive=sen, set_value=(*dot_state).dot_rzero
    widget_control, (*widget_state).sl_lmax   , sensitive=sen, set_value=(*dot_state).dot_lmax
    widget_control, (*widget_state).sl_b_thr  , sensitive=sen, set_value=(*dot_state).dot_bval
    widget_control, (*widget_state).sl_dot_thr, sensitive=sen, set_value=(*dot_state).dot_signal_thr
    widget_control, (*widget_state).sl_fa_thr , sensitive=sen*have_adt, set_value=(*dot_state).dot_anis_thr
    widget_control, (*widget_state).btn_use_ga, sensitive=sen, set_button=(1-(*dot_state).bg_image)

end

;; Subroutine name: mas_dot_gui_cleanup
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; cleans up when the user closes the DOT window. Note that
;; only the widget state pointer gets destroyed. The dot state
;; pointer stays alive for the next time the user opens the dot
;; gui. It is destroyed when a scan is removed from mas.
;;
;; Editing Information:

pro mas_dot_gui_cleanup, dot_base

    common scan_data

    project.dot_tlb = -1
    
    widget_control, dot_base, get_uvalue=widget_state
    
    ptr_free, widget_state

end

;; Subroutine name: mas_dot_init_obj
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    LMAX:      The highest order of the harmonic basis.
;;
;;    BULK:      Set this keyword to get a "BULK_RECONSTRUCTOR" version of
;;               the object
;;
;;    DIFF_TIME: The diffusion time from the pulse sequence parameters
;;               in MS (big delta)
;;
;;    DIFF_SEP:  The Diffusion pulse length in MS (small delta)
;;
;;    RADIUS:    The diffusion radius to be measured
;;
;;   DOT_OBJ:    Set this keyword to a named variable that will contain
;;               the prepared object.
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Creates an instance of the DOT reconstructor object using imaging
;; parameters from MAS's data structures.
;;
;; Editing Information:

pro mas_dot_init_obj, lmax=lmax, $
                      diff_time=diff_time, $
                      diff_sep=diff_sep, $
                      radius=radius, $
                      dot_obj=dot_obj, $
                      bval_thr=bval_thr, $
                      bulk=bulk

    common scan_data

    ci = project.ci

    ;;; IMAGING PARAMETERS ;;;;
    ;; convert to s/m^2, appears to be what Evren uses
    bvals = reform(*project.imndarray[ci].bval_array) * 1e-6
    
    ;; get the gradient directions
    gradients = fltarr(n_elements(bvals), 2)
    gradients[*,0] = (*project.imndarray[ci].angle_theta)
    gradients[*,1] = (*project.imndarray[ci].angle_phi)
    gradients = transpose(gradients)

    if (not keyword_set(bval_thr)) then begin
        bval_thr = 400.
    endif
    
    if (keyword_set(lmax)) then begin
       lmax = ( lmax > 1 ) < 10
    endif else begin
       lmax = 4
    endelse

    ;; diffusion time parameters
    if (keyword_set(diff_time)) then begin
       big_delta = float(diff_time)/1000.
    endif else begin
       big_delta = ((*project.imndArray[CI].big_delta)[0]) / 1000 
    endelse

    if (keyword_set(diff_sep)) then begin
       small_delta = float(diff_sep)/1000.
    endif else begin
       small_delta = ((*project.imndArray[CI].small_delta)[0]) / 1000 
    endelse

    td = big_delta - small_delta/3.
   
    if (keyword_set(bulk)) then begin

       dot_obj = obj_new('mas_dot_bulk_reconstructor', $
                         diff_time=td, $
                         radius=radius, $
                         b_values=bvals, $
                         bval_thr=bval_thr, $
                         gradients=gradients, $
                         l_max=lmax, $
                         /spherical, /degrees)

    endif else begin

       dot_obj = obj_new('mas_dot', $
                         diff_time=td, $
                         radius=radius, $
                         b_values=bvals, $
                         bval_thr=bval_thr, $
                         gradients=gradients, $
                         l_max=lmax, $
                         /spherical, /degrees)

    endelse
    

 end

pro mas_dot_tools

end
