;; $Id$
;; Copyright 2015 University of Florida. All Rights Reserved

;; Subroutine name: mas_3d_surface
;; Created by: Kevin Montes, 2015-01
;; Calling Information:
;;
;;    pMatrix:  pointer to a 3-d image volume. Note that this 
;;              procedure doesn't free the pointer
;;
;;    {f,p,s}_scale: voxel scaling factors for each dimension. 
;;                   mas_3d_surface will zoom each dimension accordingly. 
;;
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Allows user to display a surface rendering of a 3D volume 
;;
;; Editing Information:

pro mas_3d_surface, $
                  pMatrix, $
                  f_scale=f_scale, $
                  p_scale=p_scale, $
                  s_scale=s_scale

    common scan_data, project
    
    if (not ptr_valid(pMatrix)) then begin
        junk = dialog_message('Unable to view volume data.', /center, /error)
        return
    endif

    sz = size(*pMatrix)
    matrix = [sz[1], sz[2], sz[3]]
    
; Get minimum and maximum isovalues possible, then compute initial isovalue guesses 1 and 2    
    min_isoval = min(*pMatrix)
    max_isoval = max(*pMatrix)
    isoval_1 = (max_isoval-min_isoval)/4
    isoval_2 = 3*isoval_1
    
; Create main models to store 3D data for viewing     
    omodel = obj_new('idlgrmodel', name='MAIN_model')
    omodel_im = obj_new('idlgrmodel', name='IMG_model')
    omodel_image = obj_new('idlgrmodel', name='Image_object')

    start_mem = MEMORY(/CURRENT)  
;Create tetrahedral mesh for surface rendering and store in new omodel object   
    INTERVAL_VOLUME, *pMatrix, isoval_1, isoval_2, verts, conn
    conn = TETRA_SURFACE(verts, conn)
    omodel_surface = obj_new('IdlgrPolygon', verts, POLYGONS=conn, $
      COLOR=[200,200,200], SHADING=1, name='Surface_image')
 PRINT, 'Memory required: ', (MEMORY(/HIGHWATER) - start_mem)/(1024.0)^3
 
;Objects for axes, markers, and images created and added to MAIN model  
    omodel_axes = obj_new('idlgrmodel', name='AXES_model', lighting=0)
    
    omodel->add, omodel_im
    omodel->add, omodel_axes
    
    omodel_image->add, omodel_surface
    omodel_image->scale, f_scale, p_scale, s_scale 
    omodel_im->add, omodel_image
    
    if keyword_set(supplemental_omodel) and obj_valid(supplemental_omodel) then begin
        omodel->add, supplemental_omodel;;, position=2
    endif

;Store variables in state structure for GUI
    
    state = ptr_new({ model:omodel, $
                      tlb:0, $
                      matrix:matrix, $
                      f_scale: f_scale, $
                      p_scale: p_scale, $
                      s_scale: s_scale, $
                      image_alpha: .75, $
                      image_contrast: 0, $
                      image_datapool: ptr_new(), $
                      min_isoval: min_isoval, $
                      max_isoval: max_isoval, $
                      min_intensity: long(0), $
                      max_intensity: long(0)})

;always performed, set up axes   
    if 1 or keyword_set(show_axis) then begin
    
        oaxfnt = obj_new('idlgrfont', size=4)
        oxaxis = obj_new('idlgraxis', 0, range=[0,matrix[0]], color=[200,200,200])
        oyaxis = obj_new('idlgraxis', 1, range=[0,matrix[1]], color=[200,200,200])
        ozaxis = obj_new('idlgraxis', 2, range=[0,matrix[2]], color=[200,200,200])
        
        oxaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt
        oyaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt
        ozaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt

        omodel_axes->add, oxaxis
        omodel_axes->add, oyaxis
        omodel_axes->add, ozaxis

    endif

; Display omodel in a graphical window
    
    xobjview, omodel, BACKGROUND=[0,0,0], tlb=tlb, xsize=800, ysize=600, renderer=1
    
    (*state).tlb = tlb
    (*state).image_datapool = pMatrix
    (*state).min_intensity = min(*pMatrix)
    (*state).max_intensity = max(*pMatrix)

    mas_3d_surface_GUI, state=state
;    plot, histogram(*pMatrix)
    
end