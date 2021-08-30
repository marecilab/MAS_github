;; $Id$
;; Copyright 2008 University of Florida. All Rights Reserved

;; Subroutine name: mas_3d_slicer
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    pMatrix:  pointer to a 3-d image volume. Note that this 
;;              procedure doesn't free the pointer
;;    
;;    supplemental_omodel: An extra IDLGrModel containing some 
;;                         graphics objects to display along with
;;                         the image matrix.
;;
;;    {r,p,s}_scale: voxel scaling factors for each dimension. 
;;                   mas_3d_slicer will zoom each dimension accordingly. 
;;
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Allows user to slice up a 3d volume along the 3 orthoginal axes. 
;;
;; Editing Information:

pro mas_3d_slicer, $
                  pMatrix, $
                  supplemental_omodel=supplemental_omodel, $
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

    ;; Set up the various models
    omodel          = obj_new('idlgrmodel', name='MAIN_model', depth_test_disable=2, lighting=0)

    omodel_PS       = obj_new('idlgrmodel', name="PS_image", depth_test_disable=2, lighting=1)
    omodel_FP       = obj_new('idlgrmodel', name="FP_image", depth_test_disable=2, lighting=1)
    omodel_FS       = obj_new('idlgrmodel', name="FS_image", depth_test_disable=2, lighting=1)
    omodel_im       = obj_new('idlgrmodel', name='IMG_model', depth_test_disable=2, lighting=1)

    omodel_axes     = obj_new('idlgrmodel', name='AXES_model', lighting=0)

    if (not keyword_set(f_scale)) then f_scale = 1
    if (not keyword_set(p_scale)) then p_scale = 1
    if (not keyword_set(s_scale)) then s_scale = 1
    omodel_PS->scale, p_scale, s_scale, 1
    omodel_FS->scale, f_scale, s_scale, 1
    omodel_FP->scale, p_scale, f_scale, 1
    
    ;; order is important
    omodel_marker = obj_new('idlgrmodel', name='MARKER_model', depth_test_disable=2, lighting=1)
    omodel->add, omodel_marker

    omodel->add, omodel_im;;,       position=0    
    omodel->add, omodel_axes;;,     position=1
    
    omodel_im->add, omodel_PS
    omodel_im->add, omodel_FP
    omodel_im->add, omodel_FS
    
    omodel_FS->rotate, [1,0,0], 90
    omodel_PS->rotate, [1,0,0], 90
    omodel_PS->rotate, [0,0,1], 90

    if keyword_set(supplemental_omodel) and obj_valid(supplemental_omodel) then begin
        omodel->add, supplemental_omodel;;, position=2
    endif
    
    ;; rotate the image model as necessary
;;     case rotate_direction of
;;         0: deg = 0
;;         1: begin
;;             deg = 90
;;             omodel_im->translate, -matrix[1]/2, -matrix[0]/2, 0
;;             omodel_im->rotate, [0,0,1], deg
;;             omodel_im->translate, matrix[0]/2, matrix[1]/2, 0
;;          end

;;         2: begin
;;             deg = 180
;;             omodel_im->translate, -matrix[0]/2, -matrix[1]/2, 0
;;             omodel_im->rotate, [0,0,1], deg
;;             omodel_im->translate, matrix[0]/2, matrix[1]/2, 0
;;         end

;;         3: begin
;;             deg = 270
;;             omodel_im->translate, -matrix[1]/2, -matrix[0]/2, 0
;;             omodel_im->rotate, [0,0,1], deg
;;             omodel_im->translate, matrix[0]/2, matrix[1]/2, 0
;;         end

;;     endcase

;;  model             - the main idlgrmodel object. contains
;;                      everything else
;;  tlb               - the tlb for the main viewport. returned by
;;                      xobjview
;;  matrix            - the dimensions (from mas) of the environment
;;  active_FP         - active freq-phase slice number
;;  active_FS         - active freq-slice slice number
;;  active_PS         - active phase-slice slice number
;;  f_scale           - amount to scale freq direction based on voxel
;;                      dim.
;;  p_scale           - amount to scale phase direction based on voxel
;;                      dim.
;;  s_scale           - amount to scale slice direction based on voxel
;;                      dim.
;;  image_alpha       - alpha percent for image display
;;  image_contrast    - contrast for image display
;;  image_datapool    - pointer to the data array that contains the
;;                      image data. don't free this pointer!
;;  marker_list       - the pointer to an array of idlgrpolyline
;;                      objects that are the markers in 3D space
;;  marker_list_wid   - the widget id of the marker list. used by
;;                      add/remove/update routines
;;  min/max intensity - the intensity values used to adjust contrast
    marker_list = project.procpramarray[project.ci].marker_list

    state = ptr_new({ model:omodel, $
                      tlb:0, $
                      matrix:matrix, $
                      active_FP:-1, $
                      active_FS:-1, $
                      active_PS:-1, $
                      f_scale: f_scale, $
                      p_scale: p_scale, $
                      s_scale: s_scale, $
                      image_alpha: .75, $
                      image_contrast: 0, $
                      image_datapool: ptr_new(), $
                      marker_list: ptr_new(), $
                      marker_list_wid: 0L, $
                      min_intensity: long(0), $
                      max_intensity: long(0)})

    ;; set up axis - always performed.
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

    xobjview, omodel, background=[0,0,0], tlb=tlb, xsize=800, ysize=600, renderer=1

    (*state).tlb = tlb
    (*state).image_datapool = pMatrix
    (*state).min_intensity = min(*pMatrix)
    (*state).max_intensity = max(*pMatrix)

    mas_3d_slicer_GUI, state=state

end

;; Subroutine name: mas_3d_slicer_test
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;   matrix: Set this to a named variable that will contain the example matrix.
;;           This is necessary so that the user may free the pointer, since the
;;           procedure returns immediately.
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Exmaple procedure to demonstrate mas_3d_slicer
;;
;; Editing Information:

pro mas_3d_slicer_test, matrix=matrix

    matrix = dist(8000,8000)
    matrix = reform(matrix[*], 400,400,400)
    matrix = bytscl(matrix)
    mas_3d_slicer, ptr_new(matrix, /no_copy)

end
