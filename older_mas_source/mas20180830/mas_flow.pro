; PC-MRI data processing routine
; ;
; Subroutine name: calc_venc
; Created by: Magdoom Kulam
; Calling Information:
;
;   G      - Maximum flow-encoding gradient strength (G/cm)
;   t      - Half-duration of the flow encoding gradient (us)
;   bpt    - Flow encoding gradient shape, trapezoidal by default (Other : Sinusoid, Gauss, Shifted cosine)
;   scheme - Flow encoding scheme, Hadamard by default (Other : Four point, Six point)
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To calculate the velocity encode value (mm/s)
;
; Editing Information:

function calc_venc, G, t, ncycl, bpt, scheme

common scan_data

CI = project.ci
gamma_bar = 4258E-6                     ; Gyromagnetic ratio for H1 (MHz/G)

; Spin-echo diffusion sequence
if strcmp(project.imndArray[CI].diff,'y') eq 1 then begin
  tbig_delta = 1e3**project.imndArray[CI].big_delta
  tsmall_delta = 1e3**project.imndArray[CI].small_delta
  return, 1E7/(2*4.0*gamma_bar*G*tbig_delta*tsmall_delta)
endif

if ptr_valid(project.imndArray[CI].big_delta) then tbig_delta = 1e3**project.imndArray[CI].big_delta else tbig_delta = t
 
; Scaling factor based on the encoding scheme
case scheme of
  'Four Point' : eta = 1.0
  'Six point'  : eta = 2.0
   else        : eta = 4.0
endcase

case bpt of
'Sinusoid'       :  venc = !pi*1E7/(4*gamma_bar*G*t*tbig_delta)          ; Velocity encode value for sine (mm/s)
'ShiftedCosine'  :  venc = 1E7/(gamma_bar*G*t*tbig_delta)                ; Velocity encode value for shifted cosine (mm/s)
'Gauss'          :  venc =  9.6*1E7/(2*gamma_bar*G*t*tbig_delta)         ; Velocity encode value for Gauss (mm/s)
else             :  begin
                     max_slew_rate = project.imndarray[CI].max_slew_rate ; Maximum slew rate (G/cm.us)
                     trise = G/max_slew_rate                             ; Rise time for the bipolar trapezoidal gradient (us)
                     ; Corrected rise time to account for 4 us delay
                     if ((trise mod 4) eq 0) then begin
                        trise_corr = trise
                     endif else begin
                        trise_corr = 4*(1+floor(trise/4))                                        
                     endelse
                     tplateau = t - 2*trise_corr                                  ; Plateau time for the bipolar trapezoidal gradient (us)      
                     venc = 1E7/(2*gamma_bar*G*(tplateau+trise_corr)*tbig_delta)  ; Velocity encode value for trapezoid in mm/s
                     end
endcase
return, venc/(eta*ncycl)
end

; Subroutine name: log2physgradtransform
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Returns transformation matrix from logical to physical gradient axes
;
; Editing Information:
function log2physgradtransform

  common scan_data

  theta = (*project.imndArray[project.ci].theta)[0]
  phi = (*project.imndArray[project.ci].phi)[0]
  psi = (*project.imndArray[project.ci].psi)[0]

  theta*= !Pi/float(180)
  phi*= !Pi/float(180)
  psi*= !Pi/float(180)

  M = [[sin(phi)*cos(psi)-cos(phi)*cos(theta)*sin(psi), -cos(phi)*cos(psi)-sin(psi)*cos(theta)*sin(phi), sin(theta)*sin(psi)],$
    [-sin(phi)*sin(psi)-cos(phi)*cos(theta)*cos(psi), cos(phi)*sin(psi)-sin(phi)*cos(theta)*cos(psi), sin(theta)*cos(psi)], $
    [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]]
    
    return, M
 end
 
; Subroutine name: load_state2_data_flow
; Created by: Magdoom Kulam
; Calling Information:
;
; data - Pointer containing the processed state 1 flow data
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To prepare the velocity vector for display after rotation/flip/zoom/thresholding
;
; Editing Information:

function load_state2_data_flow, data

common scan_data

get_flow_params, gflow, tflow, ncycl, shape, fc, scheme, venc

if ~ptr_valid(data) then return, !NULL
 
xzoom = project.procpramArray[project.ci].x_zoom  ; x-zoom factor in MAS main widget
yzoom = project.procpramArray[project.ci].y_zoom  ; y-zoom factor in MAS main widget

N = size(*data,/Dimension)
temp = *data
if ptr_valid(V) then ptr_free,V
V = ptr_new(make_array(xzoom*N[0],yzoom*N[1],N[2],5,value = 0,/double))
 
for i = 0,4 do begin 
  (*V)[*,*,*,i] = rebin(temp[*,*,*,i],xzoom*N[0],yzoom*N[1],N[2],/sample)   ; Zooming
   if i ne 4 then (*V)[*,*,*,i]*=venc                                            ; Velocity scaling
endfor

mask_flow_data,V
orient_flow_data,V
rotate_flow_data, V
flip_flow_data,V

return, V

end

; Subroutine name: get_flow_params
; Created by: Magdoom Kulam
; Calling Information:
; 
; gflow  - Strength of the flow encoding gradient (mT/m)
; tflow  - Duration of the flow encoding gradient (ms)
; shape  - Shape of the flow encoding gradient, trapezoid by default (Other : Sinusoid, Gauss, Shifted cosine)
; fc     - Flow compensation (ON or OFF)
; scheme - Flow encoding scheme, Hadamard by default (Other : Four point, Six point)
; venc   - Velocity encode value (mm/s)
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To retrieve the flow related parameters for the selected scan
;
; Editing Information:

pro get_flow_params, gflow, tflow, ncycl, shape, fc, scheme, venc

common scan_data

CI = project.ci 
gflow = 10*project.imndarray[CI].gflow     
tflow = 2*1e3*project.imndarray[CI].tflow 
ncycl = project.imndArray[CI].ncycl         ; Number of cycles of flow gradient

shape = project.imndarray[CI].bptype
case shape of
's'   : shape = 'Sinusoid'     
'c'   : shape = 'ShiftedCosine' 
'g'   : shape = 'Gauss'
else  : shape = 'Trapezoid'     
endcase

fc = project.imndarray[CI].fc
case fc of
'y' :  fc = 'ON'          
'n' :  fc = 'OFF'
else : fc = '--'
endcase 

scheme = project.imndarray[CI].scheme
case scheme of
  'f'   : scheme = 'Four Point'
  's'   : scheme = 'Six point'
  'h'  : scheme = 'Hadamard'
  else : scheme = 'Test'
endcase
venc = calc_venc(float(gflow)/float(10),float(1000)*float(tflow)/float(2),ncycl,shape,scheme) 
venc*=1e3

end

; Subroutine name: flow_show_gradient_cleanup
; Created by:
; Calling Information:
;
;   tlb: top-level base
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Cleans up pointers and objects associated with the gradient direction viewer
;
; Editing Information:

pro flow_show_gradient_cleanup, tlb

  widget_control, tlb, get_uvalue=uvalue

  obj_destroy, (*uvalue).otxtmodel
  obj_destroy, (*uvalue).overts
  obj_destroy, (*uvalue).osurf
  obj_destroy, (*uvalue).pts
  ptr_free, uvalue

end

; Subroutine name: flow_show_gradient_angle_event
; Created by:
; Calling Information:
;
;   ev: the button press event
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Handles click of "Hide/Show Angles" button on the 3D viewer top
;    bar.
;
; Editing Information:

pro flow_show_gradient_angle_event, ev

  widget_control, ev.id, get_uvalue=uvalue
  (*uvalue).otxtmodel->getProperty, hide=isHidden
  (*uvalue).otxtmodel->setProperty, hide=(1 - isHidden)

  if (isHidden eq 0) then begin
    widget_control, ev.id, set_value='Show Angles'
  endif else begin
    widget_control, ev.id, set_value='Hide Angles'
  endelse

  xobjview, refresh=ev.top

end

; Subroutine name: flow_show_gradient_profile
; Created by:
; Calling Information:
;
;   GRADIENTS - user-supplied gradients in substitute for the
;               current mas scan's gradients
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To create display flow gradient profile
;
; Editing Information:

pro flow_show_gradient_profile, gradients=gradients

  common scan_data
  if ~ptr_valid(project.imndarray[project.ci].fro) && ~keyword_set(gradients) then begin
    mbox = dialog_message('Flow gradient direction not specified',/error)
    return
  endif
  
  if keyword_set(gradients) then begin

    gro = ptr_new(gradients[0,*])
    gpe = ptr_new(gradients[1,*])
    gsl = ptr_new(gradients[2,*])

  endif else begin

    gro = project.imndarray[project.ci].fro
    gpe = project.imndarray[project.ci].fpe
    gsl = project.imndarray[project.ci].fsl
  
    if not (ptr_valid(gro) and ptr_valid(gpe) and ptr_valid(gsl)) then begin
      junk = dialog_message(['Gradient profile can only be shown for flow scans.'], $
        /error, /center)
      return
    endif

  endelse

  grad = fltarr(3, n_elements(*gro))
  grad = [*gro,*gpe,*gsl]  
  grad_sphere = cv_coord(from_rect=grad, /to_sphere, /degrees)

;  catch, qhull_error1
;  if (qhull_error1 ne 0) then begin
;    catch, /cancel
;    void = dialog_message(['Unable to create visualization.', $
;      'There must be at least three unique gradient directions.'], $
;      /error, /center)
;    return
;  endif

  qhull, grad, tetrahedra
  newconn=tetra_surface(grad, tetrahedra)
  osurf = obj_new('IDLgrPolygon', grad, polygons=newconn, style=1, Color=[255,0,0])

  catch, /cancel

  catch, qhull_error
  if (qhull_error ne 0) then begin
    catch, /cancel
    goto, SKIP_QHULL
  endif
  
  SKIP_QHULL:
  overts = objarr(n_elements(*gro))

  pts = obj_new('idlgrpolygon', grad, thick=4, style=0, depth_test_disable=2)
  otxt_model = obj_new('idlgrmodel')
  ofnt = obj_new('idlgrfont', size=8)

  for i = 0, n_elements(*gro)-1 do begin
    overts[i] = obj_new('idlgrpolyline', [ [0,0,0], [grad[*,i]] ], thick=1, $
      name=strcompress('('+string(grad[0,i])+','+string(grad[1,i])+')', /remove_all), $
      linestyle=0, color=[90,90,90], alpha_channel=0.2, depth_test_disable=2)

    otxt_model->add, obj_new('idlgrtext', location=[grad[*,i]+0.05], color=[0,0,0], /onglass, font=ofnt,$
      alpha_channel=0.6, $
      strcompress(string(grad_sphere[1,i], format='(3F0.3)')+','+$
      string(grad_sphere[0,i], format='(3F0.3)'), /remove_all))

  endfor

  omodel = obj_new('idlgrmodel', depth_test_disable=2)

  omodel->add, osurf
  omodel->add, overts
  omodel->add, pts
  omodel->add, otxt_model

  xobjview, omodel, xsize=600, ysize=600, tlb=tlb

  tb = widget_info(tlb, find_by_uname="xobjview:_resetbase")

  if (widget_info(tb, /valid_id)) then begin

    btn = widget_button(tb, $
      value='Hide Angles', $
      uname='ANGLE_BTN', $
      event_pro='flow_show_gradient_angle_event', $
      uvalue=ptr_new({ otxtmodel: otxt_model, $
      overts: overts, $
      osurf: osurf, $
      pts: pts }) )
    xmanager, 'flow_show_gradient_angle', btn, cleanup='flow_show_gradient_cleanup', /no_block

  endif

end

; Subroutine name: meshgrid
; Created by: Magdoom Kulam
; Calling Information:
;
;  x,y,z - Input vectors/output matrix
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To generate a Cartesian grid using vectors x,y and/or z
;
; Editing Information:

pro meshgrid, x,y, z = z

Nx = size(x,/Dimensions)
Ny = size(y,/Dimensions)

x = rebin(x,Nx,Ny,/sample)
y = transpose(rebin(y,Ny,Nx,/sample))

if keyword_set(z) then begin  
    Nz = (size(z))[1]
    type = (size(z))[2]
    x = rebin(x,Nx,Ny,Nz, /sample)
    y = rebin(y,Nx,Ny,Nz, /sample) 
    temp = make_array(Nx,Ny,Nz,value = 0, type = type)
    for i = 0, Nz[0]-1 do temp[*,*,i] = make_array(Nx,Ny,value = z[i], type = type)
    z = reform(temp)
endif

end

; Subroutine name: rotate_3d_vector_field
; Created by: Magdoom Kulam
; Calling Information:
;
;  X      - Pointer containing the 3D vector field (fourth dimension holding the components of the spatially varying vector)
;  theta  - Angle in radians to rotate
;  ind    - Optional keyword to specify (index) where the first vector component is located in the fourth dimension 
;           of the pointer (Zero by default). The other two components are assumed to be in consecutive indices  
;  axis   - Specifies the rotation axis ('x', 'y' or 'z')
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To rotate a 3d vector field by an angle theta CCW about x,y or z (default)-axis 
;
; Editing Information:

pro rotate_3d_vector_field, X, theta, ind = ind, axis = axis

if ~keyword_set(ind) then ind = 0
; Rotation matrices
if keyword_set(axis) then begin
  case axis of
    'x' : R = [[1.0, 0.0, 0.0],[0.0, cos(theta), -sin(theta)],[0,sin(theta), cos(theta)]] 
    'y' : R = [[cos(theta), 0.0, sin(theta)],[0.0, 1.0, 0.0],[-sin(theta), 0.0, cos(theta)]] 
    'z' : R = [[cos(theta), -sin(theta), 0],[sin(theta), cos(theta), 0],[0.0, 0.0, 1.0]]  
  endcase
endif else R = [[cos(theta), -sin(theta), 0],[sin(theta), cos(theta), 0],[0, 0, 1]]  

; Rotating the vector
N = size(*X,/Dimensions)
Xraster = make_array(N[0]*N[1]*N[2],3)                        
for i = 0,2 do Xraster[*,i] = reform((*X)[*,*,*,i+ind],N[0]*N[1]*N[2],1)  ; Rastering the velocity field (3xRxPxS) into 2D (3xRPS) matrix
Xrot = R##Xraster                                                         ; Rotating the rastered array
for i = 0,2 do (*X)[*,*,*,i+ind] = reform(Xrot[*,i],N[0],N[1],N[2])       ; Reforming the rotated array back to 3D

end

; Subroutine name: flow_gui
; Created by: Magdoom Kulam
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Main flow GUI
;
; Editing Information:

pro flow_gui

common scan_data
common common_widgets

if (strcmp(project.imndarray[project.ci].flow,'y') and project.imndarray[project.ci].adim ge 4) then get_flow_params, gflow, tflow, ncycl, shape, fc, scheme, venc $
  else begin
   mbox = dialog_message('Select a phase contrast flow dataset',/error)
   return
endelse 
if ptr_valid(project.imndArray[project.ci].big_delta) then tbig_delta = *project.imndArray[project.ci].big_delta else tbig_delta = tflow/2.0

if (xregistered('flow_gui') ne 0) then return
topbase = widget_base(title = 'Flow Imaging : Scan #' + strtrim(project.ci+1,2) ,tlb_frame_attr=1,/grid_layout,/kbrd_focus_events)
mainbase = widget_base(topbase,column=3)

labelbp = widget_label(mainbase, value = 'Flow Imaging Parameters')
fbase = widget_base(mainbase,row = 2)
parbase = widget_base(fbase, row = 8, frame = 1)
gfield =    cw_field(parbase, title =   'Max. Strength     ', value = strtrim(string(gflow,format = '(F8.2)'),2)+' mT/m',xsize = 12, /noedit)
tfield1 =    cw_field(parbase, title =  'Small Delta       ', value = strtrim(string(tflow/2.0,format = '(F8.2)'),2)+ ' ms',xsize = 12,/noedit)
tfield2 =    cw_field(parbase, title =  'Big Delta         ', value = strtrim(string(tbig_delta,format = '(F8.2)'),2)+ ' ms',xsize = 12,/noedit)
ncyclfield = cw_field(parbase, title =  'No. of Cycles     ', value = ncycl, xsize = 12,/noedit)
shpfield =  cw_field(parbase, title =   'Shape             ', value = shape, xsize = 12,/noedit)
schfield =  cw_field(parbase, title =   'Encoding scheme   ', value = scheme, xsize = 12,/noedit)
fcfield =   cw_field(parbase, title =   'Flow Compensation ', value = fc, xsize = 12, /noedit)
vencfield = cw_field(parbase, title =   'Velocity Encode   ', value =  strtrim(string(venc,format = '(F10.2)'),2) + ' um/s',xsize = 12, /noedit)

ROIbase = widget_base(fbase,row = 2)
ROIlabel = widget_label(ROIbase,value = 'ROI tools')
ROIbtnbase = widget_base(ROIbase,column = 2,frame = 1)
ROIprintbtn = widget_button(ROIbtnbase, value = 'Print Stats', uvalue = 'ROIPrint',/no_release)
ROIhistbtn = widget_button(ROIbtnbase, value = 'Plot Histogram', uvalue = 'ROIHist',/no_release)
ROIProp_buttons = ['Flow Rate','Flow Anisotropy','Flow Velocity']
ROI_button_grp = cw_bgroup(ROIbtnbase, ROIProp_buttons, uvalue = 'ROI_buttons', /nonexclusive, column = 1, button_uvalue = ROIProp_buttons, /no_release,set_value = 1)

velbase =  widget_base(mainbase, column = 1)
vel_button = ['X-Velocity','Y-Velocity','Z-Velocity','Speed','Composite']
vel_button_grp = cw_bgroup(velbase,vel_button, uvalue = 'vel_button', column = 1, exclusive = 1, button_uvalue = vel_button, label_top = 'Flow Output', frame = 2,/no_release , set_value=0)
visbase = widget_base(mainbase,row = 1)
vis_buttons = ['Image','Contour','2D Vector','2D Streamlines', '3D Surface', 'Orient', 'MIP', 'Cine-MIP']
vis_button_grp = cw_bgroup(visbase,vis_buttons,column = 1, exclusive = 1, uvalue = 'vis_button',button_uvalue = vis_buttons, label_top ='Flow Visualization', frame = 2,/no_release, set_value=0)

dispbase = widget_base(mainbase,row = 5)
uwbtnbase = widget_base(dispbase, row = 1, /nonexclusive)
pc_button = widget_button(uwbtnbase, value = 'Unwrap phase', uvalue = 'UwrapPhase')
ioverlay_button = widget_button(uwbtnbase, value = 'Overlay image', uvalue = 'ImageOverlay')

ECCbase = widget_base(dispbase, row = 7,frame =1)
ECCbtnbase = widget_base(ECCbase, /nonexclusive)
ecc_button = widget_button(ECCbtnbase, value = 'Phase Error Correction' , uvalue = 'ECCbtn')
if project.roi.ni eq 1 then begin
  ROIsetlist_val = 'None' 
  ROIlist_val = ' '
endif else begin
  ROIsetlist_val = (*project.roi.pDISPLAY_NAMES)[1:project.roi.ni-1]
  ROIlist_val = make_array(n_elements(*project.roi.pROIs[1]),/string)
  for i = 0, n_elements(*project.roi.pROIs[1])-1 do begin
    (*project.roi.pROIs[1])[i]->Getproperty,name = temp
    ROIlist_val[i] = temp
  endfor
endelse  
ROIsetlist1_ecc = widget_droplist(ECCbase,value = ROIsetlist_val ,uvalue = 'ECCROIsetlist1', title =  'Select a static ROI set 1', sensitive = 0)
ROIlist1_ecc = widget_droplist(ECCbase,value = ROIlist_val, uvalue = 'ECCROIlist1', title =  'Static ROI 1', sensitive = 0)
ROIsetlist2_ecc = widget_droplist(ECCbase,value = ROIsetlist_val ,uvalue = 'ECCROIsetlist2', title =  'Select a static ROI set 2', sensitive = 0)
ROIlist2_ecc = widget_droplist(ECCbase,value = ROIlist_val, uvalue = 'ECCROIlist2', title =  'Static ROI 2', sensitive = 0)
Fit_fns = ['C','X','Y','XX','YY','XY']
Fitfn_button_grp = cw_bgroup(ECCbase,Fit_fns, uvalue = 'ECCFitfns', row = 2, /nonexclusive, button_uvalue = Fitfns, label_top = 'Select correction function' , set_value=0)
vel_list = ['Read','Phase','Slice']
ecvel_button_grp = cw_bgroup(ECCbase,vel_list, uvalue = 'ECCvels', row = 1, /nonexclusive, button_uvalue = vel_list, label_top = 'Select velocities to apply' , set_value=0)

Ithres = widget_slider(dispbase, title = 'Image Threshold (%)', uvalue = 'Image Threshold',xsize = 180)
Vthres = widget_slider(dispbase, title = 'Velocity Threshold (%)', uvalue = 'Velocity Threshold',xsize = 180)

disp_button_base = widget_base(dispbase, row = 1)
disp_button = widget_button(disp_button_base, value = 'Display', uvalue = 'Display',/no_release, xsize = 82)
movie_button = widget_button(disp_button_base, value = 'Movie', uvalue = 'Movie',/no_release, xsize = 82)
exp_button = widget_button(disp_button_base, value = 'Export .raw', uvalue = 'Export',/no_release, xsize = 82)
 
if project.procPramArray[project.ci].single_Multi_flag eq 1 then $
widget_control, movie_button, sensitive = 1 else widget_control, movie_button, sensitive = 0 
widget_control, Fitfn_button_grp, sensitive = 0
widget_control, ecvel_button_grp, sensitive = 0
 
widget_control, topbase, /realize
xmanager, 'flow_gui', topbase, cleanup='flow_gui_cleanup',/NO_BLOCK

end

; Subroutine name: flow_gui_event
; Created by: Magdoom Kulam
; Calling Information:
;
;  event 
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Event handler for the main flow gui
;
; Editing Information:

pro flow_gui_event, event

common scan_data
common common_widgets

; To catch error ahead and prevent crashing
catch,error_status
IF (error_status NE 0) THEN BEGIN
  help, calls= trace_back
   dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
    !ERROR_STATE.msg,trace_back], /ERROR, $
    TITLE='Error in WID_BASE_MAIN_event')
    widget_control, event.top, /destroy
  RETURN
ENDIF

widget_control, event.id, get_uvalue = widget
CI = project.ci
 
if (event.top eq event.id) then begin
      if project.procPramArray[CI].single_Multi_flag eq 1 then widget_control, movie_button, sensitive = 1 else widget_control, movie_button, sensitive = 0 
      widget_control, event.top, tlb_set_title = 'Flow Imaging : Scan #' + strtrim(CI+1,2)             
      if (strcmp(project.imndarray[CI].flow,'y') and project.imndarray[CI].adim ge 4) then begin        
        get_flow_params, gflow, tflow, ncycl, shape, fc, scheme, venc        
        if ptr_valid(project.imndArray[project.ci].big_delta) then tbig_delta = *project.imndArray[project.ci].big_delta else tbig_delta = tflow/2.0
        widget_control,dispbase, sensitive = 1
        widget_control,parbase, sensitive = 1
        widget_control,gfield, set_value = strtrim(string(gflow,format = '(F8.2)'),2)+' mT/m'
        widget_control,tfield1, set_value = strtrim(string(tflow/2.0,format = '(F8.2)'),2)+' ms'
        widget_control,tfield2, set_value = strtrim(string(tbig_delta,format = '(F8.2)'),2)+' ms'
        widget_control,shpfield, set_value = shape
        widget_control,schfield, set_value = scheme
        widget_control,fcfield, set_value = fc
        widget_control,ncyclfield, set_value = ncycl
        widget_control,vencfield, set_value = strtrim(string(venc,format = '(F10.2)'),2) + ' um/s'
        widget_control,pc_button, set_button = project.procPramArray[CI].flow_up_stat
        widget_control,ecc_button, set_button = project.procPramArray[CI].flow_ecc_stat
        widget_control,Ithres, set_value = project.procPramArray[CI].flow_Ithreshold
        widget_control,Vthres, set_value = project.procPramArray[CI].flow_Vthreshold
        
        if project.procPramArray[CI].flow_ecc_stat eq 1 then begin
          widget_control, ROIsetlist1_ecc, sensitive = 1
          widget_control, ROIlist1_ecc, sensitive = 1
          widget_control, ROIsetlist2_ecc, sensitive = 1
          widget_control, ROIlist2_ecc, sensitive = 1
          widget_control, Fitfn_button_grp, sensitive = 1
          widget_control, ecvel_button_grp, sensitive = 1
        endif else begin
          widget_control, ROIsetlist1_ecc, sensitive = 0
          widget_control, ROIlist1_ecc, sensitive = 0
          widget_control, ROIsetlist2_ecc, sensitive = 0
          widget_control, ROIlist2_ecc, sensitive = 0
          widget_control, Fitfn_button_grp, sensitive = 0
          widget_control, ecvel_button_grp, sensitive = 0
        endelse
        ROIset1_index = project.procPramArray[CI].flow_ecc_roi1_set_index
        ROIsetlist_val = (*project.roi.pDISPLAY_NAMES)[0:project.roi.ni-1]
        if ptr_valid(project.roi.pROIs[ROIset1_index]) then begin
          ROI1list_val = make_array(n_elements(*project.roi.pROIs[ROIset1_index]),/string)
          for i = 0, n_elements(*project.roi.pROIs[ROIset1_index])-1 do begin
            (*project.roi.pROIs[ROIset1_index])[i]->Getproperty,name = temp
            ROI1list_val[i] = temp
          endfor  
        endif else ROI1list_val = ' '
                
        ROIset2_index = project.procPramArray[CI].flow_ecc_roi2_set_index
        if ptr_valid(project.roi.pROIs[ROIset2_index]) then begin
          ROI2list_val = make_array(n_elements(*project.roi.pROIs[ROIset2_index]),/string)
          for i = 0, n_elements(*project.roi.pROIs[ROIset2_index])-1 do begin
            (*project.roi.pROIs[ROIset2_index])[i]->Getproperty,name = temp
            ROI2list_val[i] = temp
          endfor  
        endif else ROI2list_val = ' '
                 
        widget_control, ROIsetlist1_ecc, set_value = ROIsetlist_val, set_droplist_select = project.procPramArray[CI].flow_ecc_roi1_set_index
        widget_control, ROIlist1_ecc, set_value = ROI1list_val, set_droplist_select = project.procPramArray[CI].flow_ecc_roi1_index
        widget_control, ROIsetlist2_ecc, set_value = ROIsetlist_val, set_droplist_select = project.procPramArray[CI].flow_ecc_roi2_set_index
        widget_control, ROIlist2_ecc, set_value = ROI2list_val, set_droplist_select = project.procPramArray[CI].flow_ecc_roi2_index
        widget_control, Fitfn_button_grp, set_value = project.procPramArray[CI].flow_ecc_fit_fn
        widget_control, ecvel_button_grp, set_value = project.procPramArray[CI].flow_ecc_vel
      
      endif else begin
          widget_control,dispbase, sensitive = 0
          widget_control,parbase, sensitive = 0
      endelse      
     
 endif else begin
      
    case widget of                                                              
    'vis_button'        : begin
                          widget_control, event.id, get_value = vis_type 
                          if vis_type ge 2 and vis_type le 5 and vis_type ne 4 then begin
                            widget_control, vel_button_grp, set_value = 0
                            widget_control, vel_button_grp, sensitive = 0                                                        
                          endif else begin
                            widget_control, vel_button_grp, sensitive = 1                            
                          endelse                                               
                          end                         
                     
   'UwrapPhase'         : begin
                          project.procPramArray[CI].flow_up_stat = event.select
                          project.procPramArray[CI].flow_proccess_flag = 0
                          end
   
   'ECCbtn'             : begin
                          project.procPramArray[CI].flow_ecc_stat = event.select
                          project.procPramArray[CI].flow_proccess_flag = 0
                          if event.select eq 1 then begin
                            widget_control, ROIsetlist1_ecc, sensitive = 1
                            widget_control, ROIlist1_ecc, sensitive = 1
                            widget_control, ROIsetlist2_ecc, sensitive = 1
                            widget_control, ROIlist2_ecc, sensitive = 1
                            widget_control, Fitfn_button_grp, sensitive = 1
                            widget_control, ecvel_button_grp, sensitive = 1
                          endif else begin
                            widget_control, ROIsetlist1_ecc, sensitive = 0
                            widget_control, ROIlist1_ecc, sensitive = 0
                            widget_control, ROIsetlist2_ecc, sensitive = 0
                            widget_control, ROIlist2_ecc, sensitive = 0
                            widget_control, Fitfn_button_grp, sensitive = 0
                            widget_control, ecvel_button_grp, sensitive = 0
                          endelse                          
                          end
   
   'ECCROIsetlist1'      : begin
                          widget_control, ROIsetlist1_ecc, get_value = temp
                          if temp[event.index] eq 'New ROI' then project.procPramArray[CI].flow_ecc_roi1_set_index = -1 $
                           else project.procPramArray[CI].flow_ecc_roi1_set_index = event.index                           
                          project.procPramArray[CI].flow_proccess_flag = 0
                          end
                                                    
   'ECCROIlist1'         : begin
                           project.procPramArray[CI].flow_ecc_roi1_index = event.index
                           project.procPramArray[CI].flow_proccess_flag = 0
                           end
   
   'ECCROIsetlist2'      : begin
                           widget_control, ROIsetlist2_ecc, get_value = temp
                           if temp[event.index] eq 'New ROI' then project.procPramArray[CI].flow_ecc_roi2_set_index = -1 $
                           else project.procPramArray[CI].flow_ecc_roi2_set_index = event.index                           
                           project.procPramArray[CI].flow_proccess_flag = 0 
                           end

    'ECCROIlist2'        :  begin
                            project.procPramArray[CI].flow_ecc_roi2_index = event.index
                            project.procPramArray[CI].flow_proccess_flag = 0
                            end                          
                          
   'ECCFitfns'           : begin
                           widget_control, Fitfn_button_grp , get_value = select
                           project.procPramArray[CI].flow_ecc_fit_fn = select                           
                           project.procPramArray[CI].flow_proccess_flag = 0
                           end                                                    

   'ECCvels'             : begin
                           widget_control, ecvel_button_grp , get_value = select
                           project.procPramArray[CI].flow_ecc_vel = select
                           project.procPramArray[CI].flow_proccess_flag = 0
                           end
                           
   'Image Threshold'    : project.procPramArray[CI].flow_Ithreshold = event.value
                          
   'Velocity Threshold' : project.procPramArray[CI].flow_Vthreshold = event.value 
                        
   'ROIPrint'           : ROI_analyze               
                          
   'ROIHist'            : ROI_hist
                         
   'Display'            : begin     
                          ioverlay = widget_info(ioverlay_button, /button_set) 
                          if (project.procPramArray[CI].flow_proccess_flag eq 0 or ~ptr_valid(project.dataarray[CI].flow) or project.procPramArray[CI].state_1 eq 0) then analyze_flow
                          if project.procPramArray[CI].single_Multi_flag eq 0 then flow_single_display, ioverlay else flow_multi_display, ioverlay                              
                          end                    
                                                                                 
   'Movie'              : begin                                          
                          if (project.procPramArray[CI].flow_proccess_flag eq 0 or ~ptr_valid(project.dataarray[CI].flow) or project.procPramArray[CI].state_1 eq 0) then analyze_flow
                          flow_movie                           
                          end  
   'Export'             : export_velocity                       
    else                : return                                                    
    endcase
endelse 
end

; Subroutine name: orient_flow_data
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To orient the velocity according right hand coordinate system
;
; Editing Information:
pro orient_flow_data, V

  M = log2physgradtransform()
  readax = max(M##[1,0,0],/absolute)
  phaseax = max(M##[0,1,0],/absolute)
  sliceax = max(M##[0,0,1],/absolute)
    
  if signum(readax) eq -1 then (*V)[*,*,*,0]*= -1   
  if signum(phaseax) eq -1 then (*V)[*,*,*,1]*= -1
  if signum(sliceax) eq -1 then (*V)[*,*,*,2]*= -1
  
end

; Subroutine name: flip_flow_data
; Created by: Magdoom Kulam
; Calling Information:
;
;  V - Pointer containing components of the 3D velocity and position vector along with the MR image in its fourth dimension.
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To flip the flow data as specified in the main MAS GUI
;
; Editing Information:

pro flip_flow_data,V

common scan_data

flip_Direction = project.procPramArray[project.ci].flip_Direction
if flip_Direction eq 0 then return
N = size(*V,/Dimension)
case flip_direction of
  1 : axis = 'y'
  2:  axis = 'x'
endcase

rotate_3d_vector_field, V, !pi, axis = axis        ; Flipping the velocity vector

; Flipping the image
for i = 0,N[3]-1 do begin
  for j = 0,N[2]-1 do begin
    (*V)[*,*,j,i] = reverse((*V)[*,*,j,i],flip_Direction)
  endfor
endfor

end

; Subroutine name: rotate_flow_data
; Created by: Magdoom Kulam
; Calling Information:
;
;  V - Pointer containing components of the 3D velocity and position vector along with the MR image in its fourth dimension. 
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To rotate the 3D velocity and position vector to the angle specified in the main MAS GUI
;
; Editing Information:

pro rotate_flow_data, V

common scan_data
 
rotate_Direction = project.procpramArray[project.ci].rotate_Direction  ; Rotation direction in MAS main widget
if rotate_Direction eq 0 then return
N = size(*V,/Dimension)
case rotate_direction of
1 : begin
    theta = 90*!pi/180 
    N[0:1] = reverse(N[0:1])
    end
2:  theta = 180*!pi/180
3 : begin
    theta = 270*!pi/180
    N[0:1] = reverse(N[0:1])
    end
endcase 

rotate_3d_vector_field, V, theta                              ; Rotating the velocity vector

temp = *V
V = ptr_new(make_array(N[0],N[1],N[2],5,value = 0,/double))   ; Resizing the velocity pointer based on rotation

; Rotating the image
for i = 0,N[3]-1 do begin
  for j = 0,N[2]-1 do begin
      (*V)[*,*,j,i] = rotate(temp[*,*,j,i],rotate_Direction)      
  endfor
endfor

end

; Subroutine name: unwrap_phase_flow
; Created by: Magdoom Kulam
; Calling Information:
; 
; V       - Pointer containing the flow data
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To unwrap the 3D velocity using FSL's PRELUDE method (default, if installed) or Fourier based method
;
; Editing Information:

pro unwrap_phase_flow, V

common scan_data

progressbar = Obj_New('progressbar', Color='red', Text='Unwrapping Phase', /NOCANCEL)
progressbar -> Start
Nz = (size(*V,/Dimension))[2]

for i = 0, Nz - 1 do begin
    for j = 0,2 do begin
      phase_image = (*V)[*,*,i,j]*!DPI
      phase_unwrap_2D,phase_image
      (*V)[*,*,i,j] = phase_image/!DPI
      progressBar -> Update, (float(i+1)*float(j+1)/(float(3)*float(Nz)))*100.0
    endfor
endfor  

(*V)[*,*,*,3] = sqrt(total((*V)[*,*,*,0:2]^2,4))
progressbar -> Destroy

end

; Subroutine name: mask_flow_data
; Created by: Magdoom Kulam
; Calling Information:
; 
; V  - Pointer containing the flow data
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To apply image and velocity threshold mask
;
; Editing Information:

pro mask_flow_data, V

common scan_data
common common_widgets

CI = project.ci
Ithreshold = project.procPramArray[CI].flow_Ithreshold
Vthreshold = project.procPramArray[CI].flow_Vthreshold

N1 = size(*project.dataarray[CI].state1,/Dimension) ; State 1 dimensions
N2 = size(*V,/Dimension)                            ; State 2 dimensions

for i = 0, N2[2]-1 do begin
    Imask = abs((*V)[*,*,i,4]) ge max(abs((*V)[*,*,i,4]))*float(Ithreshold)/100
    if isa(ROImask) then Imask *= ROImask
    for j = 0,3 do (*V)[*,*,i,j] *= Imask*(abs((*V)[*,*,i,j]) ge max(abs((*V)[*,*,i,j]))*float(Vthreshold)/100)    
    (*V)[*,*,i,4]*=Imask
endfor

end

; Subroutine name: analyze_flow
; Created by: Magdoom Kulam
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To calculate the velocity vector from PC-MRI acquistion via complex division based 
;                        on the acquistion scheme, Hadamard by default (other : 4-point, 6-point)
;
; Editing Information:
; 1) Fixed sign error in the velocity calculation (Magdoom, 1/30/17)

pro analyze_flow

common scan_data

CI = project.ci

; Loading complex data
set_signal_type, 9
mas_load_state_1

data = *project.dataarray[CI].state1
N = size(data,/Dimension)
if ptr_valid(project.dataarray[CI].flow) then ptr_free,project.dataarray[CI].flow
project.dataarray[CI].flow = ptr_new(make_array(N[0],N[1],N[2],5,value = 0,/double))

; Computing the velocity using complex division based on the encoding scheme
scheme = project.imndarray[CI].scheme
case scheme of
'f'  : begin
       (*project.dataarray[CI].flow)[*,*,*,0] = atan(data[*,*,*,3]/data[*,*,*,0],/Phase)/!pi
       (*project.dataarray[CI].flow)[*,*,*,1] = atan(data[*,*,*,3]/data[*,*,*,1],/Phase)/!pi
       (*project.dataarray[CI].flow)[*,*,*,2] = atan(data[*,*,*,3]/data[*,*,*,2],/Phase)/!pi
       end
's'  : begin
       (*project.dataarray[CI].flow)[*,*,*,0] = atan(data[*,*,*,1]/data[*,*,*,0],/Phase)/!pi
       (*project.dataarray[CI].flow)[*,*,*,1] = atan(data[*,*,*,3]/data[*,*,*,2],/Phase)/!pi
       (*project.dataarray[CI].flow)[*,*,*,2] = atan(data[*,*,*,5]/data[*,*,*,4],/Phase)/!pi
       end  
else : begin
       (*project.dataarray[CI].flow)[*,*,*,0] = atan((data[*,*,*,0]*data[*,*,*,3])/(data[*,*,*,2]*data[*,*,*,1]),/Phase)/!pi
       (*project.dataarray[CI].flow)[*,*,*,1] = atan((data[*,*,*,0]*data[*,*,*,2])/(data[*,*,*,3]*data[*,*,*,1]),/Phase)/!pi
       (*project.dataarray[CI].flow)[*,*,*,2] = atan((data[*,*,*,0]*data[*,*,*,1])/(data[*,*,*,2]*data[*,*,*,3]),/Phase)/!pi
       end
endcase
(*project.dataarray[CI].flow)[*,*,*,3] = sqrt(total((*project.dataarray[CI].flow)[*,*,*,0:2]^2,4))
(*project.dataarray[CI].flow)[*,*,*,4] = mean(abs(data),dimension = 4)

if project.procPramArray[CI].flow_up_stat eq 1 then unwrap_phase_flow,project.dataarray[CI].flow
if project.procPramArray[CI].flow_ecc_stat eq 1 then ecc,project.dataarray[CI].flow
if ptr_valid(project.dataarray[CI].flow) then project.procPramArray[CI].flow_proccess_flag = 1

end

; Subroutine name: ROIextract_vel
; Created by: Magdoom Kulam
; Calling Information:
; 
; V     - State 1 processed flow data 
; Vroi  - 4 by N pointer array containing the components of velocity vector & speed in each rois in the active ROI set
; index - Optional keyword to select a particular ROI in the active ROI set 
; 
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To extract the velocity in the active ROI set
; 
; Editing Information

pro ROIextract_vel, Vroi, index = index

common scan_data

CI = project.ci
N1 = size(*project.dataarray[CI].state1,/Dimension) ; State 1 dimensions
rotate_Direction = project.procpramArray[CI].rotate_Direction
islice = project.procpramArray[CI].sdim_start
V = project.dataarray[CI].flow
N2 = size(*V,/Dimension)                            ; State 2 dimensions
get_flow_params, gflow, tflow, ncycl, shape, fc, scheme, venc

if ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
  Vroi = ptrarr(4,n_elements(*project.roi.pROIs[project.roi.ci]),/allocate_heap)
  for i = 0, n_elements(*project.roi.pROIs[project.roi.ci])-1 do begin 
      if isa(index) then if i ne index then continue
      ROImask = congrid(rotate(((*project.roi.pROIs[project.roi.ci])[i] -> ComputeMask(Dimensions= [N1[0],N1[1]]))/255, rotate_direction),N2[0],N2[1],interp = 0) 
      for j = 0,2 do *Vroi[j,i] = venc*((*V)[*,*,islice,j])[where(ROImask ne 0)]
      *Vroi[3,i] = sqrt(*Vroi[0,i]^2+*Vroi[1,i]^2+*Vroi[2,i]^2)
  endfor 
endif else Vroi = !NULL

end

; Subroutine name: ROI_hist
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To plot histogram of velocity (only non-zero values) in the active ROI
; 
; Editing Information

pro ROI_hist

common scan_data
common common_widgets

widget_control, vel_button_grp, get_value = vel_type
get_flow_params, gflow, tflow, ncycl, shape, fc, scheme, venc

CI = project.ci
if (project.procPramArray[CI].flow_proccess_flag eq 0 or ~ptr_valid(project.dataarray[CI].flow)) then analyze_flow
ROIextract_vel, Vroi, index = 0
if ~isa(Vroi) then begin
  mbox = dialog_message('Please make an ROI to display its histogram',/error)
  return
endif

if vel_type eq 4 then begin
  for i=0,3 do begin
    Vselect = *Vroi[i,0]
    if n_elements(where(Vselect ne 0)) gt 1 then hist = histogram(Vselect[where(Vselect ne 0)], LOCATIONS = xbin, BINSIZE = 0.05*venc) else begin $
      mbox = dialog_message([strtrim(vel_button(i),2) + '- Velocity is zero in the ROI'])
      continue
    endelse
    if i eq 3 then range = [0,sqrt(3)*venc] else range = [-venc,venc]
    if i eq 0 then current_kw = 0 else current_kw = 1
    hplot = BARPLOT(xbin, hist, XRANGE = range, TITLE='ROI Histogram', CURRENT = current_kw, layout = [3,2,i+1], $
                XTITLE = vel_button(i) + ' '+'Velocity (um/s)', YTITLE = 'Frequency', AXIS_STYLE=1, COLOR='red')           
  endfor
endif else begin
  Vselect = *Vroi[vel_type,0]
  if n_elements(where(Vselect ne 0)) gt 1 then hist = histogram(Vselect[where(Vselect ne 0)], LOCATIONS = xbin, BINSIZE = 0.05*venc) else begin $
      mbox = dialog_message([strtrim(vel_button(vel_type),2) + '- Velocity is zero in the ROI'])
      if (ptr_valid(Vroi))[0] then ptr_free,Vroi
      return
  endelse
  if vel_type eq 3 then range = [0,2*venc] else range = [-venc,2*venc]
  IMAGE_STATISTICS, Vselect, MAXIMUM=ma, MINIMUM=mi, MEAN=me, STDDEV=s
  hplot = BARPLOT(xbin, hist, XRANGE = range, TITLE = 'ROI Histogram', $
                XTITLE = vel_button(vel_type) + ' '+'Velocity (um/s)', YTITLE = 'Frequency', AXIS_STYLE = 1, COLOR = 'red')
  t1 = TEXT(venc, max(hist), 'Max   :  '+ strtrim(string(1e3*ma,format = '(F7.2)'),2) +' $\mu$m/s', /DATA)
  t2 = TEXT(venc, 9*max(hist)/(10.0), 'Min    :  '+strtrim(string(1e3*mi,format = '(F7.2)'),2) +' $\mu$m/s', /DATA)
  t3 = TEXT(venc, 8*max(hist)/(10.0), 'Mean :  '+strtrim(string(1e3*me,format = '(F7.2)'),2) +' $\mu$m/s',/DATA)
  t4 = TEXT(venc, 7*max(hist)/(10.0), 'SD      :  '+strtrim(string(1e3*s,format = '(F7.2)'),2) +' $\mu$m/s',/DATA)    
endelse
if (ptr_valid(Vroi))[0] then ptr_free,Vroi

end

; Subroutine name: ROI_analyze
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To calculate the volumetric flow rate (mL/min) and flow anisotropy index for the active ROI
;
; Editing Information

pro ROI_analyze

common scan_data
common common_widgets

CI = project.ci
if (project.procPramArray[CI].flow_proccess_flag eq 0 or ~ptr_valid(project.dataarray[CI].flow)) then analyze_flow
get_flow_params, gflow, tflow, ncycl, shape, fc, scheme, venc
ROIextract_vel, Vroi

if ~isa(Vroi) then begin
  mbox = dialog_message('Please make an ROI to analyze',/error)
  return
endif

topbase = widget_base(title = 'Flow ROI Statistics' ,tlb_frame_attr=1,/grid_layout, frame = 1)
mainbase = widget_base(topbase)
widget_control, ROI_button_grp, get_value = ROIproperty

rotate_Direction = project.procpramArray[project.ci].rotate_Direction
if (rotate_Direction eq 0 or rotate_Direction eq 2) then begin
    dX = 10*project.imndarray[CI].f_voxsz     ; X resolution (mm)
    dY = 10*project.imndarray[CI].p_voxsz     ; Y resolution (mm)
endif else begin
    dX = 10*project.imndarray[CI].p_voxsz     ; X resolution (mm)
    dY = 10*project.imndarray[CI].f_voxsz     ; Y resolution (mm)
endelse    
dZ = project.imndarray[CI].s_voxsz         ; Z resolution (mm)

rowtxt1 = '  '  
rowtxt2 = 'ROIName'
for i = 0, n_elements(*project.roi.pROIs[project.roi.ci])-1 do begin
  (*project.roi.pROIs[project.roi.ci])[i]->Getproperty,name = ROIname
  temp = strtrim(ROIname,2)+' ' 
  if ROIproperty[0] eq 1 then begin
      if i eq 0 then begin
        rowtxt1 += '  FlowRate(uL/min) '
        rowtxt2 += ' YZ plane  XZ plane  XY plane ' 
      endif
      FR = 0.06*[total(*Vroi[0,i])*dY*dZ, total(*Vroi[1,i])*dX*dZ, total(*Vroi[2,i])*dX*dY]
      temp += strjoin(string(FR,format='(f10.4)'))+ '  '
  endif
  if ROIproperty[1] eq 1 then begin
    T = [mean(*Vroi[0,i]^2),mean(*Vroi[1,i]^2),mean(*Vroi[2,i]^2)]/(mean(*Vroi[0,i]^2)+mean(*Vroi[1,i]^2)+mean(*Vroi[2,i]^2))
    if i eq 0 then begin
      rowtxt1 += '  FlowAnisotropy  '
      rowtxt2 += '  X  Y  Z '
    endif  
    temp += strjoin(string(T,format='(f10.2)'))+'  '
  endif
  if ROIproperty[2] eq 1 then begin
    if i eq 0 then begin
      rowtxt1 += '  Flow velocity(um/s) '
      rowtxt2 += ' X  Y  Z '
    endif
    Velocity = [mean(*Vroi[0,i]), mean(*Vroi[1,i]), mean(*Vroi[2,i])]
    temp += strjoin(string(Velocity,format='(f10.4)'))+ '  '
  endif
  if i eq 0 then rowtxt3 = temp else rowtxt3 = [[rowtxt3],[temp]]
endfor 
txt = [[rowtxt1],[rowtxt2],[rowtxt3]]
rowlabel = widget_text(mainbase, value = txt, ysize = n_elements(*project.roi.pROIs[project.roi.ci])+2)
widget_control, topbase, /realize

if (ptr_valid(Vroi))[0] then ptr_free,Vroi

end

; Subroutine name: ECC_function
; Created by: Magdoom Kulam
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Function to be fit for eddy current correction
;
; Editing Information:

pro ECC_Function, XY,A,V,PDER
common scan_data

X = XY[*,0]
Y = XY[*,1]
Fn = [[replicate(1.0, N_ELEMENTS(X))],[X],[Y],[X*X],[Y*Y],[X*Y]]

fit = project.procPramArray[project.ci].flow_ecc_fit_fn
index = where(fit eq 1)
V = 0
PDER = []
for i = 0,n_elements(index)-1 do begin
  V+=A[i]*Fn[*,index(i)]
  PDER = [[PDER],[Fn[*,index(i)]]]
endfor 

end

; Subroutine name: ECC
; Created by: Magdoom Kulam
; Calling Information:
; 
; V  - Pointer containing the flow data
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To correct the eddy current artifacts using a static ROI (first region in the selected ROI set in flow widget)
;
; Editing Information:

pro ECC,V
common scan_data
common common_widgets

CI = project.ci
N = size(*V,/Dimension)
islice = project.procpramArray[CI].sdim_start
eccmask1 = make_Array(N[0],N[1],value = 0)
eccmask2 = make_Array(N[0],N[1],value = 0)
ROIset1_index = project.procPramArray[CI].flow_ecc_roi1_set_index
ROI1_index = project.procPramArray[CI].flow_ecc_roi1_index

ROIset2_index = project.procPramArray[CI].flow_ecc_roi2_set_index
ROI2_index = project.procPramArray[CI].flow_ecc_roi2_index

if ptr_valid(project.roi.pROIs[ROIset1_index]) then eccmask1 = (*project.roi.pROIs[ROIset1_index])[ROI1_index] -> ComputeMask(Dimensions= [N[0],N[1]]) 
if ptr_valid(project.roi.pROIs[ROIset2_index]) then eccmask2 = (*project.roi.pROIs[ROIset2_index])[ROI2_index] -> ComputeMask(Dimensions= [N[0],N[1]])
if ~ptr_valid(project.roi.pROIs[ROIset1_index]) and ~ptr_valid(project.roi.pROIs[ROIset2_index]) then begin
  mbox = dialog_message(['Please make an static ROI', 'for eddy current correction'],/error,/center)
  if ptr_valid(V) then ptr_free,V
  return
end

Fn_ind = where(project.procPramArray[CI].flow_ecc_fit_fn eq 1,/NULL)
if n_elements(Fn_ind) eq 0 then begin
  mbox = dialog_message(['Please select atleast one term in the', 'eddy current correction function to fit'],/error,/center)
 if ptr_valid(V) then ptr_free,V  
 return
end 
V_ind = where(project.procPramArray[CI].flow_ecc_vel eq 1, /NULL)

x = findgen(N[0])/float(N[0]-1)
y = findgen(N[1])/float(N[1]-1)
meshgrid,x,y

indices1 = where(eccmask1 gt 0,/NULL)
indices2 = where(eccmask2 gt 0,/NULL) 

if isa(indices1) then XY = [[x[indices1]],[y[indices1]]]
if isa(indices2) then XY = [XY,[[x[indices2]],[y[indices2]]]]
Fn = make_array(N[0],N[1],6)
Fn[*,*,0:5] = [[make_array(N[0],N[1],value = 1, /double)],[X],[Y],[X*X],[Y*Y],[X*Y]]
fit_result = make_array(n_elements(Fn_ind)+1,3,N[2],value=0,/double)

progressbar = Obj_New('progressbar', Color='red', /NOCANCEL)
progressbar -> Start

for i = 0,N[2]-1 do begin 
    if project.procPramArray[CI].single_Multi_flag eq 0 and i ne islice then continue                                              
    for j = 0,n_elements(V_ind)-1 do begin  
        progressbar->SetProperty, Text='Performing ECC for slice#' +strtrim(i,2)         
        temp = (*V)[*,*,i,V_ind[j]]
        if isa(indices1) then Z = temp[indices1]
        if isa(indices2) then Z = [Z,temp[indices2]]        
        A = replicate(1.0,n_elements(Fn_ind))
        yfit = CURVEFIT(XY, Z, weights, A, FUNCTION_NAME='ECC_Function', CHISQ = chi,/double, status = status)
        SS_tot  = variance(Z,/double)
        rsquared     = 1 - chi/SS_tot 
        if rsquared gt 0 then for k = 0,n_elements(Fn_ind)-1 do (*V)[*,*,i,V_ind[j]] -= A[k]*Fn[*,*,Fn_ind[k]]               
        fit_result[*,V_ind[j],i] = [A, rsquared]
        if project.procPramArray[CI].single_Multi_flag eq 0 then progressBar -> Update, (float(j+1)/(float(3)))*100.0 $
        else progressBar -> Update, (float(i+1)*float(j+1)/(float(3)*float(N[2])))*100.0                         
    endfor         
endfor            
(*V)[*,*,*,3] = sqrt((*V)[*,*,*,0]^2+(*V)[*,*,*,1]^2+(*V)[*,*,*,2]^2)
progressbar -> Destroy
display_ecc_stats, fit_result

end

; Subroutine name: display_ecc_stats
; Created by: Magdoom Kulam
; Calling Information:
; 
; A - Matrix containing eddy current correction fit parameters
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display the eddy current fitted parameters
;
; Editing Information:

pro display_ecc_stats, A

common scan_data

CI = project.ci

FOVx = 10*project.imndarray[CI].f_fov               ; X-FOV (mm)
FOVy = 10*project.imndarray[CI].p_fov               ; Y-FOV (mm)  

ECC_vel = project.procPramArray[CI].flow_ecc_vel
V_ind = where(ECC_vel eq 1, /NULL)

if project.procPramArray[project.ci].single_Multi_flag eq 0 then islice = project.procpramArray[project.ci].sdim_start else $
  islice = indgen(project.imndArray[project.ci].sdim)

topbase = widget_base(title = 'Eddy Current Correction ROI Fit Statistics' ,tlb_frame_attr=1,/grid_layout, frame = 1)

mainbase = widget_base(topbase)
Fn_ind = where(project.procPramArray[project.ci].flow_ecc_fit_fn eq 1,/NULL)
F = ['Intercept (rad) ','   Read (rad/mm)  ','   Phase (rad/mm)  ',' Read-Read (rad/mm^2)',' Phase-Phase (rad/mm^2)',' Read-Phase (rad/mm^2)']
scaling_factor = !Pi*[1, 1/float(FOVx), 1/float(FOVy), 1/float(FOVx^2), 1/float(FOVy^2), 1/float(FOVx*FOVy)]
iscale = [scaling_factor[Fn_ind], 1.0]

rowtxt1 = '        '
rowtxt2 = '        ' 
if ECC_vel[0] eq 1 then begin
   rowtxt1 += '      Read - velocity      ' 
   rowtxt2 += strjoin(F[Fn_ind],' ') + ' Rsquared    '
endif 
if ECC_vel[1] eq 1 then begin
  rowtxt1 += '      Phase - velocity      ' 
  rowtxt2 += strjoin(F[Fn_ind],' ') + ' Rsquared    '
endif 
if ECC_vel[2] eq 1 then begin
  rowtxt1 += '      Slice - velocity      '
  rowtxt2 += strjoin(F[Fn_ind],' ') + ' Rsquared    '
endif
for i = 0, n_elements(islice)-1 do begin
  temp = 'Slice' + strtrim(islice[i],2) +' '
  for j = 0,n_elements(V_ind)-1 do begin
    for k = 0, n_elements(Fn_ind) do temp += strtrim(string(iscale[k]*A[k,V_ind[j],islice[i]],format='(f7.3)'),1) +'   '
  endfor
    if i eq 0 then rowtxt3 = temp else rowtxt3 = [rowtxt3,temp]
endfor

txt = [rowtxt1,rowtxt2,rowtxt3]
rowlabel = widget_text(mainbase, value = txt, ysize = n_elements(islice)+2 )

widget_control, topbase, /realize

end    

; Subroutine name: flow_orient_map
; Created by: Magdoom Kulam
; Calling Information:
; 
; V - Pointer containing flow data
; islice - Optional parameter to specify the slice (index) to display
; 
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display the RGB flow orientation map for single or multiple slices 
;
; Editing Information:

pro flow_orient_map, V, islice = islice, venc

N = size(*V,/Dimensions)

x = float(255)*indgen(N[0])/float(N[0]-1)
y = float(255)*indgen(N[1])/float(N[1]-1)
meshgrid,x,y

top = make_array(3,N[0],N[1])
LL = make_array(3,N[0],N[1])
LR = make_array(3,N[0],N[1])

top[0,*,*] = bytscl(x)
top[1,*,*] = make_array(N[0],N[1],value=127,/integer)
top[2,*,*] = bytscl(y)

LL[0,*,*] = make_array(N[0],N[1],value=127,/integer)
LL[1,*,*] = bytscl(y)
LL[2,*,*] = bytscl(x)

LR[0,*,*] = bytscl(x)
LR[1,*,*] = bytscl(y)
LR[2,*,*] = make_array(N[0],N[1],value=127,/integer)
  
; Change to the Z buffer device
SET_PLOT, 'Z'
DEVICE, SET_PIXEL_DEPTH=24, SET_RESOLUTION=[N[0],N[1]]

; Establish 3-D scaling as (0,1) cube:
SCALE3, XRANGE=[0,1], YRANGE=[0,1], ZRANGE=[0,1]
; Define vertices of cube. Vertices 0-3 are bottom, 4-7 are top:
verts = [[0,0,0], [1,0,0], [1,1,0], [0,1,0], [0,0,1], [1,0,1], [1,1,1], [0,1,1]]

; Fill lower left face:
POLYFILL, verts[*, [3,0,4,7]], /T3D, PATTERN=LL,IMAGE_COORD=[[0,0], [N[0]/2,0], [N[0]/2,N[1]/2], [0,N[1]/2]]

; Fill lower right face:
POLYFILL, verts[*, [0,1,5,4]], /T3D, PATTERN=LR, IMAGE_COORD=[[N[0]/2,0], [N[0]-1,0], [N[0]-1,N[1]/2], [N[0]/2,N[1]/2]]

; Fill top face:
POLYFILL, verts[*, [4,5,6,7]], /T3D, PATTERN=top, IMAGE_COORD = [[N[0]/2,N[1]/2], [N[0]-1,N[1]/2], [N[0]-1,N[1]-1], [N[0]/2,N[1]-1]]

; Draw edges of cube in black:
PLOTS, verts[*, [0,4]], /T3D, COLOR=0

; Edges of top face:
PLOTS, verts[*, [4,5,6,7,4]], /T3D, COLOR=0
cube = TVRD(/TRUE)

; Change back to the original device and display
if !version.os_family eq 'Windows' then SET_PLOT, 'WIN' else SET_PLOT,'X'

im = make_array(N[0],N[1],3)
case isa(islice) of
  0 : begin
      WINDOW, YSIZE = floor(sqrt(N[2]))*N[0], XSIZE = (1+floor(N[2]/floor(sqrt(N[2]))))*N[1]
      for i=0,N[2]-1 do begin
        for j=0,2 do im[*,*,j] = bytscl((*V)[*,*,i,j], min = -venc, max = +venc)   ; Assigning RGB channels
        TV, im,i,TRUE = 3 ; Display the image
      endfor
      TV, cube,N[2], TRUE=1
      end
      
  1 : begin
      WINDOW, XSIZE = 2*N[0], YSIZE = N[1]
      for i=0,2 do im[*,*,i] = bytscl((*V)[*,*,islice,i], min = -venc, max = +venc)   ; Assigning RGB channels
      TV, im,0,0,TRUE = 3 ; Display the image 
      TV, cube,N[0],0, TRUE=1
      end
endcase
close,/all

end

; Subroutine name: flow_single_display
; Created by: Magdoom Kulam
; Calling Information:
; 
;    
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display individual components of the velocity vector for a single slice
;
; Editing Information:

pro flow_single_display, ioverlay
common scan_data
common common_widgets

widget_control, vel_button_grp, get_value = vel_type
widget_control, vis_button_grp, get_value = vis_type

CI = project.ci
V = load_state2_data_flow(project.dataarray[CI].flow)
if ~ptr_valid(V) then return

islice = project.procpramArray[CI].sdim_start
rotate_Direction = project.procpramArray[project.ci].rotate_Direction
if (rotate_Direction eq 0 or rotate_Direction eq 2) then begin
    dX = 10*project.imndarray[CI].f_voxsz     ; X resolution (mm)
    dY = 10*project.imndarray[CI].p_voxsz     ; Y resolution (mm)
endif else begin
    dX = 10*project.imndarray[CI].p_voxsz     ; X resolution (mm)
    dY = 10*project.imndarray[CI].f_voxsz     ; Y resolution (mm)
endelse    

widget_control,vencfield, get_value = vencstring
venc = (strsplit(vencstring,/extract))[0]

N = size(*V,/Dimensions)
if (vel_type eq 4) then begin
  composite,V, islice, dX, dY, venc, ioverlay
  return
endif 
if vel_type eq 3 then range = [0, sqrt(3)*venc] else range = [-venc, +venc]
  
case vis_type of
0:  begin
    if ioverlay eq 1 then begin
      iimage,(*V)[*,*,islice,4], rgb_table = 0, aspect_ratio =  dY/dX
      iimage, (*V)[*,*,islice,vel_type], identifier = iid, aspect_ratio =  dY/dX, rgb_table=33, /insert_colorbar, title = vel_button(vel_type) + '(um/s)', min_value = range[0], max_value = range[1], /overplot, transparency = 50
    endif else iimage, (*V)[*,*,islice,vel_type], identifier = iid, aspect_ratio =  dY/dX, rgb_table=33, /insert_colorbar, title = vel_button(vel_type) + '(um/s)', min_value = range[0], max_value = range[1]
         
    end
1:  begin   
  if ioverlay eq 1 then begin
    h1 = image((*V)[*,*,islice,4], aspect_ratio =  dY/dX)
    h2 = contour((*V)[*,*,islice,vel_type], aspect_ratio =  dY/dX, n_levels =50, c_label_show = 0, RGB_TABLE=33, title = vel_button(vel_type) + '(um/s)',min_value = range[0], max_value = range[1], /overplot, transparency = 50)
    endif else h2 = contour((*V)[*,*,islice,vel_type], aspect_ratio =  dY/dX, n_levels =50, c_label_show = 0, RGB_TABLE=33, title = vel_button(vel_type) + '(um/s)',min_value = range[0], max_value = range[1], axis_style = 0)    
    end
2:  begin
    x = findgen(N[0])
    y = findgen(N[1])
    if ioverlay eq 1 then begin
      h1 = image((*V)[*,*,islice,4], aspect_ratio = dY/dX)
      h2 = vector((*V)[*,*,islice,0], (*V)[*,*,islice,1], x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1, arrow_style = 0, /current,/overplot)  
    endif else begin
      h2 = vector((*V)[*,*,islice,0], (*V)[*,*,islice,1], x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1, arrow_style = 0)
      ax =h2.axes
      ax[0].hide = 1
      ax[1].hide = 1
    endelse      
    end   
    
3:  begin
    x = FINDGEN(N[0])
    y = FINDGEN(N[1]) 
    if ioverlay eq 1 then begin
      h1 = image((*V)[*,*,islice,4], aspect_ratio = dY/dX)
      h2 = streamline(reform((*V)[*,*,islice,0],N[0],N[1]), reform((*V)[*,*,islice,1],N[0],N[1]), x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1,/current,/overplot,thick = 2)  
    endif else begin
      h2 = streamline(reform((*V)[*,*,islice,0],N[0],N[1]), reform((*V)[*,*,islice,1],N[0],N[1]), x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1, thick = 2)
      ax = h2.axes
      ax[0].hide = 1
      ax[1].hide = 1
    endelse
    end
   
4: ivolume,(*V)[*,*,*,vel_type],RENDER_QUALITY = 2, composite_function=1, rgb_table0 =33,/ZERO_OPACITY_SKIP,RENDER_EXTENTS = 1,/INTERPOLATE    
     
5:  flow_orient_map, V, islice = islice, venc
6:  mas_mip,(*V)[*,*,*,vel_type], 0
7:  mas_mip,(*V)[*,*,*,vel_type], 1
endcase

if ptr_valid(V) then ptr_free,V
close,/all

end

; Subroutine name: composite
; Created by: Magdoom Kulam
; Calling Information:
;  
;  V      - Pointer containing the flow data
;  islice - Selected slice index to display
;  dX   - Resolution in 'x'
;  dY   - Resolution in 'y' 
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display all the components of the velocity vector for a single slice
;
; Editing Information:

pro composite, V, islice, dX, dY, venc, ioverlay

common common_widgets

widget_control, vis_button_grp, get_value = vis_type
N = size(*V,/Dimensions)
 
case vis_type of

0:  begin
    iimage
    for i=1,4 do begin 
      if i eq 4 then range = [0, sqrt(3)*venc] else range = [-venc, +venc]
      if ioverlay eq 1 then begin
        iimage,(*V)[*,*,islice,4], layout = [3,2,i], aspect_ratio = dY/dX, margin = 0.15 , /current       
        iimage,(*V)[*,*,islice,i-1], RGB_TABLE=33, layout = [3,2,i], title = vel_button[i-1], font_size=15, aspect_ratio = dY/dX, min_value = range[0], max_value = range[1], /overplot, transparency = 50       
      endif else iimage,(*V)[*,*,islice,i-1], /current, RGB_TABLE=33, layout = [3,2,i], title = vel_button[i-1], font_size=15, aspect_ratio = dY/dX, min_value = range[0], max_value = range[1]
     endfor
     if ioverlay eq 0 then iimage,(*V)[*,*,islice,4], layout = [3,2,5], aspect_ratio = dY/dX, margin = 0.15, title = 'MR image', font_size = 15, /current
      a1 = ARROW([0.8,0.9], [0.2,0.2], /CURRENT)
      h1 = TEXT(0.91,0.19,'X')
      a2 = ARROW([0.8,0.8], [0.2,0.35], /CURRENT)
      h2 = TEXT(0.79,0.38,'Y') 
      a3 = ARROW([0.8,0.75], [0.2,0.1], /CURRENT)
      h3 = TEXT(0.73,0.05,'Z')    
     end
     
1:  begin
    h = image(make_array(5,5),/NODATA)
    for i=1,4 do begin      
      if i eq 4 then range = [0, sqrt(3)*venc] else range = [-venc, +venc]
      if ioverlay eq 1 then begin
        h1 = image((*V)[*,*,islice,4], /current, layout = [3,2,i], aspect_ratio = dY/dX)
        h2 = contour((*V)[*,*,islice,i-1], n_levels = 50, c_label_show = 0, RGB_TABLE=33, layout = [3,2,i], title = vel_button[i-1], font_size=15, aspect_ratio = dY/dX, min_value = range[0], max_value = range[1],transparency = 50, /overplot)
      endif else h2 = contour((*V)[*,*,islice,i-1], /current, axis_style = 0, n_levels = 50, c_label_show = 0, RGB_TABLE=33, layout = [3,2,i], title = vel_button[i-1], font_size=15, aspect_ratio = dY/dX, min_value = range[0], max_value = range[1])
     endfor
     if ioverlay eq 0 then iimage,(*V)[*,*,islice,4], layout = [3,2,5], aspect_ratio = dY/dX, margin = 0.15, title = 'MR image', font_size = 15, /current
       a1 = ARROW([0.8,0.9], [0.2,0.2], /CURRENT)
      h1 = TEXT(0.91,0.19,'X')
      a2 = ARROW([0.8,0.8], [0.2,0.35], /CURRENT)
      h2 = TEXT(0.79,0.38,'Y') 
      a3 = ARROW([0.8,0.75], [0.2,0.1], /CURRENT)
      h3 = TEXT(0.73,0.05,'Z')    
     end
else: return

endcase               
close,/all

end

; Subroutine name: flow_multi_display
; Created by: Magdoom Kulam
; Calling Information:
;  
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display individual components of the velocity vector for all slices
;
; Editing Information:

pro flow_multi_display, ioverlay

common scan_data
common common_widgets

widget_control, vel_button_grp, get_value = vel_type
widget_control, vis_button_grp, get_value = vis_type
if vel_type eq 4 then return

CI = project.ci
V = load_state2_data_flow(project.dataarray[CI].flow)
if ~ptr_valid(V) then return

rotate_Direction = project.procpramArray[CI].rotate_Direction
if rotate_Direction eq 0 or rotate_Direction eq 2 then begin
  dX = 10*project.imndarray[CI].f_voxsz     ; X resolution (cm)
  dY = 10*project.imndarray[CI].p_voxsz     ; Y resolution (cm)
endif else begin
  dX = 10*project.imndarray[CI].p_voxsz     ; X resolution (cm)
  dY = 10*project.imndarray[CI].f_voxsz     ; Y resolution (cm)
endelse

widget_control,vencfield, get_value = vencstring
venc = (strsplit(vencstring,/extract))[0]
if vel_type eq 3 then range = [0, sqrt(3)*venc] else range = [-venc, +venc]

N = size(*V,/Dimensions)
case vis_type of

0:  begin       
    iimage 
    for i = 1,N[2] do begin       
      if ioverlay eq 1 then begin
       iimage,(*V)[*,*,i-1,4], aspect_ratio = dY/dX, /current, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i], margin = 0.1
       iimage,(*V)[*,*,i-1,vel_type], aspect_ratio = dY/dX, RGB_TABLE=33, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i], min_value = range[0], max_value = range[1], /overplot, transparency = 50
     endif else   iimage,(*V)[*,*,i-1,vel_type], aspect_ratio = dY/dX, /current, RGB_TABLE=33, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i], min_value = range[0], max_value = range[1]
    endfor
    end
1:  begin
    h = image(make_array(5,5),/NODATA)
    for i = 1,N[2] do begin
      if ioverlay eq 1 then begin
        h1 = image(rebin((*V)[*,*,i-1,4],N[0],N[1]), /current, aspect_ratio = dY/dX, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i], margin = 0.1)
        h2 = contour((*V)[*,*,i-1,vel_type], /overplot, aspect_ratio = dY/dX, n_levels =50, c_label_show = 0, RGB_TABLE=33, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i], min_value = range[0], max_value = range[1], transparency = 50)
      endif else h2 = contour((*V)[*,*,i-1,vel_type], /current, axis_style = 0, aspect_ratio = dY/dX, n_levels =50, c_label_show = 0, RGB_TABLE=33, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i], min_value = range[0], max_value = range[1])     
    endfor
    end   
2:  begin
    for i = 0,N[2]-1 do begin
      x = findgen(N[0])
      y = findgen(N[1])
      if ioverlay eq 1 then begin
        h1 = image((*V)[*,*,i,4], aspect_ratio = dY/dX,layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i+1],/current)
        h2 = vector((*V)[*,*,i,0], (*V)[*,*,i,1], x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1, arrow_style = 0, /current,/overplot,layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i+1])        
      endif else begin
        h2 = vector((*V)[*,*,i,0], (*V)[*,*,i,1], /current, x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1, arrow_style = 0,layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i+1])
        ax =h2.axes
        ax[0].hide = 1
        ax[1].hide = 1
      endelse
    endfor
    end     
3:  begin
    x = findgen(N[0])
    y = findgen(N[1])
    for i = 0,N[2]-1 do begin 
      if ioverlay eq 1 then begin
        h1 = image((*V)[*,*,i,4], aspect_ratio = dY/dX, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i+1],/current)
        h2 = streamline(reform((*V)[*,*,i,0],N[0],N[1]), reform((*V)[*,*,i,1],N[0],N[1]), x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1,thick = 3,/current,/overplot, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i+1])  
      endif else begin
        h2 = streamline(reform((*V)[*,*,i,0],N[0],N[1]), reform((*V)[*,*,i,1],N[0],N[1]), x,y,aspect_ratio = dY/dX, RGB_TABLE=33,AUTO_COLOR=1,thick = 3,/current, layout=[1+floor(N[2]/floor(sqrt(N[2]))),floor(sqrt(N[2])),i+1])
        ax =h2.axes
        ax[0].hide = 1
        ax[1].hide = 1
      endelse     
    endfor
    end  
4: ivolume,(*V)[*,*,*,vel_type],RENDER_QUALITY = 2, composite_function=1, rgb_table0 =33,/ZERO_OPACITY_SKIP,RENDER_EXTENTS = 1,/INTERPOLATE    
5: flow_orient_map, V, venc
6: mas_mip,(*V)[*,*,*,3], 0
7: mas_mip,(*V)[*,*,*,3], 1

endcase

if ptr_valid(V) then ptr_free,V
close,/all

end

; Subroutine name: flow_movie
; Created by: Magdoom Kulam
; Calling Information:
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display the velocity maps as a movie
;
; Editing Information:

pro flow_movie

common scan_data
common common_widgets

widget_control, vel_button_grp, get_value = vel_type
widget_control, vis_button_grp, get_value = vis_type
CI = project.ci
V = load_state2_data_flow(project.dataarray[CI].flow)
if ~ptr_valid(V) then return 

widget_control,vencfield, get_value = vencstring
venc = (strsplit(vencstring,/extract))[0]
if vel_type eq 3 then range = [0, sqrt(3)*venc] else range = [-venc, +venc]

N = size(*V,/Dimensions)

top_base = widget_base(title = vel_button[vel_type]+'-Velocity',/grid_layout)
if vel_type eq 4 then anim = cw_animate(top_base,2*N[0],2*N[1],N[2], /NO_KILL) else anim = cw_animate(top_base,N[0],N[1],N[2], /NO_KILL)
if vis_type eq 0 or vis_type eq 7 then widget_control, /realize, top_base

case vis_type of
0 : begin
    DEVICE, DECOMPOSED = 0
    loadCT,13
    for i=0,N[2]-1 do begin
        case vel_type of   
        3   : TV, bytscl((*V)[*,*,i,vel_type], min = range[0], max = range[1])       
        4   :  begin
               TV, bytscl((*V)[*,*,i,2], min = -venc, max = +venc),0,0
               TV, bytscl((*V)[*,*,i,3], min = 0, max = sqrt(3)*venc),N[0],0
               TV, bytscl((*V)[*,*,i,0], min = -venc, max = +venc),0,N[1]
               TV, bytscl((*V)[*,*,i,1], min = -venc, max = +venc),N[0],N[1] 
               end
        5    : begin
               rotate_Direction = project.procpramArray[CI].rotate_Direction
               if rotate_Direction eq 0 or rotate_Direction eq 2 then begin
                  dX = project.imndarray[CI].f_voxsz     ; X resolution (mm)
                  dY = project.imndarray[CI].p_voxsz     ; Y resolution (mm)
               endif else begin
                  dX = project.imndarray[CI].p_voxsz     ; X resolution (mm)
                  dY = project.imndarray[CI].f_voxsz     ; Y resolution (mm)
               endelse
               dZ = project.imndarray[CI].s_voxsz         ; Z resolution (mm)
               flow_rate = abs(100*0.06*((*V)[*,*,i,0]*dY*dZ + (*V)[*,*,i,1]*dX*dZ + (*V)[*,*,i,2]*dX*dY)/(dX*dY*dZ))
               TV, bytscl(flow_rate, min = 0)
               end
        else : TV, bytscl((*V)[*,*,i,vel_type], min = range[0], max = range[1])    
        endcase
        cw_animate_load, anim, frame = i, Window = !D.Window
    endfor
    loadCT,0  
    end
else : return
endcase    

cw_animate_run, anim,10
XMANAGER, 'CW_ANIMATE Demo', top_base, EVENT_HANDLER = 'anim_ehandler'
if ptr_valid(V) then ptr_free,V
close,/all

end

; Subroutine name: anim_ehandler
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Dummy event handler for flow animation (required for the movie to loop indefinitely)
;
; Editing Information:

PRO anim_ehandler, EV
WIDGET_CONTROL, /DESTROY, EV.TOP  
END

pro export_velocity

 common scan_data

 CI = project.ci
 V = load_state2_data_flow(project.dataarray[CI].flow)
 if ~ptr_valid(V) then return
 Vel = (*V)[*,*,*,0:2]
 Vreformed = transpose(Vel,[3,0,1,2])
  
 fpath = DIALOG_PICKFILE(Title = 'Select Directory') 
 OpenW, lun,fpath + '.raw', /get_lun
 WriteU,lun,Vreformed
 close, lun
 free_lun, lun
 
 topbase = widget_base(title = 'File export information' ,tlb_frame_attr=1,/grid_layout, frame = 1)
 mainbase = widget_base(topbase)
 rowtxt1 = 'Velocity vectors exported as a binary file to ' + fpath + '.raw'
 rowtxt2 = 'Data type       : ' + strtrim(typename(Vreformed),2)
 rowtxt3 = 'Data order      : U, V, W, X, Y, Z'
 rowtxt4 = 'File dimensions :' + strcompress(strjoin(size(Vreformed,/DIMENSIONS)))
 txt = [[rowtxt1],[rowtxt2],[rowtxt3],[rowtxt4]]
 rowlabel = widget_text(mainbase, value = txt, ysize = 5)
 widget_control, topbase, /realize
 
end

; Subroutine name: flow_gui_cleanup
; Created by: Magdoom Kulam
; Calling Information:
; 
; top_base
;  
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Cleanup procedure before exiting the flow GUI
; 
; Editing Information:

pro flow_gui_cleanup, topbase

    widget_control, topbase, get_uvalue=state
    if (ptr_valid(state)) then ptr_free, state    
    return
   
end