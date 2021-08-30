
; Subroutine name: mas_dce_processing
; Created by: Garrett Astary
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Creates widget base and widgets for implementing a signal to concentration conversion of DCE images

; Editing Information:
; Modified by Magdoom

pro mas_dce_processing

  COMMON scan_data
  COMMON common_widgets
  COMMON dce_data, Dose, m1, m2, a1, a2, aif_flag, aifroi_set,aifroi_name, aifroi_caudate_set, aifroi_caudate_name, $
      plroi_set, plroi_name, r1, r2, vtr_num, dce_num, mask_num, ktransroi_set, ktransroi_name, Ns_caudate, Ns_vessel, $
       t1_fit_type, conc_NM, mask_thresh, n_Pst, n_skip, smooth_ktrans_flag, plot_type, dce_model_fit_flag
    
  ;Population averaged AIF parameters from Bankson paper (Ng et al. (2015) Dependence of DCE-MRI biomarker values on analysis algorithm. PLoS ONE 10(7): e0130168)
  Dose = double(0.24)
  m1 = double(.0333)
  m2 = double(.001)
  a1 = double(3.2)
  a2 = double(2.1)
  r1 = double(2.8)
  r2 = double(4.0)
  
  CI = project.CI
  Ns = project.imndarray[CI].sdim;
  ;creates top base all other widgets are built on
  dce_window_base = widget_base(column=2, TITLE='DCE Processing', xoffset=420)

  ;creating parameter base that will hold widgets for user to input parameters
  dce_parameters_base = widget_base(dce_window_base, /align_left, /base_align_center,frame=2, row=15)

  ;the T1 base holds the checkbox to enable the user to manually enter the T10 into dce_t_one_text
  ;this option may be removed, as of right now it is only useful for hydrogel infusion experiments where
  ;T10 is assumed to be uniform throughout the region of interest

  dce_textbase = widget_base(dce_parameters_base, row=1)
  dce_r1_label = widget_label(dce_textbase, value='r1 (1/mM.s)')
  dce_r1_text = widget_text(dce_textbase, value=strtrim(r1,2), xsize=5, /editable)
  dce_r2_label = widget_label(dce_textbase, value='r2 (1/mM.s)')
  dce_r2_text = widget_text(dce_textbase, value=strtrim(r2,2), xsize=5, /editable)

  dce_filepath_base = widget_base(dce_parameters_base, /base_align_left, row=1)
  dce_droplist_label=widget_label(dce_filepath_base, value='Select DCE-MRI Data:')

  NI = project.ni
  scan_names = strarr(NI)
  scan_names_mask = strarr(NI+1)
  FOR scan=0,NI-1 DO BEGIN
    scan_name_temp = project.scan_list[scan]
    scan_name_split = strsplit(scan_name_temp,PATH_SEP(),/EXTRACT)
    sz_split = size(scan_name_split)
    scan_name = PATH_SEP() + scan_name_split[sz_split[1]-2] + PATH_SEP() + scan_name_split[sz_split[1]-1] + PATH_SEP()
    scan_names[scan] = scan_name
  ENDFOR
  scan_names_mask = [' ',scan_names]
  kep_droplist_dce = widget_droplist(dce_filepath_base, value=scan_names, uvalue = dce_droplist)
  ;builds a base, label and droplist widget for selecting the appropriate VTR/VTI data set from the scan list
  ;This requires the user to load the VTR/VTI data set before running the dce processing package
  kep_filepath_base = widget_base(dce_parameters_base, /base_align_left, row=1)
  kep_droplist_label=widget_label(kep_filepath_base, value='Select T1 Data:')
  kep_droplist_vtr = widget_droplist(kep_filepath_base, value=scan_names, uvalue = vtr_droplist)

  mask_filepath_base = widget_base(dce_parameters_base, /base_align_left, row=1)
  mask_droplist_label= widget_label(mask_filepath_base, value='Select Mask:')
  kep_droplist_mask = widget_droplist(mask_filepath_base, value=scan_names_mask, uvalue = mask_droplist)

  ktransroi_set_base = widget_base(dce_parameters_base, row=1)
  ktransroi_set_label = widget_label(ktransroi_set_base, value = 'Select Ktrans ROI Set')
  ktransroi_set_droplist = widget_droplist(ktransroi_set_base,value = *project.roi.pdisplay_names, uname='ktransroi_set_droplist')
  ktransroi_label = widget_label(ktransroi_set_base, value='Ktrans ROI')
  ktransroi_droplist = widget_droplist(ktransroi_set_base,value='      ', uvalue=' ', sensitive=0, /dynamic_resize)

  ;A AIF ROI will be used to calculate subject specific AIF
  ;The code is designed to allow the user to select a previously defined ROI as the AIF ROI
  ;The user must select the ROI set their ROI is located in and then select the appropriate ROI

  dce_aif_base = widget_base(dce_parameters_base, row=3)
  aifroi_label_base = widget_base(dce_aif_base, column=4)  
  vessroi = widget_button(aifroi_label_base, value='Display vessel ROI',ysize=30)
  caudateroi = widget_button(aifroi_label_base, value='Display caudate ROI',ysize=30)
  aifroi_calc = widget_button(aifroi_label_base, value='AIF from ROI',ysize=30)
  aifroi_calc_data = widget_button(aifroi_label_base, value='AIF from Population',ysize=30)
  aifroi_set_base = widget_base(dce_aif_base, row=1)
  aifroi_set_label = widget_label(aifroi_set_base, value = 'Select Vessel ROI Set')
  aifroi_set_droplist = widget_droplist(aifroi_set_base,value = *project.roi.pdisplay_names, uname='aifroi_set_droplist')
  aifroi_label = widget_label(aifroi_set_base, value='Vessel ROI')
  aifroi_droplist = widget_droplist(aifroi_set_base,value='      ', uvalue=' ', /dynamic_resize)
  aifroi_slice = widget_label(aifroi_set_base, value='Slice #')
  aifroi_slice_droplist = widget_droplist(aifroi_set_base,value=strtrim(indgen(Ns),2), uvalue=strtrim(indgen(Ns),2))
  aifroi_set_caudate_base = widget_base(dce_aif_base, row=1)
  aifroi_set_caudate_label = widget_label(aifroi_set_caudate_base, value = 'Select Caudate ROI Set')
  aifroi_set_caudate_droplist = widget_droplist(aifroi_set_caudate_base,value = *project.roi.pdisplay_names, uname='aifroi_set_caudate_droplist')
  aifroi_caudate_label = widget_label(aifroi_set_caudate_base, value='Caudate ROI')
  aifroi_caudate_droplist = widget_droplist(aifroi_set_caudate_base,value='      ', uvalue=' ', /dynamic_resize)
  aifroi_caudate_slice = widget_label(aifroi_set_caudate_base, value='Slice #')
  aifroi_caudate_slice_droplist = widget_droplist(aifroi_set_caudate_base,value=strtrim(indgen(Ns),2), uvalue=strtrim(indgen(Ns),2))

  fit_thresh_base = widget_base(dce_parameters_base, row=1)
  fit_thresh_label = widget_label(fit_thresh_base, value='Fitting Threshold')
  fit_thresh_text = widget_text(fit_thresh_base, value='5.0', xsize=6, /editable)
  fit_thresh_label2 = widget_label(fit_thresh_base, value='(% times max signal)')

  num_pst_base = widget_base(dce_parameters_base, row = 1)
  num_pst_label = widget_label(num_pst_base, value = 'Select Post Scans to Fit')
  num_pst_slider = widget_slider(num_pst_base, minimum = 1, maximum=2, value=1)

  smooth_ktrans_base = widget_base(num_pst_base, row =1,/nonexclusive)
  smooth_ktrans_enable = widget_button(smooth_ktrans_base, value='Smooth time series data')

  ;creating functions base that will hold widgets that actually execute the desired functions
  dce_functions_base = widget_base(dce_window_base, /align_right, /base_align_center, frame=2, row=8)

  ;pixel-by-pixel tracer concentration calculation
  cc_buttons_base = widget_base(dce_functions_base, /base_align_center,row=2)
  methods = ["Newton's Method", 'Assumption Method']
  conc_button_grp = cw_bgroup(cc_buttons_base, methods, exclusive = 1, column = 1, button_uvalue = methods, label_top = 'Concentration Calculation', frame = 2,/no_release , set_value=1)
  cc_display_base = widget_base(cc_buttons_base, /base_align_center, row=1)
  sl_display_button = widget_button(cc_display_base, value='Display slice')
  cc_display_button = widget_button(cc_display_base, value='Display [CA]')
  ; cc_curves_button = widget_button(cc_display_base,value='Print C(t)')

  dce_pl_label_base = widget_base(dce_functions_base, row=2)
  plots = ['Patlak Plot', 'Logan Plot', 'Extended-Patlak Plot']
  pl_label = widget_label(dce_pl_label_base, value='Graphical Analysis')
  dce_plroi_base = widget_base(dce_functions_base, row=3, frame = 1)
  plroi_set_label = widget_label(dce_plroi_base, value = 'Select ROI Set')
  plroi_set_droplist = widget_droplist(dce_plroi_base,value = *project.roi.pdisplay_names, uname='plroi_set_droplist')
  plroi_display = widget_button(dce_plroi_base, value='Display ROI',ysize=30)
  plroi_label = widget_label(dce_plroi_base, value='ROI')
  plroi_droplist = widget_droplist(dce_plroi_base,value='      ', uvalue=' ', sensitive=0, /dynamic_resize)
  plot_droplist = widget_droplist(dce_plroi_base,value = plots, uname='plot_set_droplist')
  torigin_base = widget_base(dce_plroi_base, row = 1)
  torigin_label = widget_label(torigin_base, value = 'Scans to Skip')
  torigin_slider = widget_slider(torigin_base, minimum = 0, maximum=2, value=0)
  pl_display_button = widget_button(dce_plroi_base, value='Display')

  aif_parameters_base = widget_base(dce_functions_base, row=4, frame = 1)
  aif_label = widget_label(aif_parameters_base, value='AIF Parameters')
  aif_m_base = widget_base(aif_parameters_base, row=1)
  kep_m1_label = widget_label(aif_m_base, value='m1 (1/min):')
  kep_m1_text = widget_text(aif_m_base, value=strtrim(m1,2), xsize=5, /editable)
  kep_m1_base = widget_base(aif_m_base, row=1)
  kep_m2_label = widget_label(aif_m_base, value='m2 (1/min):')
  kep_m2_text = widget_text(aif_m_base, value=strtrim(m2,2), xsize=5, /editable)
  aif_a_base = widget_base(aif_parameters_base, row=1)
  kep_a1_label = widget_label(aif_a_base, value='a1 (kg/L):')
  kep_a1_text = widget_text(aif_a_base, value=strtrim(a1,2), xsize=5, /editable)
  kep_a2_label = widget_label(aif_a_base, value='a2 (kg/L):')
  kep_a2_text = widget_text(aif_a_base, value=strtrim(a2,2), xsize=5, /editable)
  aif_dose_base = widget_base(aif_parameters_base, row=1)
  kep_dose_label = widget_label(aif_dose_base, value='Dose (mmol/kg.bw):')
  kep_dose_text = widget_text(aif_dose_base, value=strtrim(Dose,2), xsize=5, /editable)
  
  dce_kep_label_base =  widget_base(dce_functions_base, row = 2)
  kep_button = ['Model','Ktrans','EES volume fraction','Plasma volume fraction']
  kep_button_grp = cw_bgroup(dce_kep_label_base,kep_button, column = 2, exclusive = 1, button_uvalue = kep_button, label_top = 'Model Fit', frame = 2,/no_release , set_value=0)
  kep_display_base = widget_base(dce_kep_label_base, /base_align_center, row=1)
  kep_display_button = widget_button(dce_kep_label_base, value='Display')
  kep_export_button = widget_button(dce_kep_label_base, value='Export')
  
  aif_flag = long(1)
  aifroi_set = long(0)
  aifroi_name = long(0) 
  aifroi_caudate_set = long(0)
  aifroi_caudate_name = long(0)
  plroi_set = long(0)
  plroi_name = long(0)
  ktransroi_set = long(0)
  ktransroi_name = long(0)
  
  vtr_num = long(CI+1)
  dce_num = long(CI)
  mask_num = long(0)
  Ns_caudate = long(9)
  Ns_vessel = long(1)
  t1_fit_type = long(1)
  dce_model_fit_flag = long(0)
  conc_NM = 'Assumption Method'
  mask_thresh = double(0.05) 
  n_Pst = long(project.imndArray[CI].adim - project.imndArray[CI].n_Pre)
  n_skip = 0
  smooth_ktrans_flag = long(0)
  plot_type = long(0)
    
  widget_control, dce_window_base, /realize

  xmanager, 'mas_dce_processing',dce_window_base,/NO_BLOCK

end

; Subroutine name: time_array
; Created by: Magdoom Kulam
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Time array for DCE data

; Editing Information:
function time_array

  common scan_data
  common dce_data

  n_Pre = project.imndArray[dce_num].n_Pre
  time = *project.imndArray[dce_num].time.array
  time_length = *project.imndArray[dce_num].time.length
  time_zero = time[n_Pre] - time_length[n_Pre]/2
  time -= time_zero
  X = dblarr(n_Pst)
  X = time[n_Pre:n_Pst+n_Pre-1]/60

  return, X
end

; Subroutine name: AIF_function
; Created by: Magdoom Kulam
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Bi-exponential AIF function

; Editing Information:
function AIF_function,p,X=x,Y=y,ERR=err

  COMMON dce_data
  model = Dose*(p[0]*exp(-p[1]*X)+p[2]*exp(-p[3]*X))
  return, (y-model)/err

end

; Subroutine name: linear_fit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Linear fit model for Patlak and Logan plots
;
; Editing Information:
function linear_fit,p,X=x,Y=y,ERR=err

  model = p[0]*X+p[1]
  return, (y-model)/err

end
; Subroutine name: tofts_model
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model3 for DCE data (Extended Tofts Model)
;
; Editing Information:
function tofts_model,p,X=x,Y=y,ERR=err

  t = time_Array()
  S = 1/(t[1]-t[0])
  K = make_array(n_elements(t),value = 1.0,/double) 
  model = p[0]*convol(X,K,S,center = 0,/edge_zero) + p[1]*X

  return, (y-model)/err

end

; Subroutine name: std_model
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model3 for DCE data (Extended Tofts Model)
;
; Editing Information:
function std_model,p,X=x,Y=y,ERR=err

  t = time_Array()
  S = 1/(t[1]-t[0])
  model = p[0]*convol(X,exp(-p[1]*t),S,center = 0,/edge_zero) + p[2]*X

  return, (y-model)/err

end

; Subroutine name: dce_calculate_T10
; Created by: Garrett Astary
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Will perform pixel-by-pixel T10 calculation for desired slice and adim. T10 map will be used in
; concentration calculation.

; Editing Information:
pro dce_calculate_T10, p_T10arr, mask = mask, slice = slice
	COMMON scan_data
  COMMON dce_data
  
	s_T10 = size(p_T10arr)
	
	if ~n_elements(slice) then slice = project.procPramArray[vtr_num].sdim_start
	
	if t1_fit_type EQ 0 then begin
	 if isa(mask) then mas_cfit2_IMG_FIT_T1_VarTR, fit_result = fit_result,threshold_pct = 0.1, for_slice = slice, /progressbar, mask = mask, ci = vtr_num $
	   else mas_cfit2_IMG_FIT_T1_VarTR, fit_result = fit_result,threshold_pct = 0.1, for_slice = slice, /progressbar, ci = vtr_num
	 *p_T10arr = (*fit_result).param_maps[*,*,1]
	endif else begin
	if isa(mask) then mas_cfit2_IMG_FIT_T1_VarTI_POWELL, fit_result = fit_result,threshold_pct = 0.1, for_slice = slice, /progressbar, mask = mask, ci = vtr_num $
	 else mas_cfit2_IMG_FIT_T1_VarTI_POWELL, fit_result = fit_result,threshold_pct = 0.1, for_slice = slice, /progressbar, ci = vtr_num
	*p_T10arr = (*fit_result).param_maps[*,*,1]
	endelse 
	
end

; Subroutine name: dce_calculate_conc
; Created by: Garrett Astary
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Will perform pixel-by-pixel T10 calculation for desired slice and adim. T10 map will be used in
; concentration calculation.

; Editing Information:
; Magdoom 
; - Included calculation for spoiled gradient echo (SPGR) (05/03/15)
; - Updated convergence criteria for Newton's method from niterations = 20 to residual < 0.1 (06/01/15)

pro dce_calculate_conc, output = output, adim = adim, T10_map = T10_map, slice = slice, mask = mask

  COMMON scan_data
  COMMON common_widgets
  common dce_data
  
  fdim = project.ImndArray[dce_num].fdim*project.procPramArray[dce_num].Freq_interp
  pdim = project.ImndArray[dce_num].pdim*project.procPramArray[dce_num].Phase_interp
  alpha = project.ImndArray[dce_num].alpha*!pi/180
  n_Pre = project.imndArray[dce_num].n_Pre
  if ~n_elements(adim) then adim = project.procPramArray[dce_num].adim_start
  if ~n_elements(slice) then slice = project.procPramArray[dce_num].sdim_start
  if ~n_elements(mask) then mask = make_array(fdim,pdim, value = 1.0)
  
  IF ~n_elements(T10_map) then begin   
      p_T10arr = ptr_new(dblarr(fdim,pdim))
      dce_calculate_T10, p_T10arr
      T10_map = *p_T10arr  
  ENDIF
  neg_index = where(T10_map le 0)
  nan_index = where(~finite(T10_map))
  T10_map[neg_index] = 0
  T10_map[nan_index] = 0
  
  ;calculate value of pre-images to get S(0)
  pre_avg = mask*reform(mean((*project.dataArray[dce_num].state1)[*,*,slice,0:n_Pre-1],dimension = 4))
  p_image = mask*(*project.dataArray[dce_num].state1)[*,*,slice,adim]
  ;calculate the concentration pixel-by-pixel

  TR = project.imndArray[dce_num].recov_time/1000.0
  TE = project.imndArray[dce_num].echo_time/1000.0
  conc_array = dblarr(fdim,pdim)
  
  case conc_NM of
   "Newton's Method" : begin
                       FOR xx=0, fdim-1 DO BEGIN
                          FOR yy=0, pdim-1 DO BEGIN
                              if p_image[xx,yy] ne 0 then begin
                                  k1 = -TE*r2
                                  k2 = -TR/T10_map[xx,yy]
                                  k3 = TR*r1
                                  S0 = pre_avg[xx,yy]
                                  SC = p_image[xx,yy]
                                  K = (SC/S0)*(1-exp(-TR/T10_map[xx,yy]))
                                  F = (S0*(1-cos(alpha)*exp(k2))-SC*cos(alpha)*(1-exp(k2)))/(S0*(1-cos(alpha)*exp(k2)) - SC*(1-exp(k2)))
                                  x = (alog(F)/TR - 1/T10_map[xx,yy])/r1 
                                  if x lt 0 or ~finite(x) then x = 0               
                                  x0 = x
                                  maxiter = 20
                                  residual = 1
                                  niter = 0
                                  while residual gt 0.1 do begin
                                      f = (1-cos(alpha)*exp(k2))*(exp(k1*x) - exp(k2)*exp(x*(k1-k3))) - K*(1-cos(alpha)*exp(k2)*exp(-k3*x))
                                      fp = (1-cos(alpha)*exp(k2))*(k1*exp(k1*x) - (k1-k3)*exp(k2)*exp(x*(k1-k3)))-K*k3*cos(alpha)*exp(k2)*exp(-k3*x)
                                      x1 = x0 - f/fp
                                      residual = abs(x1-x0)
                                      x0 = x1
                                      niter++
                                      if niter gt maxiter then break
                                 endwhile
                                 if residual le 0.1 then conc_array[xx,yy] = x1 else conc_array[xx,yy] = x
                                 endif
                             ENDFOR
                          ENDFOR
                          end
    'Assumption Method' : begin
                          F = (pre_avg*(1-cos(alpha)*exp(-TR/T10_map))-p_image*cos(alpha)*(1-exp(-TR/T10_map)))/(pre_avg*(1-cos(alpha)*exp(-TR/T10_map)) - p_image*(1-exp(-TR/T10_map)))
                          conc_array = (alog(F)/TR - 1/T10_map)/r1                          
                          end
  endcase
  neg_index = where(conc_array le 0)
  nan_index = where(~finite(conc_array))
  conc_array[neg_index] = 0
  conc_array[nan_index] = 0
  output = conc_array
  
end

; Subroutine name: mas_dce_processing_event
; Created by: Garrett Astary
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Handles button and text field events generated by GUI.

; Editing Information:


pro mas_dce_processing_event, event

 ;bring in global variables
    COMMON scan_data
    COMMON common_widgets
    COMMON dce_data
   
	CASE event.id OF

	    dce_r1_text: begin
	                 widget_control, dce_r1_text, get_value = r1
	                 r1 = (double(r1))[0]
	                 dce_model_fit_flag = 0
	                 end

	    dce_r2_text: begin
	                 widget_control, dce_r2_text, get_value= r2
	                 r2 = (double(r2))[0]
	                 dce_model_fit_flag = 0
	                 end
	     
      kep_droplist_dce: BEGIN
                        dce_num = widget_info(kep_droplist_dce, /DROPLIST_SELECT)
                        n_Pst = project.imndArray[dce_num].adim - project.imndArray[dce_num].n_Pre
                        widget_control, num_pst_slider, set_slider_max = n_Pst
                        widget_control, num_pst_slider, set_value = n_Pst
                        widget_control, torigin_slider, set_slider_max = n_Pst
                        dce_model_fit_flag = 0
                        END  
	    
	  kep_droplist_vtr: BEGIN
	                    vtr_num = widget_info(kep_droplist_vtr, /DROPLIST_SELECT)
                      if ptr_valid(project.imndArray[vtr_num].inversion_time_ptr) then t1_fit_type = 1 else t1_fit_type = 0
                      dce_model_fit_flag = 0
                      END
      
      kep_droplist_mask: BEGIN
                         widget_control, kep_droplist_mask, get_value = val
                         mask_index = widget_info(kep_droplist_mask, /DROPLIST_SELECT)
                         if strcmp(val[mask_index],' ') then mask_num = -1 else mask_num = mask_index-1
                         END
	    
     fit_thresh_text: BEGIN
	                    widget_control, fit_thresh_text, get_value=fit_thresh
                      mask_thresh = fit_thresh/100.0
                      dce_model_fit_flag = 0 
                      END

     num_pst_slider: begin
                     widget_control, num_pst_slider, get_value = n_Pst
                     dce_model_fit_flag = 0
                     end
                     
    torigin_slider: begin
                    widget_control, torigin_slider, get_value = n_skip
                    dce_model_fit_flag = 0
                    end                

	  aifroi_set_droplist: BEGIN
	                       ;we check to see which item in droplist is selected
                         aifroi_set = widget_info(aifroi_set_droplist, /DROPLIST_SELECT)
                         IF aifroi_set EQ 0 THEN widget_control, aifroi_droplist,set_value='      ', sensitive=0 $
                         ELSE BEGIN
                         ;collect ROI names in the selected set
                         rois = (*project.roi.pRois[aifroi_set])
                         IF     (size(rois))[0] EQ 0 THEN num_roi = 1 $
                         ELSE num_roi =  (size(rois))[1]
                         roi_names = strarr(num_roi)
                         FOR temp=0 , num_roi-1 DO BEGIN
                              rois[temp] -> GETPROPERTY, name= name
                              roi_names[temp] = name
                         END
                         widget_control, aifroi_droplist, set_value=roi_names, sensitive= 1
                         ENDELSE
                         dce_model_fit_flag = 0
                         END
	    
	  aifroi_droplist: begin
	                   aifroi_name = widget_info(aifroi_droplist, /DROPLIST_SELECT)
	                   dce_model_fit_flag = 0
	                   end

	  aifroi_set_caudate_droplist: BEGIN
                        	      ;we check to see which item in droplist is selected
                        	      aifroi_caudate_set = widget_info(aifroi_set_caudate_droplist, /DROPLIST_SELECT)
                        	      IF aifroi_caudate_set EQ 0 THEN widget_control, aifroi_caudate_droplist,set_value='      ', sensitive=0 $
                        	      ELSE BEGIN
                        	        ;collect ROI names in the selected set
                        	        rois = (*project.roi.pRois[aifroi_caudate_set])
                        	        IF     (size(rois))[0] EQ 0 THEN num_roi = 1 $
                        	        ELSE num_roi =  (size(rois))[1]
                        	        roi_names = strarr(num_roi)
                        	        FOR temp=0 , num_roi-1 DO BEGIN
                        	          rois[temp] -> GETPROPERTY, name= name
                        	          roi_names[temp] = name
                        	        END
                        	        widget_control, aifroi_caudate_droplist, set_value=roi_names, sensitive= 1
                        	      ENDELSE
                                dce_model_fit_flag = 0
                        	    END

	    aifroi_caudate_droplist: begin
	                             aifroi_caudate_name = widget_info(aifroi_caudate_droplist, /DROPLIST_SELECT)
	                             dce_model_fit_flag = 0
	                             end
	    
	    aifroi_caudate_slice_droplist: begin
	                                   Ns_caudate = event.index
	                                   dce_model_fit_flag = 0
	                                   end
	                                   
	    aifroi_slice_droplist : begin
	                            Ns_vessel = event.index
	                            dce_model_fit_flag = 0
	                            end
	                            
	    plot_droplist: plot_type = widget_info(plot_droplist, /DROPLIST_SELECT) 
	 
	    plroi_set_droplist: BEGIN
                  	      ;we check to see which item in droplist is selected
                  	      plroi_set = widget_info(plroi_set_droplist, /DROPLIST_SELECT)
                  	      IF plroi_set EQ 0 THEN widget_control, plroi_droplist,set_value='      ', sensitive=0 $
                  	      ELSE BEGIN
                  	        ;collect ROI names in the selected set
                  	        rois = (*project.roi.pRois[plroi_set])
                  	        IF     (size(rois))[0] EQ 0 THEN num_roi = 1 $
                  	        ELSE num_roi =  (size(rois))[1]
                  	        roi_names = strarr(num_roi)
                  	        FOR temp=0 , num_roi-1 DO BEGIN
                  	          rois[temp] -> GETPROPERTY, name= name
                  	          roi_names[temp] = name
                  	        END
                  	        widget_control, plroi_droplist, set_value=roi_names, sensitive= 1
                  	      ENDELSE
                  
                  	    END

	   plroi_droplist: plroi_name = widget_info(plroi_droplist, /DROPLIST_SELECT)	
	    
	   ktransroi_set_droplist: BEGIN
                    	      ;we check to see which item in droplist is selected
                    	      ktransroi_set = widget_info(ktransroi_set_droplist, /DROPLIST_SELECT)
                    	      IF ktransroi_set EQ 0 THEN widget_control, ktransroi_droplist,set_value='      ', sensitive=0 $
                    	      ELSE BEGIN
                    	        ;collect ROI names in the selected set
                    	        rois = (*project.roi.pRois[ktransroi_set])
                    	        IF     (size(rois))[0] EQ 0 THEN num_roi = 1 $
                    	        ELSE num_roi =  (size(rois))[1]
                    	        roi_names = strarr(num_roi)
                    	        FOR temp=0 , num_roi-1 DO BEGIN
                    	          rois[temp] -> GETPROPERTY, name= name
                    	          roi_names[temp] = name
                    	        END
                    	        widget_control, ktransroi_droplist, set_value=roi_names, sensitive= 1
                    	      ENDELSE
                    
                    	    END

	  ktransroi_droplist: ktransroi_name = widget_info(ktransroi_droplist, /DROPLIST_SELECT)
	     
    aifroi_calc: dce_calculate_AIF
    
    aifroi_calc_data: dce_calculate_AIF, flag = 1
    
    vessroi   : begin
              fdim = project.ImndArray[dce_num].fdim*project.procPramArray[dce_num].Freq_interp
              pdim = project.ImndArray[dce_num].pdim*project.procPramArray[dce_num].Phase_interp
              if ptr_valid(project.roi.pROIs[aifroi_set]) then aifroimask = (*project.roi.pROIs[aifroi_set])[aifroi_name] -> ComputeMask(Dimensions= [fdim,pdim]) $
              else aifroimask = make_array(fdim,pdim,value = 1.0)
              curve = roberts(aifroimask)
              array = project.procPramArray[dce_num].adim_start
              mask_bet = (*project.dataArray[mask_num].state1)[*,*,Ns_vessel]              
              im = (*project.dataArray[dce_num].state1)[*,*,Ns_vessel,array]
              ind = where(mask_bet ne 0)
              rc = array_indices(mask_bet,ind)
              bb = [(rc[where(rc eq min(rc[0,*]))])[0], (rc[where(rc eq max(rc[0,*]))])[0], (rc[where(rc eq min(rc[1,*]))])[0], (rc[where(rc eq max(rc[1,*]))])[0]]
              dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
              resolution*=0.393701
              dims = long([1.75,2]/float(resolution))
              h1 = image(rotate((im*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), DIMENSIONS = dims,transparency = 0, MARGIN = 0,FONT_NAME = 'Times',FONT_STYLE = 'bf', FONT_SIZE = 12, max_value = 0.5)
              h2 = image(rotate(curve[bb[0]:bb[1],bb[2]:bb[3]],1), DIMENSIONS = dims,/overplot,transparency = 50, MARGIN = 0, RGB_TABLE = 13)
              h2.Save, "VesselROI.tiff", RESOLUTION=300
              end

     caudateroi : begin
                  fdim = project.ImndArray[dce_num].fdim*project.procPramArray[dce_num].Freq_interp
                  pdim = project.ImndArray[dce_num].pdim*project.procPramArray[dce_num].Phase_interp
                  if ptr_valid(project.roi.pROIs[aifroi_caudate_set]) then aifroimask_caudate = (*project.roi.pROIs[aifroi_caudate_set])[aifroi_caudate_name] -> ComputeMask(Dimensions= [fdim,pdim])$
                    else aifroimask_caudate = make_array(fdim,pdim,value = 1.0)
                  curve = roberts(aifroimask_caudate) 
                  array = project.procPramArray[dce_num].adim_start
                  mask_bet = (*project.dataArray[mask_num].state1)[*,*,Ns_caudate]
                  im = (*project.dataArray[dce_num].state1)[*,*,Ns_caudate,array]
                  ind = where(mask_bet ne 0)
                  rc = array_indices(mask_bet,ind)
                  bb = [(rc[where(rc eq min(rc[0,*]))])[0], (rc[where(rc eq max(rc[0,*]))])[0], (rc[where(rc eq min(rc[1,*]))])[0], (rc[where(rc eq max(rc[1,*]))])[0]]
                  dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
                  resolution*=0.393701
                  dims = long([1.75,2]/float(resolution))
                  h1 = image(rotate((im*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), DIMENSIONS = dims,transparency = 0, MARGIN = 0,FONT_NAME = 'Times',FONT_STYLE = 'bf', FONT_SIZE = 15)
                  h2 = image(rotate(curve[bb[0]:bb[1],bb[2]:bb[3]],1), DIMENSIONS = dims,/overplot,transparency = 50, MARGIN = 0, RGB_TABLE = 13)
                  h2.Save, "CaudateROI.tiff", RESOLUTION=300
                  end
                  
    plroi_display : begin
                    fdim = project.ImndArray[dce_num].fdim*project.procPramArray[dce_num].Freq_interp
                    pdim = project.ImndArray[dce_num].pdim*project.procPramArray[dce_num].Phase_interp
                    if ptr_valid(project.roi.pROIs[plroi_set]) then mask = (*project.roi.pROIs[plroi_set])[plroi_name] -> ComputeMask(Dimensions= [fdim,pdim]) else $
                    mask = make_array(fdim,pdim,value = 1.0)
                    curve = roberts(mask)
                    array = project.procPramArray[dce_num].adim_start
                    Ns = project.procPramArray[dce_num].sdim_start
                    mask_bet = (*project.dataArray[mask_num].state1)[*,*,Ns]
                    im = (*project.dataArray[dce_num].state1)[*,*,Ns,array]
                    dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)                
                    ind = where(mask_bet ne 0)
                    rc = array_indices(mask_bet,ind)
                    bb = [(rc[where(rc eq min(rc[0,*]))])[0], (rc[where(rc eq max(rc[0,*]))])[0], (rc[where(rc eq min(rc[1,*]))])[0], (rc[where(rc eq max(rc[1,*]))])[0]]
                    resolution*=0.393701
                    dims = long([1.75,2]/float(resolution))
                    h1 = image(rotate((im*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), DIMENSIONS = dims,transparency = 0, MARGIN = 0,FONT_NAME = 'Times',FONT_STYLE = 'bf', FONT_SIZE = 12, max_value = 0.5)
                    h2 = image(rotate(curve[bb[0]:bb[1],bb[2]:bb[3]],1), DIMENSIONS = dims,/overplot,transparency = 50, MARGIN = 0, RGB_TABLE = 13)
                    h2.Save, "GraphROI.tiff", RESOLUTION=300
                  end                        
    pl_display_button: dce_plots
     
	  smooth_ktrans_enable: begin
	                        smooth_ktrans_flag = event.select
	                        dce_model_fit_flag = 0
	                        end
	 
		kep_m1_text: begin
		             widget_control, kep_m1_text, get_value = m1
		             m1 = (double(m1))[0]
		             dce_model_fit_flag = 0
                 end
		           

	  kep_m2_text: begin
	               widget_control, kep_m2_text, get_value=m2
	               m2 = (double(m2))[0]
	               dce_model_fit_flag = 0
                 end
	     
	  kep_a1_text: begin
	               widget_control, kep_a1_text, get_value=a1
	               a1 = (double(a1))[0]
	               dce_model_fit_flag = 0
	               end
	             
	      
	  kep_a2_text: begin
	               widget_control, kep_a2_text, get_value=a2
	               a2 = (double(a2))[0]
	               dce_model_fit_flag = 0
	               end
	             
	      
	  kep_dose_text: begin
	                 widget_control, kep_dose_text, get_value=Dose
	                 Dose = (double(Dose))[0]
	                 dce_model_fit_flag = 0
	                 end
	      
    conc_button_grp: begin
                     conc_NM = event.value
                     dce_model_fit_flag = 0
                     end
    sl_display_button : begin
                        array = project.procPramArray[dce_num].adim_start
                        slice = project.procPramArray[dce_num].sdim_start
                        mask_bet = (*project.dataArray[mask_num].state1)[*,*,slice]
                        im = (*project.dataArray[dce_num].state1)[*,*,slice,array]
                        ind = where(mask_bet ne 0)
                        rc = array_indices(mask_bet,ind)
                        bb = [(rc[where(rc eq min(rc[0,*]))])[0], (rc[where(rc eq max(rc[0,*]))])[0], (rc[where(rc eq min(rc[1,*]))])[0], (rc[where(rc eq max(rc[1,*]))])[0]]
                        dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
                        resolution*=0.393701
                        dims = long([2,2]/float(resolution))
                        h = image(rotate((im*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), dimensions = dims, MARGIN = 0)       
                        h.Save, "Slice.tiff", RESOLUTION=300
                        end                
	  cc_display_button: BEGIN
                       dce_calculate_conc, output=conc                      
                       array = project.procPramArray[dce_num].adim_start
                       slice = project.procPramArray[dce_num].sdim_start
                       mask_bet = (*project.dataArray[mask_num].state1)[*,*,slice]  
                       ind = where(mask_bet ne 0)
                       rc = array_indices(mask_bet,ind)
                       bb = [(rc[where(rc eq min(rc[0,*]))])[0], (rc[where(rc eq max(rc[0,*]))])[0], (rc[where(rc eq min(rc[1,*]))])[0], (rc[where(rc eq max(rc[1,*]))])[0]]
                       dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
                       resolution*=0.393701
                       dims = long([2.5,2]/float(resolution))                       
                       h = image(rotate((conc*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), rgb_table = 33,min_value = 0, max_value = 0.5, dimensions = dims,FONT_SIZE=9,FONT_NAME='Times' , MARGIN = [0.0,0.0,0.3,0.0])
                       c = colorbar(target = h, orientation = 1, title = '[Gd-DTPA] (mM)',FONT_SIZE=8,TEXTPOS = 1,/NORMAL,POSITION = [0.75,0.1,0.8,0.9],FONT_NAME='Times')
                       c.scale,0.65,0.65,1
                       h.Save, "Concentration.tiff", RESOLUTION=600
                       END
	   kep_display_button: BEGIN
                         fdim = project.imndArray[dce_num].fdim
                         pdim = project.imndArray[dce_num].pdim
                         slice = project.procPramArray[dce_num].sdim_start
                         widget_control,kep_button_grp, get_value = model_out_type
                         
                         ; Brain extraction mask
                         if mask_num ne -1 then mask_bet = (*project.dataArray[mask_num].state1)[*,*,slice] else mask_bet = make_array(fdim,pdim,value = 1.0)
                         ; ROI mask
                         if ptr_valid(project.roi.pROIs[ktransroi_set]) then mask_ktrans = (*project.roi.pROIs[ktransroi_set])[ktransroi_name] -> ComputeMask(Dimensions= [fdim,pdim]) else $
                           mask_ktrans = make_array(fdim,pdim,value = 1.0)
                           
                         if dce_model_fit_flag eq 0 then dce_model_selection
                         model_select = (*project.dataarray[dce_num].dce_model_fit)[*,*,0] 
                         ktrans = (*project.dataarray[dce_num].dce_model_fit)[*,*,1] 
                         kep = (*project.dataarray[dce_num].dce_model_fit)[*,*,2] 
                         vp = (*project.dataarray[dce_num].dce_model_fit)[*,*,3] 
                         phi = ktrans/kep
                         phi(where(~finite(phi))) = 0    
                         
                         ind = where(mask_bet ne 0)
                         rc = array_indices(mask_bet,ind)
                         bb = [(rc[where(rc eq min(rc[0,*]))])[0], (rc[where(rc eq max(rc[0,*]))])[0], (rc[where(rc eq min(rc[1,*]))])[0], (rc[where(rc eq max(rc[1,*]))])[0]]
                         dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
                         resolution*=0.393701
                         dims = long([2,2]/float(resolution))
                                   
                         case model_out_type of
                                0 : begin
                                    h = image(rotate((model_select*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), rgb_table = 33, dimensions = dims,FONT_SIZE=9,FONT_NAME='Times' , MARGIN =0)
                                    h.Save, "Model.tiff", RESOLUTION=600
                                    end
                                1 : begin
                                    h = image(1e3*rotate((ktrans*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), rgb_table = 33, min_value = 0, max_value = 25, dimensions = dims, MARGIN = [0.0,0.0,0.3,0.0])
                                    c = colorbar(target = h, orientation = 1, title = '$K^{trans} (10^{-3}min^{-1})$',FONT_SIZE=8, TEXTPOS = 1,/NORMAL,POSITION = [0.75,0.1,0.8,0.9],FONT_NAME='Times')
                                    c.scale,0.65,0.65,1
                                    h.Save, "Ktrans.tiff", RESOLUTION=600                             
                                    end
           
                                2 : begin
                                  h = image(rotate((phi*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), rgb_table = 33, min_value = 0, max_value = 0.1, dimensions = dims,FONT_SIZE=9,FONT_NAME='Times', MARGIN = [0.0,0.0,0.3,0.0])
                                  c = colorbar(target = h, orientation = 1, title = '$v_e$',FONT_SIZE=8, TEXTPOS = 1,/NORMAL,POSITION = [0.75,0.1,0.8,0.9],FONT_NAME='Times')
                                  c.scale,0.65,0.65,1
                                  h.Save, "Ve.tiff", RESOLUTION=600
                                    
                                    end    
                                3  : begin
                                  h = image(rotate((vp*mask_bet)[bb[0]:bb[1],bb[2]:bb[3]],1), min_value = 0, max_value = 0.05, rgb_table = 33, dimensions = dims,FONT_SIZE=9,FONT_NAME='Times', MARGIN = [0.0,0.0,0.3,0.0])
                                  c = colorbar(target = h, orientation = 1, title = '$v_p$',FONT_SIZE=8, TEXTPOS = 1,/NORMAL,POSITION = [0.75,0.1,0.8,0.9],FONT_NAME='Times')
                                  c.scale,0.65,0.65,1
                                  h.Save, "Vp.tiff", RESOLUTION=600                                 
                                     end                               
                         endcase                                                                   
                        END
	kep_export_button : begin
	                    if dce_model_fit_flag eq 0 then dce_model_selection
	                    dce_dat_export
                      end
  else              : break                    
	ENDCASE
	
end

; Subroutine dce_param_file
; Created by: Magdoom Kulam
; Calling Information: 
; Purpose: This function exports a text file with model fit parameters
; Bugs or Important Comments to Developers:

pro dce_param_file, output = output

  COMMON common_widgets
  COMMON scan_data
  COMMON dce_data
  
  CI = project.ci
  fdim = project.ImndArray[CI].fdim*project.procPramArray[CI].Freq_interp
  pdim = project.ImndArray[CI].pdim*project.procPramArray[CI].Phase_interp
  sdim = project.imndArray[CI].sdim
  adim = project.imndArray[CI].adim

  param_string = strarr(21)
  param_string[0] = 'Imaging Parameters' + string(13B) + string(10B)
  param_string[1] = 'Dynamics Data: ' + project.scan_list[CI]
  param_string[2] = 'VTR Data: ' + project.scan_list[vtr_num]
  param_string[3] = 'Freq Interpolation: ' + strtrim(string(project.procPramArray[CI].Freq_interp),2)+ string(13B) + string(10B)
  param_string[4] = 'Phase Interpolation: ' + strtrim(string(project.procPramArray[CI].Phase_interp),2) + string(13B) + string(10B)
  param_string[5] = 'Freq Shift: ' +strtrim(string(project.imndArray[CI].fdim_shift),2) + string(13B) + string(10B)
  param_string[6] = 'Phase Shift: ' +strtrim(string(project.imndArray[CI].pdim_shift),2) + string(13B) + string(10B)
  param_string[7] = ' '  + string(13B) + string(10B)
  param_string[8] = 'Processing Parameters'  + string(13B) + string(10B)
  param_string[9] = 'Fitting Threshold: ' + strtrim(string(mask_thresh*100),2) + '%'  + string(13B) + string(10B)
  param_string[10] = 'Final Post Scan: ' + strtrim(string(n_Pst),2)+ string(13B) + string(10B)
  param_string[11] = 'Time Series Smoothing Enabled: ' +  strtrim(string(smooth_ktrans_flag),2)
  param_string[12] = ' '+ string(13B) + string(10B)
  param_string[13] = 'Fit Parameters'
  IF t1_fit_type EQ 1 THEN BEGIN
    param_string[14] = 'T1 Fit Type: Inversion Recovery'
  ENDIF ELSE BEGIN
    param_string[14] = 'T1 Fit Type: Saturation Recovery'
  ENDELSE

  IF conc_NM EQ 0 THEN BEGIN
    param_string[15] = 'Concentration Calculation Type: Newtons Method'
  ENDIF ELSE IF conc_NM EQ 1 THEN  BEGIN
    param_string[15] = 'Concentration Calculation Type: Assumption Method'
  ENDIF

  param_string[16] = 'm1: ' + strtrim(string(m1),2)  + string(13B) + string(10B)
  param_string[17] = 'm2: ' + strtrim(string(m2),2)
  param_string[18] = 'a1: ' + strtrim(string(a1),2)   + string(13B) + string(10B)
  param_string[19] = 'a2: ' + strtrim(string(a2),2)   + string(13B) + string(10B)
  param_string[20] = 'Dose: ' + strtrim(string(dose),2)  + string(13B) + string(10B)

  output = param_string

end

; Subroutine dce_dat_export
; Created by: Magdoom Kulam
; Calling Information: 
; that holds 7 columns of information: x    y    z   Model   Ktrans   Phi  Vp
; The first 4 columns are in integer format and the last 3 columns are in floating point format.

; Bugs or Important Comments to Developers:


pro dce_dat_export
  
  COMMON dce_data
  common scan_data
  
  fpath = DIALOG_PICKFILE( Title = 'Select Directory')
  OpenW, lun,fpath + '.dat', /get_lun
  
  model_select = (*project.dataarray[dce_num].dce_model_fit)[*,*,0]
  ktrans = (*project.dataarray[dce_num].dce_model_fit)[*,*,1]
  kep = (*project.dataarray[dce_num].dce_model_fit)[*,*,2]
  vp = (*project.dataarray[dce_num].dce_model_fit)[*,*,3]
  phi = ktrans/kep
  phi(where(~finite(phi))) = 0 
  
  sz_ktrans = size(model_select,/dimensions)
  x_ktrans = sz_ktrans(0)
  y_ktrans = sz_ktrans(1)
  z = project.procPramArray[dce_num].sdim_start
  
  FOR x=0,x_ktrans-1 DO BEGIN
      FOR y=0,y_ktrans-1 DO BEGIN
            model_temp = model_select[x,y]
            ktrans_temp = ktrans[x,y]
            phi_temp = phi[x,y]
            vp_temp = vp[x,y]
            temp = [x+1, y+1, z+1, model_temp, ktrans_temp, phi_temp, vp_temp]
            PrintF, lun, format = '(4I6,3F)', temp
       ENDFOR
   ENDFOR   
   close, lun
   free_lun, lun
   
   dce_param_file, output = param_string_array
   OpenW, lun,fpath + '.txt', /get_lun
   PrintF, lun, param_string_array
   close, lun
   free_lun, lun
   
end



; Subroutine name: dce_plots
; Created by: Magdoom Kulam
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Patlak/Logan/Extendend Patlak plots of the chosen ROI

; Editing Information:

pro dce_plots
  common scan_data
  COMMON common_widgets
  COMMON dce_data

  n_Pre = project.imndArray[dce_num].n_Pre
  t = time_array()
  fdim = project.ImndArray[dce_num].fdim*project.procPramArray[dce_num].Freq_interp
  pdim = project.ImndArray[dce_num].pdim*project.procPramArray[dce_num].Phase_interp
  adim = project.ImndArray[dce_num].adim
  
  if ptr_valid(project.roi.pROIs[plroi_set]) then mask = (*project.roi.pROIs[plroi_set])[plroi_name] -> ComputeMask(Dimensions= [fdim,pdim]) else $
    mask = make_array(fdim,pdim,value = 1.0)

  Ct = dblarr(n_Pst)
  Cp = dblarr(n_Pst)
  xaxis = dblarr(n_Pst)
  yaxis = dblarr(n_Pst)
  p_T10arr = ptr_new(dblarr(fdim,pdim))
  dce_calculate_T10, p_T10arr
  T10_map = *p_T10arr
  neg_index = where(T10_map le 0)
  nan_index = where(~finite(T10_map))
  T10_map[neg_index] = 0
  T10_map[nan_index] = 0
  
  progressBar = Obj_New('progressbar', Color='red', Text='Calculating Concentration' ,/NOCANCEL)
  progressBar -> Start
  
  FOR time=n_Pre,n_Pre+n_Pst-1 DO BEGIN
    dce_calculate_conc, output = output, adim = time, mask = mask, T10_map = T10_map
    time_dimension = time - n_Pre
    Cp[time_dimension] = Dose*(a1*exp(-m1*t[time_dimension])+a2*exp(-m2*t[time_dimension]))
    Ct[time_dimension] = mean(output[where(mask gt 0)],/nan)    
    progressBar -> Update, (float(time_dimension)/float(n_Pst))*100.0
  endfor
  progressBar -> Destroy
  
  K = make_array(n_elements(t),value = 1.0,/double)
  S = 1/(t[1]-t[0])
  case plot_type of
    0: begin
       xaxis = convol(Cp,K,S,center = 0,/edge_zero)/Cp 
       yaxis = Ct/Cp   
       xlabel = '$\int_0^t C_p(\tau) d\tau/C_t(t)$'
       ylabel = '$C_t(t)/C_p(t)$'    
       end
   1:  begin
       xaxis = convol(Cp,K,S,center = 0,/edge_zero)/Ct
       yaxis = convol(Ct,K,S,center = 0,/edge_zero)/Ct
       xlabel = '$\int_0^t C_p(\tau) d\tau/C_t(t)$'
       ylabel = '$\int_0^t C_t(\tau) d\tau/C_t(t)$'
       end
   2:  begin
       fargs   = {x:Cp,y:Ct,err:1}
       pinfo    = replicate({value:1,limited:[1,0],limits:[0.0,0.0]},3)
       pinfo[2].limited[1] = 1
       pinfo[2].limits[1] = 1.0
       params  = MPFIT('std_model', functargs = fargs, parinfo=pinfo,bestnorm = chisquared,status=status)
       SS_tot  = n_elements(Cp)*variance(Ct,/double)
       rsquared     = 1 - chisquared/SS_tot
       xaxis = convol(Cp,exp(-params[1]*t),S,center = 0,/edge_zero)/Cp
       yaxis = Ct/Cp
       xlabel = '$\int_0^t C_p(\tau) e^{-k_{ep}(t-\tau)}d\tau/C_t(t)$'
       ylabel = '$C_t(t)/C_p(t)$'
       end
  endcase
  
  fargs2   = {x:xaxis[n_skip:time_dimension],y:yaxis[n_skip:time_dimension],err:1}
  pinfo2    = replicate({value:1.0},2)
  params2  = MPFIT('linear_fit', functargs = fargs2, parinfo=pinfo2,bestnorm = chisquared,status=status)
  SS_tot   = (time_dimension-n_skip)*variance(yaxis[n_skip:time_dimension],/double)
  rsquared = 1 - chisquared/SS_tot
  yfit = params2[0]*xaxis + params2[1]
  
  dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
  resolution*=0.393701
  dims = long([3,4]/float(resolution))
  plot_title = ['Patlak', 'Logan', 'Extended-Patlak']
  h1 = plot(xaxis,yaxis,"b2+",xtitle = xlabel,ytitle = ylabel,FONT_SIZE = 11, DIMENSIONS = dims, FONT_NAME = 'Times', MARGIN = [0.22,0.12,0.05,0.05])
  h2 = plot(xaxis,yfit,"r2",/current,/overplot,FONT_SIZE = 11, DIMENSIONS = dims, FONT_NAME = 'Times', MARGIN = [0.22,0.12,0.05,0.05])
  h1.name = plot_title[plot_type]
  h2.name = 'Linear Fit'
  leg = LEGEND(TARGET=[h1,h2],FONT_SIZE = 11, POSITION = [170,350],/DEVICE)
  if params2[1] gt 0 then sign = '+' else sign = ' '
  t1 = TEXT(70,285, '$y = $'+strtrim(string(1e3*params2[0],format='(F8.2)'),2)+'E-3 x'+sign +strtrim(string(params2[1],format='(F8.2)'),2) , /DEVICE, FONT_SIZE=11, FONT_NAME='Times')
  t2 = TEXT(70,265, '$ r^2 = $'+strtrim(string(rsquared,format='(F5.2)'),2), /DEVICE, FONT_SIZE=11, FONT_NAME='Times')
  h2.Save, "Logan.tiff", RESOLUTION=1200
  end

  ; Subroutine name: AIF_data
  ; Created by: Magdoom Kulam
  ; Calling Information:

  ; Bugs or Important Comments to Developers:


  ; Purpose of subroutine: Time (min) and counts per minute average radiological blood samples in 200-250 g rat

  ; Editing Information:  
  function AIF_data

  data = [[0.033333333,  0.052359205],$
          [0.1 , 3.320762767],$
          [0.166666667,  5.167747002],$
          [0.25,  3.47968972],$
          [0.45 , 2.197912041],$
          [0.6 , 1.685093418],$
          [0.75 , 1.412886425],$
          [1.008333333 , 1.14829854],$
          [2.0 , 0.676451302],$
          [3.0 , 0.554787325],$
          [5.0 , 0.492346101],$
          [7.5 ,  0.379544803],$
          [10.0 , 0.333632819],$
          [15.0 , 0.265640913],$
          [20.0 , 0.218385957]]

  return,data  
  end
  
  ; Subroutine name: dce_calculate_AIF
  ; Created by: Magdoom Kulam
  ; Calling Information:

  ; Bugs or Important Comments to Developers:


  ; Purpose of subroutine: Will perform AIF calculation based on the blood vessel ROI drawn for the selected slice

  ; Editing Information:

pro dce_calculate_AIF, flag = flag

  COMMON scan_data
  COMMON common_widgets
  COMMON dce_data

  fdim = project.ImndArray[dce_num].fdim*project.procPramArray[dce_num].Freq_interp
  pdim = project.ImndArray[dce_num].pdim*project.procPramArray[dce_num].Phase_interp
  alpha = project.ImndArray[dce_num].alpha*!pi/180
  n_Pre = project.imndArray[dce_num].n_Pre
  if Dose eq 0 then begin
    msg = dialog_message('Enter a valid dose',/error)
    return
  endif
  
  if ~n_elements(flag) then begin
    if ptr_valid(project.roi.pROIs[aifroi_set]) then aifroimask = (*project.roi.pROIs[aifroi_set])[aifroi_name] -> ComputeMask(Dimensions= [fdim,pdim]) $
    else begin
      msg = dialog_message('Select a valid blood vessel ROI',/error)
      return
    endelse 
  endif
  
  if ptr_valid(project.roi.pROIs[aifroi_caudate_set]) then aifroimask_caudate = (*project.roi.pROIs[aifroi_caudate_set])[aifroi_caudate_name] -> ComputeMask(Dimensions= [fdim,pdim]) $
  else begin
    msg = dialog_message('Select a valid caudate ROI',/error)
    return
  endelse
  
  ;Calculating the T10 Map that will be used to calculate the concentration at each slice
  if ~n_elements(flag) then T10_map_vessel = MAKE_ARRAY(fdim,pdim, /DOUBLE, VALUE = 2.7)
  p_T10arr = ptr_new(dblarr(fdim,pdim))
  dce_calculate_T10, p_T10arr, mask = aifroimask_caudate, slice = Ns_caudate
  T10_map_caudate = *p_T10arr
  neg_index = where(T10_map_caudate le 0)
  nan_index = where(~finite(T10_map_caudate))
  T10_map_caudate[neg_index] = 0
  T10_map_caudate[nan_index] = 0
  
  ;CALCULATING S(0) USING PRE IMAGES
  if ~n_elements(flag) then pre_avg_vessel = reform(mean((*project.dataArray[dce_num].state1)[*,*,Ns_vessel,0:n_Pre-1],dimension = 4))
  pre_avg_caudate = reform(mean((*project.dataArray[dce_num].state1)[*,*,Ns_caudate,0:n_Pre-1],dimension = 4))

  ;creating mask based on user-defined threshold
  if ~n_elements(flag) then begin
    max_pre = max(pre_avg_vessel)
    mask_threshold = mask_thresh*max_pre
    mask_vessel = pre_avg_vessel gt mask_threshold
    mask_vessel*=aifroimask    
  endif

  max_pre = max(pre_avg_caudate)
  mask_threshold = mask_thresh*max_pre
  mask_caudate = pre_avg_caudate gt mask_threshold
  mask_caudate*= aifroimask_caudate

  ;CALCULATING CONCENTRATION MAPS;;;;;;;;;;;;;;;;;;;;;;;;;
  ;create the 3D array that will hold the concentration maps
  if ~n_elements(flag) then Conc_3D_Array_Vessel = dblarr(fdim,pdim,n_Pst)
  Conc_3D_Array_Caudate = dblarr(fdim,pdim,n_Pst)

  progressBar = Obj_New('progressbar', Color='red', Text='Calculating Concentration ' ,/NOCANCEL)
  progressBar -> Start

  FOR time=n_Pre, n_Pre + n_Pst - 1 DO BEGIN
    time_dimension = time - n_Pre 
    if ~n_elements(flag) then begin
      dce_calculate_conc, output = conc_array_vessel, adim = time, T10_map = T10_map_Vessel, slice = Ns_vessel, mask = mask_vessel
      Conc_3D_Array_Vessel[*,*,time_dimension] = conc_array_vessel
    endif
    dce_calculate_conc, output = conc_array_caudate, adim = time, T10_map = T10_map_caudate, slice = Ns_caudate, mask = mask_caudate
    Conc_3D_Array_Caudate[*,*,time_dimension] = conc_array_caudate
    print, time_dimension
    progressBar -> Update,float(time)*100.0/float(n_Pst)
  ENDFOR

  progressBar -> Destroy

  IF smooth_ktrans_flag EQ 1 THEN BEGIN
    filter_size = 4
    if ~n_elements(flag) then begin
      Conc_3D_Array_Vessel_smooth  = dblarr(fdim,pdim,n_Pst)
      Conc_3D_Array_Vessel_smooth[*,*,0:filter_size-2] = Conc_3D_Array_Vessel[*,*,0:filter_size-2]
    endif
    Conc_3D_Array_Caudate_smooth  = dblarr(fdim,pdim,n_Pst)
    Conc_3D_Array_Caudate_smooth[*,*,0:filter_size-2] = Conc_3D_Array_Caudate[*,*,0:filter_size-2]
    ; Moving average filter
    FOR t = filter_size-1,n_Pst-1 DO begin
      if ~n_elements(flag) then Conc_3D_Array_Vessel_smooth[*,*,t] = mean(Conc_3D_Array_Vessel[*,*,t-filter_size+1:t],dimension = 3)  
      Conc_3D_Array_Caudate_smooth[*,*,t] = mean(Conc_3D_Array_Caudate[*,*,t-filter_size+1:t],dimension = 3)  ; Moving average filter
    ENDFOR
    if ~n_elements(flag) then Conc_3D_Array_Vessel = Conc_3D_Array_Vessel_smooth
    Conc_3D_Array_Caudate = Conc_3D_Array_Caudate_smooth
  ENDIF

  X = time_array()
  if n_elements(flag) gt 0 then begin
      AIFdata = AIF_data()
      Xact = AIFdata[0,*]
      Cpact = AIFdata[1,*]
      Cp = interpol(Cpact,Xact,X,/SPLINE)
  endif else  Cp = dblarr(n_Pst)
  C_cp = dblarr(n_Pst)
  for i=0,n_Pst-1 do begin
    if ~n_elements(flag) then begin
      temp1 = Conc_3D_Array_Vessel[*,*,i]
      Cp[i] = mean(temp1[where(temp1 ne 0,/NULL)],/double)
    endif
    temp2 = Conc_3D_Array_Caudate[*,*,i]    
    C_cp[i] = mean(temp2[where(temp2 ne 0,/NULL)],/double)
  endfor
  ind = where(X ge 5 and X le 10)
  eta = 100*int_tabulated(X[ind],C_cp[ind])/int_tabulated(X[ind],Cp[ind])
  Cp*=eta;

  A = [0.5,0.05,0.5,0.5]
  pinfo = replicate({value:0,limited:[0,0],limits:[0.0,0.0]},4)
  pinfo[*].value = A
  pinfo[*].limited[0] = [1,1,1,1]
  pinfo[*].limits[0]  =   [0,0,0,0]
  fargs = {x:X,y:Cp,err:1}
  params = mpfit('AIF_Function', functargs = fargs, parinfo=pinfo,bestnorm = chisquared,status=status)
  a1 = params[0]
  m1 = params[1]
  a2 = params[2]
  m2 = params[3]
  yfit = Dose*(params[0]*exp(-params[1]*X)+params[2]*exp(-params[3]*X))
  SS_res = chisquared
  SS_tot = (N_pst-1)*variance(Cp,/double)
  rsquared = 1- SS_res/SS_tot

  widget_control, kep_a1_text, set_value = strtrim(string(a1,format='(F8.3)'),2)
  widget_control, kep_a2_text, set_value = strtrim(string(a2,format='(F8.3)'),2)
  widget_control, kep_m1_text, set_value = strtrim(string(m1,format='(F8.3)'),2)
  widget_control, kep_m2_text, set_value = strtrim(string(m2,format='(F8.3)'),2)
  
  dimensions = GET_SCREEN_SIZE(RESOLUTION=resolution)
  resolution*=0.393701
  dims = long([3,4]/float(resolution))
  plot1 = plot(X,Cp,'b2+', DIMENSIONS = dims,NAME = "Actual", XTITLE = 'Time (mins)', YTITLE = '[Gd-DTPA]$_{blood}$ (mM)', FONT_SIZE = 11, FONT_NAME='Times', MARGIN = [0.18,0.12,0.05,0.05])
  plot2 = plot(X,yfit,'r2', DIMENSIONS = dims,/overplot, NAME = "Fitted",FONT_SIZE = 11, FONT_NAME='Times', MARGIN = [0.18,0.12,0.05,0.05])
  leg = LEGEND(TARGET=[plot1,plot2],FONT_SIZE = 11, POSITION = [270,350],/DEVICE)
  t1 = TEXT(55,270, $
    '$ C_p(t)=$'+strtrim(string(Dose,format='(F8.2)'),2)+'('+strtrim(string(a1,format='(F8.2)'),2)+'e!E-'+strtrim(string(m1,format='(F8.2)'),2)+'t!N+'+strtrim(string(a2,format='(F8.2)'),2)+'e!E-'+strtrim(string(m2,format='(F8.2)'),2)+'t!N)', $
    /DEVICE, FONT_SIZE=11, FONT_NAME='Times')
  t2 = TEXT(70,245, '$ r^2 = $'+strtrim(string(rsquared,format='(F8.2)'),2), FONT_SIZE=11, FONT_NAME='Times', /DEVICE)
  plot2.Save, "AIF.tiff", RESOLUTION=1200
  
  end
  
  pro dce_model_selection

  COMMON scan_data
  COMMON common_widgets
  COMMON dce_data

  fdim = project.ImndArray[dce_num].fdim*project.procPramArray[dce_num].Freq_interp
  pdim = project.ImndArray[dce_num].pdim*project.procPramArray[dce_num].Phase_interp
  n_Pre = project.imndArray[dce_num].n_Pre
  slice = project.procPramArray[dce_num].sdim_start
 
  ; Brain extraction mask
  if mask_num ne -1 then begin
    mask_bet = (*project.dataArray[mask_num].state1)[*,*,slice]
  endif else mask_bet = make_array(fdim,pdim,value = 1.0)
  ; ROI mask
  if ptr_valid(project.roi.pROIs[ktransroi_set]) then mask_ktrans = (*project.roi.pROIs[ktransroi_set])[ktransroi_name] -> ComputeMask(Dimensions= [fdim,pdim]) else $
    mask_ktrans = make_array(fdim,pdim,value = 1.0)

  ; Calculating T10 map
  p_T10arr = ptr_new(dblarr(fdim,pdim))
  dce_calculate_T10, p_T10arr
  T10_map = *p_T10arr
  
  ;create the 3D array that will hold the concentration maps
  Conc_3D_Array = dblarr(fdim,pdim,n_Pst)

  progressBar = Obj_New('progressbar', Color='red', Text='Calculating Concentration ',/NOCANCEL)
  progressBar -> Start

  FOR time=n_Pre,n_Pre + n_Pst-1 DO BEGIN
    dce_calculate_conc, output = conc_array, adim = time, mask = mask_ktrans*mask_bet, T10_map = T10_map
    time_dimension = time - n_Pre
    print, time_dimension+1
    Conc_3D_Array[*,*,time_dimension] = conc_array
    progressBar -> Update, (float(time_dimension)/float(n_Pst))*100.0
  ENDFOR

  progressBar -> Destroy
  
  IF smooth_ktrans_flag EQ 1 THEN BEGIN
    filter_size = 4
    Conc_3D_Array_smooth  = dblarr(fdim,pdim,n_Pst)
    Conc_3D_Array_smooth[*,*,0:filter_size-2] = Conc_3D_Array[*,*,0:filter_size-2]
    FOR t = filter_size-1,n_Pst-1 DO Conc_3D_Array_smooth[*,*,t] = mean(Conc_3D_Array[*,*,t-filter_size+1:t],dimension = 3)  ; Moving average filter
    Conc_3D_Array = Conc_3D_Array_smooth
  ENDIF

  kep = dblarr(fdim,pdim)
  ktrans = dblarr(fdim,pdim)
  vp = dblarr(fdim,pdim)
  model_select = dblarr(fdim,pdim)
  
  t = time_Array()
  Ct = dblarr(n_Pst)
  Cp = Dose*(a1*exp(-m1*t)+a2*exp(-m2*t))
 
  progressBar = Obj_New('progressbar', Color='red', Text='Fitting Model',/NOCANCEL)
  progressBar -> Start
  AIC = dblarr(4)
  i_pst = 0  
  for i = 0,fdim-1 do begin
    for j = 0,pdim-1 do begin      
      Ct = reform(Conc_3D_Array[i,j,*])          
      if array_equal(Ct,0) then continue   
      
      AIC[0] = n_Pst*alog(total(Ct^2)/n_Pst) 
        
      fargs   = {x:Cp,y:Ct,err:1}
      pinfo1    = replicate({value:1,limited:[1,1],limits:[0.0,1.0],fixed:0},2)
      pinfo1[1].fixed = 1
      pinfo1[1].value = 0
      params1  = MPFIT('linear_fit', functargs = fargs, parinfo=pinfo1,bestnorm = chisquared1,status=status,QUIET = 1)
      AIC[1] = n_Pst*alog(chisquared1/n_Pst) + 2 + 2*1*2/float(n_Pst-1-1) 
     
      pinfo2    = replicate({value:1,limited:[1,0],limits:[0.0,0.0]},2)
      pinfo2[1].limited[1] = 1
      pinfo2[1].limits[1] = 1.0
      params2  = MPFIT('tofts_model', functargs = fargs, parinfo=pinfo2,bestnorm = chisquared2,status=status,QUIET = 1)            
      AIC[2] = n_Pst*alog(chisquared2/n_Pst) + 4 + 2*2*3/float(n_Pst-2-1)
      
      pinfo3    = replicate({value:1,limited:[1,0],limits:[0.0,0.0]},3)
      pinfo3[2].limited[1] = 1
      pinfo3[2].limits[1] = 1.0
      params3  = MPFIT('std_model', functargs = fargs, parinfo=pinfo3,bestnorm = chisquared3,status=status,QUIET = 1)
      AIC[3] = n_Pst*alog(chisquared3/n_Pst) + 6 + 2*3*4/float(n_Pst-3-1)
       
      rAIC = AIC - min(AIC)
      ind = min(where(abs(rAIC) le 2))
      w = exp(-rAIC/float(2))/total(exp(-rAIC/float(2)))
;       AICmin = min(AIC,ind)
      case ind of
      0 : begin
          model_select[i,j] = 0
;          Ktrans[i,j] = 0
;          Kep[i,j] = 0
;          vp[i,j] = 0
          end
      1  : begin
           model_select[i,j] = 1
;           Ktrans[i,j] = 0
;           Kep[i,j] = 0
;           vp[i,j] = params1[0]
           end    
       2  : begin
            model_select[i,j] = 2
;            Ktrans[i,j] = params2[0]
;            Kep[i,j] = 0
;            vp[i,j] = params2[1]
            end  
       3  : begin
            model_select[i,j] = 3
;            Ktrans[i,j] = params3[0]
;            Kep[i,j] = params3[1]
;            vp[i,j] = params3[2]
            end                   
        endcase
        Ktrans[i,j] = params2[0]*w[2] + params3[0]*w[3]
        Kep[i,j] = params3[1]*w[3]
        vp[i,j] = params1[0]*w[1] + params2[1]*w[2] +  params3[2]*w[3]
        
      i_pst++
      progressBar -> Update, (float(i_pst)/float(fdim*pdim))*100.0
    endfor
  endfor
  progressBar -> Destroy
  
  if ptr_valid(project.dataarray[dce_num].dce_model_fit) then ptr_free,project.dataarray[dce_num].dce_model_fit
  project.dataarray[dce_num].dce_model_fit = ptr_new(make_array(fdim,pdim,4,value = 0,/double))
 
 (*project.dataarray[dce_num].dce_model_fit)[*,*,0] = model_select
 (*project.dataarray[dce_num].dce_model_fit)[*,*,1] = Ktrans
 (*project.dataarray[dce_num].dce_model_fit)[*,*,2] = Kep
 (*project.dataarray[dce_num].dce_model_fit)[*,*,3] = vp
 dce_model_fit_flag = 1
 
end

