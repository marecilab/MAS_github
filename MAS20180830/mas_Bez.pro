; MAS code to compute switched gradient eddy current magnetic field (Bez) using gradient echo

; Subroutine name: Bez_gui
; Created by: Magdoom Kulam
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Main eddy current GUI
;
; Editing Information:

pro Bez_gui
  common scan_data
  common common_widgets
  common Bez_params, T, axis, Xs, Ys, Zs, gBez, tBez, tasc, tplateau, tdesc, BW, beta, shape, scheme, eta, scale_sinc

  if ptr_valid(project.imndarray[project.ci].tgrdelay) then begin
    shape = project.imndarray[project.ci].tgrshape
    case shape of
      't' : get_Bez_params, T, axis, shape, gBez, tBez, BW, tasc = tasc, tplateau = tplateau, tdesc = tdesc, scheme
      'a' : get_Bez_params, T, axis, shape, gBez, tBez, BW, tasc = tasc, tplateau = tplateau, tdesc = tdesc, scheme
      'f' : get_Bez_params, T, axis, shape, gBez, tBez, BW, beta = beta, scheme
      'i' : get_Bez_params, T, axis, shape, gBez, tBez, BW, beta = beta, scheme, eta = eta, scale_sinc = scale_sinc
      else : get_Bez_params, T, axis, shape, gBez, tBez, BW, scheme
    endcase
  endif else begin
    mbox = dialog_message('Select a eddy current measurement dataset',/error)
    return
  endelse
  Nz   = project.imndarray[project.ci].sdim                     ; Number of slices

  if (xregistered('Bez_gui') ne 0) then return
  topbase = widget_base(title = 'Gradient Eddy Current Field Mapping : Scan #' + strtrim(project.ci+1,2) ,tlb_frame_attr=1,/grid_layout,/kbrd_focus_events)
  mainbase = widget_base(topbase,column=3)

  labelbp      = widget_label(mainbase, value = 'Acquisition parameters')
  fbase        = widget_base(mainbase)
  parbase_Bez  = widget_base(fbase, row = 7, frame = 1)
  gaxis        = cw_field(parbase_Bez, title = 'Axis                  ', value = axis ,xsize = 14, /noedit)
  schfield_Bez = cw_field(parbase_Bez, title = 'Encoding Scheme       ', value = scheme,xsize = 14,/noedit)
  gfield_Bez   = cw_field(parbase_Bez, title = 'Peak Gradient Strength', value = strtrim(string(gBez[0],format = '(F8.2)'),2)+' mT/m',xsize = 14, /noedit)
  tfield_Bez   = cw_field(parbase_Bez, title = 'Gradient Duration     ', value = strtrim(string(tBez,format = '(F8.2)'),2)+ ' ms',xsize = 14,/noedit)
  shpfield_Bez = cw_field(parbase_Bez, title = 'Gradient Shape        ', value = shape,xsize = 14,/noedit)
  BW_Bez       = cw_field(parbase_Bez, title = 'Excitation Bandwidth  ', value = strtrim(string(BW,format = '(F8.2)'),2)+ ' KHz',xsize = 14,/noedit)
  Tbase_Bez  = widget_base(parbase_Bez, column = 2)
  Ttext_Bez = widget_label(Tbase_Bez, value =  'Field Sampling Time    ')
  Tlist_Bez = widget_list(Tbase_Bez, value = strtrim(string(T,format = '(F8.3)') + ' ms',2), uvalue = 'Time List',Ysize = 3, Xsize = 10)

  visbase_Bez =  widget_base(mainbase, row = 4)
  vis_button = ['Field Map','Time Series','Spectra','Tensor','Model Fit']
  vis_button_grp_Bez = cw_bgroup(visbase_Bez,vis_button, uvalue = 'vis_button', column = 1, exclusive = 1, button_uvalue = vis_button, label_top = 'Response', frame = 2,/no_release , set_value=0)
  model_button = ['Eddy Current','Oscillations','Full model']
  model_button_grp_Bez = cw_bgroup(visbase_Bez,model_button, uvalue = 'model_button', column = 1, exclusive = 1, button_uvalue = model_button, label_top = 'Model type', frame = 2,/no_release , set_value=0)
  axes = ['X','Y','Z']
  Axlist_Bez = widget_droplist(visbase_Bez, title = 'Axis', value = axes, uvalue = 'Axes',/DYNAMIC_RESIZE)

  Bezbase   =  widget_base(mainbase, row = 6)
  uwbtnbase = widget_base(Bezbase, /nonexclusive)
  pc_button_Bez = widget_button(uwbtnbase, value = 'Unwrap phase', uvalue = 'UnwrapPhase')
  sm_button_Bez = widget_button(uwbtnbase, value = 'Smooth Time Series', uvalue = 'Smoothing')
  ap_button_Bez = widget_button(uwbtnbase, value = 'Apodize', uvalue = 'Apodize')
  zf_button_Bez = widget_button(uwbtnbase, value = 'Zero fill (x2)', uvalue = 'Zero fill')
  Bez_constraint = cw_fslider(Bezbase, title = 'Amplitude Constraint (%)', uvalue = 'Amplitude Constraint',xsize = 165, /edit)
  Nfit_Bez = widget_slider(Bezbase, title = 'Points to fit', uvalue = 'FitPoints',xsize = 165, min = 4, max = n_elements(T))
  Ithres_Bez = widget_slider(Bezbase, title = 'Image Threshold (%)', uvalue = 'Image Threshold',xsize = 165)

  disp_button_base = widget_base(Bezbase, column = 2)
  movie_button_Bez = widget_button(disp_button_base, value = 'Movie', uvalue = 'Movie',/no_release, xsize = 82)
  disp_button_Bez = widget_button(disp_button_base, value = 'Display', uvalue = 'Display',/no_release, xsize = 82)
  export_button_Bez = widget_button(disp_button_base, value = 'Export', uvalue = 'Export',/no_release, xsize = 82)

  if project.procPramArray[project.ci].single_Multi_flag eq 1 then $
    widget_control, movie_button_Bez, sensitive = 1 else widget_control, movie_button_Bez, sensitive = 0
  if scheme ne 'Hadamard' and scheme ne 'Four Point' and scheme ne 'Six Point' then widget_control,Axlist_Bez, sensitive  = 0 else $
    widget_control,Axlist_Bez, sensitive  = 1
  widget_control,model_button_grp_Bez , sensitive = 0
  widget_control,zf_button_Bez, sensitive = 0
  widget_control,ap_button_Bez, sensitive = 0

  widget_control, topbase, /realize
  xmanager, 'Bez_gui', topbase, cleanup='Bez_gui_cleanup',/NO_BLOCK

end

; Subroutine name: Bez_gui_event
; Created by: Magdoom Kulam
; Calling Information:
;
;  event
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Event handler for the main eddy current gui
;
; Editing Information:

pro Bez_gui_event, event

  common scan_data
  common common_widgets
  common Bez_params

  ; To catch error ahead and prevent crashing
  catch,error_status
  IF (error_status NE 0) THEN BEGIN
    help, calls= trace_back
    dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
      !ERROR_STATE.msg,trace_back], /ERROR, TITLE='Error in WID_BASE_MAIN_event')
    widget_control, event.top, /destroy
    RETURN
  ENDIF

  widget_control, event.id, get_uvalue = widget
  CI = project.ci

  if (event.top eq event.id) then begin
    widget_control, event.top, tlb_set_title = 'Gradient Eddy Current Field Mapping : Scan #' + strtrim(CI+1,2)
    if ptr_valid(project.imndarray[CI].tgrdelay) then begin
      shape = project.imndarray[project.ci].tgrshape
      case shape of
        't' : get_Bez_params, T, axis, shape, gBez, tBez, BW, tasc = tasc, tplateau = tplateau, tdesc = tdesc, scheme
        'a' : get_Bez_params, T, axis, shape, gBez, tBez, BW, tasc = tasc, tplateau = tplateau, tdesc = tdesc, scheme
        'f' : get_Bez_params, T, axis, shape, gBez, tBez, BW, beta = beta, scheme
        'i' : get_Bez_params, T, axis, shape, gBez, tBez, BW, beta = beta, scheme, eta = eta, scale_sinc = scale_sinc
        else : get_Bez_params, T, axis, shape, gBez, tBez, BW, scheme
      endcase
      widget_control,Bezbase, sensitive = 1
      widget_control,parbase_Bez, sensitive = 1
      widget_control,visbase_Bez, sensitive = 1

      widget_control,gaxis, set_value = axis
      widget_control,schfield_Bez, set_value = scheme
      widget_control,gfield_Bez, set_value = strtrim(string(gBez[project.procPramArray[CI].Bez_axis],format = '(F8.2)'),2)+' mT/m'
      widget_control,tfield_Bez, set_value = strtrim(string(tBez,format = '(F8.2)'),2)+' ms'
      widget_control,shpfield_Bez, set_value = shape
      widget_control,BW_Bez, set_value = strtrim(string(BW,format = '(F8.2)'),2)+' KHz'

      time_list = strtrim(string(T,format='(F8.3)') +' ms',2)
      widget_control,Tlist_Bez, set_value = time_list, set_list_select = project.procPramArray[CI].Bez_time_index
      widget_control,Axlist_Bez, set_droplist_select = project.procPramArray[CI].Bez_axis
      widget_control,vis_button_grp_Bez, get_value = vis_type
      if scheme ne 'Hadamard' and scheme ne 'Four Point' and scheme ne 'Six Point' then widget_control,Axlist_Bez, sensitive  = 0 else widget_control,Axlist_Bez, sensitive  = 1

      widget_control,pc_button_Bez, set_button = project.procPramArray[CI].Bez_up_stat
      widget_control,ap_button_Bez, set_button = project.procPramArray[CI].Bez_apodize
      widget_control,zf_button_Bez, set_button = project.procPramArray[CI].Bez_zf_stat
      widget_control,sm_button_Bez, set_button = project.procPramArray[CI].Bez_smooth
      widget_control,Bez_constraint, set_value = project.procPramArray[CI].Bez_alpha_max
      widget_control,Ithres_Bez, set_value = project.procPramArray[CI].Bez_Ithreshold
      widget_control,Nfit_Bez, set_slider_max = n_elements(T)
      widget_control,Nfit_Bez, set_value = project.procPramArray[CI].Bez_Nfit
      sort_dir = project.procPramArray[CI].sort_dir
      if project.procPramArray[CI].single_Multi_flag eq 1 then begin
        if sort_dir eq 1 then widget_control,Tlist_Bez, sensitive = 0 else widget_control,Tlist_Bez, sensitive = 1
        if vis_type eq 0 then widget_control, movie_button_Bez, sensitive = 1 else widget_control, movie_button_Bez, sensitive = 0
      endif else widget_control, movie_button_Bez, sensitive = 0

    endif else begin
      widget_control,Bezbase, sensitive = 0
      widget_control,visbase_Bez, sensitive = 0
      widget_control,parbase_Bez, sensitive = 0
    endelse

  endif else begin

    case widget of

      'UnwrapPhase'          : begin
        project.procPramArray[CI].Bez_up_stat = event.select
        project.procPramArray[CI].Bez_proccess_flag = 0
      end

      'Amplitude Constraint' : project.procPramArray[CI].Bez_alpha_max = event.value

      'FitPoints'            : project.procPramArray[CI].Bez_Nfit = event.value

      'Image Threshold'      : project.procPramArray[CI].Bez_Ithreshold = event.value

      'Apodize'              : project.procPramArray[CI].Bez_apodize = event.select

      'Smoothing'            : project.procPramArray[CI].Bez_smooth = event.select

      'Zero fill'            : project.procPramArray[CI].Bez_zf_stat = event.select

      'Time List'            : project.procPramArray[CI].Bez_time_index = event.index

      'Axes'                 : project.procPramArray[CI].Bez_axis = event.index

      'Export'               : export_Bez

      'vis_button'           : begin
        widget_control, event.id, get_value = vis_type

        if vis_type eq 0 or vis_type eq 3 then  begin
          widget_control,Tlist_Bez, sensitive = 1
          widget_control, movie_button_Bez, sensitive = 1
          widget_control,sm_button_Bez, sensitive = 0
          widget_control,zf_button_Bez, sensitive = 0
          widget_control,ap_button_Bez, sensitive = 0
        endif else begin
          widget_control,Tlist_Bez, sensitive = 0
          widget_control, movie_button_Bez, sensitive = 0
          widget_control,sm_button_Bez, sensitive = 1
          widget_control,zf_button_Bez, sensitive = 1
          widget_control,ap_button_Bez, sensitive = 1
        endelse

        if vis_type eq 4 then widget_control,model_button_grp_Bez, sensitive = 1 $
        else widget_control,model_button_grp_Bez, sensitive = 0

        if vis_type eq 3 then widget_control,Axlist_Bez, sensitive = 0 $
        else widget_control,Axlist_Bez, sensitive = 1
      end

      'Display'             : begin
        if (project.procPramArray[CI].Bez_proccess_flag eq 0 or ~ptr_valid(project.dataarray[CI].Bez) or project.procPramArray[CI].state_1 eq 0) then Calculate_Bez, project.ci
        widget_control, vis_button_grp_Bez, get_value = vis_type
        widget_control, model_button_grp_Bez, get_value = model_time
        case vis_type of
          0 : if project.procPramArray[CI].single_Multi_flag eq 0 then Bez_single_display else Bez_multi_display
          1 : TimeSeries
          2 : Spectrum
          3 : Tensor_Bez
          4 : Plot_Bez_fit, model_time
        endcase
      end

      'Movie'              : begin
        if (project.procPramArray[CI].Bez_proccess_flag eq 0 or ~ptr_valid(project.dataarray[CI].Bez) or project.procPramArray[CI].state_1 eq 0) then Calculate_Bez,project.ci
        Bez_movie
      end

      else             : return
    endcase
  endelse
end

;+
; NAME:
;   PHUNWRAP
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;   UPDATED VERSIONs can be found on my WEB PAGE:
;      http://cow.physics.wisc.edu/~craigm/idl/idl.html
;
; PURPOSE:
;   Unwrap phase jumps to recover cycle counts
;
; MAJOR TOPICS:
;   Mathematics
;
; CALLING SEQUENCE:
;   CYCLES = PHUNWRAP(PHASE, TOLERANCE=, MAXVAL=)
;
; DESCRIPTION:
;
;   PHUNWRAP unwraps a sequence of phases to produce a new series of
;   cycle counts.  Phase jumps due to crossing over the PHASE=0
;   boundary are removed by adding an integral number of cycles.  The
;   algorithm is based on the MATLAB "unwrap" function.
;
;   NOTE: the unwrapping process can be ambiguous if there is a phase
;   jump of more than a half cycle in the series.  For example, if the
;   phase changes by ~0.5 cycles, it is not possible to distinguish
;   whether there wasa +0.5 cycle or -0.5 cycle jump.  The most
;   accurate unwrapping can be performed if the PHASE series is nearly
;   continuous and does not have rapid phase changes.
;
;   Users can select the tolerance used to determine the phase jump.
;   Users can also select the definition of "1 cycle" by changing
;   MAXVAL.  By default, MAXVAL is 2*!DPI, which correspondes to 1
;   cycle = 2*!DPI radians, but other values of 1 (cycle), or 360
;   (degrees) are possible.
;
; INPUTS:
;
;   PHASE - phase series to be unwrapped.  Values should range from 0
;           to MAXVAL.  The ordering of the series is important.
;
; RETURNS:
;
;   A new series, expressed in cycles, with cycle jumps larger than
;   TOLERANCE removed.
;
; OPTIONAL KEYWORDS:
;
;   TOLERANCE - phase jump tolerance.  If the phase from one sample to
;               the next changes by more than TOLERANCE, then a single
;               cycle jump is assumed to have occurred.
;               DEFAULT: 0.5*MAXVAL
;
;   MAXVAL - Maximum value for phase. Common values are: 2*!DPI
;            (radians; DEFAULT); 1 (cycle); 360 (degrees), but any
;            positive value may be used.
;
; EXAMPLE:
;
;  ;; Set up some fake data
;  x = dindgen(100)/10d
;  y = x/2
;  ph = y MOD 1.0    ;; Mock phases
;
;  cycles = phunwrap(ph, maxval=1)
;
; MODIFICATION HISTORY:
;   Written and documented, CM, July 2003
;   Handle the case of unsigned integer input, CM, Feb 2006
;
;  $Id: phunwrap.pro,v 1.3 2006/03/28 14:19:53 craigm Exp $
;
;-
; Copyright (C) 2003, 2006, Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

function phunwrap, ph, tolerance=tol0, maxval=maxval0

  common phunwrap_common, idlver
  if n_elements(idlver) EQ 0 then begin
    idlver = !version.release
  endif


  if n_elements(maxval0) EQ 0 then maxval = 2d*!dpi else maxval = maxval0(0)
  if n_elements(tol0) EQ 0 then tol = 0.5*maxval else tol = tol0(0)*maxval

  if n_elements(ph) LT 2 then return, ph

  sz = size(ph)
  tp = sz(sz(0)+1)

  ;; First order difference
  case tp of
    12: dph = [0, long(ph)-long(ph(1:*))]
    13: dph = [0, long64(ph)-long64(ph(1:*))]
    15: dph = [0, long64(ph)-long64(ph(1:*))]
    else: dph = [0, ph - ph(1:*)]
  endcase

  p = maxval * (fix((dph GT tol) EQ 1) - fix((dph LT (-tol)) EQ 1))
  if idlver GT 5.25 then begin
    ;; Use built-in version if available
    r = total(p, /cumulative)
  endif else begin
    ;; .. if not, then use the lame FOR loop
    r = p
    for i = 1L, n_elements(r)-1 do $
      r(i) = r(i) + r(i-1)
  endelse

  return, ph+r
end

; Subroutine name: gradient_specs
; Created by: Magdoom Kulam
; Calling Information:
;
; GradCoil  - Gradient coil name
; Pmax      - Maximum power dissipation (kW)
; Gmax      - Maximum gradient strength (G/cm)
; Irms_max  - Maximum RMS current (A)
; DeltaT    - Maximum change in coil temperature (degree celsius)
; alpha     - Temperature coefficient of resistance for the coil (1/Celsius)
; C         - Gradient strength coeff (G/cm/A)
; C0        - Strength of B0 shim coil (uT/A)
; Imax_B0   - Maximum current on B0 shim coil
; R0        - Resistance of gradient coil (millOhms)
; R         - Resistance adjusted for change in temperature (millOhms)
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Specification for gradient coil
;
; Editing Information:

function gradient_specs, GradCoil

  gammabar = 42.58      ;Proton gyromagnetic ratio MHz/T
  
  GradCoilParam = {Pmax:0.0, Gmax:0.0, Irms_max:0.0, DTemp:0.0, C:fltarr(3), C0:0.0, Imax_B0:0.0, R:fltarr(3)}
  case GradCoil of
    'RRI_240_120_S7' : begin
      GradCoilParam.Pmax = 1.5
      GradCoilParam.Gmax = 100
      GradCoilParam.Irms_max = 75
      GradCoilParam.DTemp = 45
      GradCoilParam.C = [0.383,0.43,0.313] 
      GradCoilParam.C0 = 6000/gammabar
      GradCoilParam.Imax_B0 = 10
      alpha = 0.00393
      R0 = [150.3,166.5,121.6]
    end

    'RRI_400_260_S7' : begin
      GradCoilParam.Pmax = 1.5
      GradCoilParam.Gmax = 15.5
      GradCoilParam.Irms_max = 75
      GradCoilParam.DTemp = 45
      GradCoilParam.C = [0.0785,0.0782,0.0778]
      GradCoilParam.C0 = 3812/gammabar
      GradCoilParam.Imax_B0 = 2
      alpha = 0.00393
      R0 = [201,185,316]
    end

    'RRI_200_115_S14' : begin
      GradCoilParam.Pmax = 1
      GradCoilParam.Gmax = 67
      GradCoilParam.Irms_max = 75
      GradCoilParam.DTemp = 45
      GradCoilParam.C = [0.231,0.233,0.230]
      GradCoilParam.C0 = 4232/gammabar
      GradCoilParam.Imax_B0 = 2
      alpha = 0.00393
      R0 = [98,121,122]
    end
  endcase
  GradCoilParam.R = R0*(1+alpha*GradCoilParam.DTemp)

  return, GradCoilParam
end


; Subroutine name: FresnelS
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Evaluate the Fresnel sine integral for a given 'x'
;
; Editing Information:

function FresnelS,x

  erf1 = erfz(dcomplex(1,1)*x*sqrt(!pi)/2.0)
  erf2 = erfz(dcomplex(1,-1)*x*sqrt(!pi)/2.0)
  y = dcomplex(1,1)*dcomplex(real_part(erf1) + imaginary(erf2), imaginary(erf1)-real_part(erf2))/4.0
  return, y

end

; Subroutine name: Si
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Evaluate the Sine integral function for a given 'x'
;
; Editing Information:

function Si,x
  
  if x eq 0 then return, 0.0 else begin
    y = x*findgen(100)/99
    f = sin(y)/y
    f[where(finite(f) eq 0)] = 1

    return, int_tabulated(y,f)
    
  endelse

end

; Subroutine name: Ei
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Evaluate the exponential integral function for any complex number, z
;
; Editing Information:


function Ei,z

  gamma = double(0.57721)                 ; Euler–Mascheroni constant 
  if ~isa(z,/complex) then z = dcomplex(z,0.0)
  
  if abs(real_part(z)) le 40 and abs(imaginary(z)) le 35 then begin
    ; Taylor expansion for small arguments
    E0 = -gamma - alog(-z)
    E1 = -gamma - alog(-z) - z
    i = 2
    while abs(real_part(E0-E1)) gt 1e-3 or abs(imaginary(E0-E1)) gt 1e-3 do begin
      E0 = E1
      E1 = E0 - z^i/double(i*factorial(i))
      i++
    endwhile
  endif else begin
    ; Asymptotic expansion for large arguments
    E0 = -exp(z)/z
    E1 = E0 - exp(z)/z^2
    i = 2
    while abs(real_part(E0-E1)) gt 1e-3 or abs(imaginary(E0-E1)) gt 1e-3 do begin
      E0 = E1
      E1 = E0 + (-1)^i*exp(z)*factorial(i)/(-z)^(i+1)
      i++
    endwhile      
  endelse
  
  ; Branch cut
  if imaginary(z) lt 0 then f = -E1 - dcomplex(0,!Pi) $
    else if imaginary(z) ge 0  then f = -E1 + dcomplex(0,!Pi) $
    else f = -E1 
  return, f

end

; Subroutine name: GammaU
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Evaluate the incomplete gamma function, Gamma(0) for any complex, z as defined in Mathematica
;
; Editing Information:

function GammaU,z

  return, -Ei(-z) + 0.5*(alog(-z)-alog(-1/z))-alog(z)

end

; ERFZ  Error function for complex inputs
;   Z may be complex and of any size.
;   Accuracy is better than 12 significant digits.
;
;
;   Ref: Abramowitz & Stegun section 7.1
;   equations 7.1.9, 7.1.23, and 7.1.29
;
;   Tested under version 5.3.1
;
;   See also erf, erfc, erfcx, erfinc, erfcore

;   Main author Paul Godfrey <pgodfrey@intersil.com>
;   Small changes by Peter J. Acklam <jacklam@math.uio.no>
;   09-26-01
;   Modified by Magdoom for IDL to evaluate for a scalar complex argument 
;   09-16-17

function erfz,zz

i = dcomplex(0,1)
twopi = 2*!Pi;
sqrtpi=1.772453850905516027298;

f = 0
ff=f;

az = abs(zz);
p1 = 0
p2 =0

if az le 8 then begin
  z=zz
  nn = 32;
  x = real_part(z);
  y = imaginary(z);
  k1 = 2/!Pi*exp(-x*x);
  k2 = exp(-i*2*x*y);

  s1 = erf(x);

  s2 = 0
  if x ne 0 then s2 = k1/ (4.0*x)*(1 - k2) else s2 = i / !pi * y
  f = s1 + s2;
  if y ne 0 then begin
    xk = x
    yk = y;
  endif else begin
    xk = 0
    yk = 0
  endelse
  s5 = 0;
  for n = 1, nn do begin
    s3 = exp(-n*n/4.0)/(n*n + 4.0*xk*xk);
    s4 = 2*xk - k2*(2*xk*cosh(n*yk) - i*n*sinh(n*yk));
    s5 = s5 + s3*s4;
  endfor

  s6 = k1*s5;
  f = f + s6;
endif else begin  
  z=zz;
  if real_part(z) lt 0 then z = -z

  nmax=193;
  s=1;
  y=2*z*z;
  for n=nmax,1,-2 do s=1-n*s/y

  f=1.0-s*exp(-z*z)/(sqrtpi*z);
  if real_part(z) lt 0 then f = -f
  if real_part(z) eq 0 then f = f-1
  
endelse

return,f
end


; Subroutine name: erfi
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Evaluate the imaginary error function for a given 'x'
;
; Editing Information:
function erfi,x

  i = dcomplex(0,1)
  y = erfz(i*x)/i
  return, y

end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_X,p,X=x,Y=y,ERR=err

  model = p[0]*X[*,0]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_Y,p,X=x,Y=y,ERR=err

  model = p[0]*X[*,1]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_Z,p,X=x,Y=y,ERR=err

  model = p[0]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_C,p,X=x,Y=y,ERR=err

  model = p[0]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_XY,p,X=x,Y=y,ERR=err

  model = p[0]*X[*,0] + p[1]*X[*,1]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_XZ,p,X=x,Y=y,ERR=err

  model = p[0]*X[*,0] + p[1]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_YZ,p,X=x,Y=y,ERR=err

  model = p[0]*X[*,1] + p[1]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_CX,p,X=x,Y=y,ERR=err

  model = p[0] + p[1]*X[*,0]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_CY,p,X=x,Y=y,ERR=err

  model = p[0] + p[1]*X[*,1]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_CZ,p,X=x,Y=y,ERR=err

  model = p[0] + p[1]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_CXY,p,X=x,Y=y,ERR=err

  model = p[0] + p[1]*X[*,0] + p[2]*X[*,1]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_CXZ,p,X=x,Y=y,ERR=err

  model = p[0] + p[1]*X[*,0] + p[2]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_CYZ,p,X=x,Y=y,ERR=err

  model = p[0] + p[1]*X[*,1] + p[2]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_XYZ,p,X=x,Y=y,ERR=err

  model = p[0]*X[*,0] + p[1]*X[*,1] + p[2]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_space_mpfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Procedure used in Craig B. Markwardt's MPFIT function to fit the field to linear function
;
; Editing Information:

function Bez_fit_space_CXYZ,p,X=x,Y=y,ERR=err

  model = p[0] + p[1]*X[*,0] + p[2]* X[*,1] + p[3]*X[*,2]
  return, (y-model)/err
end

; Subroutine name: Bez_fit_time_trapz
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Time dependent eddy current distortion model for trapezoidal test pulse
;
; Editing Information:
function Bez_fit_time_trapz_eddy,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  n = n_elements(p)
  model = 0
  for i=0,n-1,2 do model += p[i]*p[i+1]*gBez[ax]*((1-exp(-tdesc/p[i+1]))/tdesc-(1-exp(-tasc/p[i+1]))*exp(-(tdesc+tplateau)/p[i+1])/tasc)*exp(-X/p[i+1])
  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_trapz_fullmodel
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model for eddy current distortion and oscillation following a trapezoidal test pulse
;
; Editing Information:
function Bez_fit_time_trapz_oscillation,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0
  for j=0,n-1,4 do begin
    alpha = p[j]
    tau = p[j+1]
    f = p[j+2]
    phi0 = p[j+3]
    model += imaginary(tasc^(-1.0)*tdesc^(-1.0)*exp(1.0)^((I*2.0)*f*!Pi*X+(-1.0)*X*tau^(-1.0)+I*phi0)*(tasc*((-1.0)+exp(1.0)^((I*2.0)*tdesc*f*!Pi+(-1.0)*tdesc*tau^(-1.0)))+tdesc*(exp(1.0)^(I*(tdesc+tplateau)*tau^(-1.0)*(I+2.0*f*!Pi*tau))+(-1.0)*exp(1.0)^(I*(tasc+tdesc+tplateau)*tau^(-1.0)*(I+2.0*f*!Pi*tau))))*alpha*tau*((-1.0)+(I*2.0)*f*!Pi*tau)^(-1.0)*gBez[ax])
  endfor
  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_sine
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Time dependent eddy current distortion model for half-sinusoidal test pulse
;
; Editing Information:
function Bez_fit_time_sine_eddy,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  n = n_elements(p)
  model = 0
  for i=0,n-1,2 do model += !Pi*p[i]*p[i+1]*gBez[ax]*tBez*(1+exp(-tBez/p[i+1]))*exp(-X/p[i+1])/(tBez^2+!Pi^2*p[i+1]^2)
  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_sine_fullmodel
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model for eddy current distortion and oscillation following a half-sinusoidal test pulse
;
; Editing Information:
function Bez_fit_time_sine_oscillation,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0
  for j=0,n-1,4 do begin
    alpha = p[j]
    tau = p[j+1]
    f = p[j+2]
    phi0 = p[j+3]
    model += imaginary(I*exp(1.0)^((-1.0)*tau^(-1.0)*(X+tBez+(I*(-2.0))*f*!Pi*X*tau+(I*(-1.0))*tau*phi0))*(exp(1.0)^((I*2.0)*f*!Pi*tBez)+exp(1.0)^(tBez*tau^(-1.0)))*!Pi*tBez*alpha*tau*(I+2.0*f*!Pi*tau)*(I*tBez+!Pi*((-1.0)+2.0*f*tBez)*tau)^(-1.0)*(!Pi*tau+tBez*(I+2.0*f*!Pi*tau))^(-1.0)*gBez[ax])
  endfor
  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_scosine
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Time dependent eddy current distortion model for half-shifted cosine test pulse
;
; Editing Information:
function Bez_fit_time_scosine_eddy,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  n = n_elements(p)
  model = 0
  for i=0,n-1,2 do model += 2*!Pi^2*p[i]*p[i+1]^2*gBez[ax]*(1-exp(-tBez/p[i+1]))*exp(-X/p[i+1])/(tBez^2+4*!Pi^2*p[i+1]^2)
  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_scosine_fullmodel
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model for eddy current distortion and oscillation following a half-shifted cosine test pulse
;
; Editing Information:
function Bez_fit_time_scosine_oscillation,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0
  for j=0,n-1,4 do begin
    alpha = p[j]
    tau = p[j+1]
    f = p[j+2]
    phi0 = p[j+3]
    model += imaginary((-2.0)*exp(1.0)^((-1.0)*tau^(-1.0)*(t+tBez+(I*(-2.0))*f*!Pi*t*tau+(I*(-1.0))*tau*phi0))*((-1.0)*exp(1.0)^((I*2.0)*f*!Pi*tBez)+exp(1.0)^(tBez*tau^(-1.0)))*!Pi^2.0*alpha*tau^2.0*(I*tBez+2.0*!Pi*((-1.0)+f*tBez)*tau)^(-1.0)*(I*tBez+2.0*!Pi*(1.0+f*tBez)*tau)^(-1.0)*gBez[ax])
  endfor
  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_Gauss
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Time dependent eddy current distortion model for Gaussian test pulse (width of the Gaussian function set as on VnmrJ)
;
; Editing Information:
function Bez_fit_time_Gauss_eddy,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  sigma = tBez/(9.6*sqrt(2*!Pi))

  n = n_elements(p)
  model = 0
  for i=0,n-1,2 do model += p[i]*gBez[ax]*exp(-tBez^2/(8*sigma^2))*(-2*(1-exp(-tBez/p[i+1]))+$
    exp(-tBez/p[i+1])*exp((2*sigma^2+tBez*p[i+1])^2/(8*sigma^2*p[i+1]^2))*sqrt(2*!Pi)*sigma*$
    (erfz((2*sigma^2+tBez*p[i+1])/(2*sqrt(2)*sigma*p[i+1]))- $
    erfz((2*sigma^2-tBez*p[i+1])/(2*sqrt(2)*sigma*p[i+1]))))*exp(-X/p[i+1])/2.0

  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_Gauss_fullmodel
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model for eddy current distortion and oscillation following a half-shifted cosine test pulse
;
; Editing Information:
function Bez_fit_time_Gauss_oscillation,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  sigma = tBez/(9.6*sqrt(2*!Pi))
  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0
  for j=0,n-1,4 do begin
    alpha = p[j]
    tau = p[j+1]
    f = p[j+2]
    phi0 = p[j+3]
    model += imaginary(2.0^(-1.0/2.0)*exp(1.0)^((1.0/4.0)*((I*8.0)*f*!Pi*X+(-1.0)*tBez^2.0*sigma^(-2.0)+(-8.0)*f^2.0*!Pi^2.0*sigma^2.0+(-2.0)*tau^(-2.0)*(sigma^2.0+(2.0*X+tBez)*tau)+(I*4.0)*phi0))*alpha*tau^(-1.0)*(2.0^(1.0/2.0)*exp(1.0)^(2.0*f^2.0*!Pi^2.0*sigma^2.0+(1.0/8.0)*sigma^(-2.0)*tau^(-2.0)*((-2.0)*sigma^2.0+tBez*tau)^2.0)*(exp(1.0)^((I*2.0)*f*!Pi*tBez)+(-1.0)*exp(1.0)^(tBez*tau^(-1.0)))*tau+exp(1.0)^((1.0/4.0)*tBez^2.0*sigma^(-2.0)+sigma^2.0*tau^(-2.0)+I*f*!Pi*tau^(-1.0)*((-2.0)*sigma^2.0+tBez*tau))*!Pi^(1.0/2.0)*sigma*(1.0+(I*(-2.0))*f*!Pi*tau)*(erfz((1.0/2.0)*2.0^(-1.0/2.0)*(tBez*sigma^(-1.0)+2.0*sigma*((I*(-2.0))*f*!Pi+tau^(-1.0))))+erfz((1.0/2.0)*2.0^(-1.0/2.0)*sigma^(-1.0)*tau^(-1.0)*(tBez*tau+sigma^2.0*((-2.0)+(I*4.0)*f*!Pi*tau)))))*gBez[ax])
  endfor
  return, (y-model)/err

end


; Subroutine name: Bez_fit_time_FresnelT
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Time dependent eddy current distortion model for Fresnel test pulse with Tukey window
;
; Editing Information:
; Model expression commented does not account for windowing
function Bez_fit_time_FresnelT_eddy,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0

  for j=0,n-1,2 do begin
    alpha = p[j]
    tau = p[j+1]

    model+=(-((alpha*(exp((I*tBez)/(4.0*!Pi*BW*tau^2.0))*tBez*$
      (exp((4.0*tBez^2.0*BW*(beta + 2.0*I*!Pi*BW*tau + 1)*beta^2.0 + 3.0*I*!Pi*tau +$
      2.0*tBez*(6*I*!Pi*BW*tau*beta + beta))/(4.0*tBez*beta^2.0*BW*tau))*(4.0 - 4.0*I)*tBez^2.0*$
      Erfi(((1/2.0 + I/2.0)*Sqrt(tBez/BW))/(Sqrt(2.0*!Pi)*tau))*beta^2.0 -$
      I^(2.0*1/4.0)*Sqrt(2.0)*exp((I*!Pi*(2.0*tBez*beta*BW + 1)^2.0)/(2.0*tBez*beta^2.0*BW) + (tBez*(beta + 1))/tau +$
      (I*tBez)/(2.0*!Pi*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 - I^(2.0*3.0/4.0)*Sqrt(2.0)*$
      exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 4.0/beta + 1/(tBez*beta^2.0*BW)))*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 - 2.0*I^(2.0*1/4.0)*Sqrt(2.0)*$
      exp((2.0*tBez^2.0*(4.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 2.0*!Pi*(beta + 1)*BW*tau + I)*beta^2.0 +$
      2.0*!Pi*tBez*tau*(6*I*!Pi*BW*tau + 1)*beta + 3.0*I*!Pi^2.0*tau^2.0)/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(-I + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      2.0*I^(2.0*3.0/4.0)*Sqrt(2.0)*exp((4.0*tBez^2.0*BW*(beta + 2.0*I*!Pi*BW*tau + 1)*beta^2.0 + 3.0*I*!Pi*tau +$
      2.0*tBez*(6*I*!Pi*BW*tau*beta + beta))/(4.0*tBez*beta^2.0*BW*tau))*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      I^(2.0*1/4.0)*exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) +$
      (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*(tBez*beta*(-I + 2.0*!Pi*BW*tau) - !Pi*tau))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*$
      tau))*Sqrt(2.0)*beta^2.0 + I^(2.0*1/4.0)*exp((I*!Pi*(2.0*tBez*beta*BW + 1)^2.0)/(2.0*tBez*beta^2.0*BW) +$
      (tBez*(beta + 1))/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I + 2.0*!Pi*BW*tau)))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*$
      tau))*Sqrt(2.0)*beta^2.0 + I^(2.0*3.0/4.0)*exp((tBez*(beta + 1))/tau +$
      I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(tBez*beta*(I + 2.0*!Pi*BW*tau) - !Pi*tau))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*$
      Sqrt(2.0)*beta^2.0 + I^(2.0*3.0/4.0)*exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 4.0/beta +$
      1/(tBez*beta^2.0*BW)))*tBez^2.0*Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0)*beta^2.0 +$
      I^(2.0*3.0/4.0)*exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) +$
      1/(beta*BW*tau))*tBez^2.0*Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(-I + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0)*beta^2.0 +$
      I^(2.0*1/4.0)*exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) +$
      (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0)*beta^2.0 +$
      2.0*I^(2.0*3.0/4.0)*exp((4.0*tBez^2.0*BW*(beta + 2.0*I*!Pi*BW*tau + 1)*beta^2.0 + 3.0*I*!Pi*tau +$
      2.0*tBez*(6*I*!Pi*BW*tau*beta + beta))/(4.0*tBez*beta^2.0*BW*tau))*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(-I + 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      Sqrt(2.0)*beta^2.0 + 2.0*I^(2.0*1/4.0)*$
      exp((2.0*tBez^2.0*(4.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 2.0*!Pi*(beta + 1)*BW*tau + I)*beta^2.0 +$
      2.0*!Pi*tBez*tau*(6*I*!Pi*BW*tau + 1)*beta + 3.0*I*!Pi^2.0*tau^2.0)/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      Sqrt(2.0)*beta^2.0 - I^(2.0*1/4.0)*Sqrt(2.0)*exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) +$
      (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi*tBez*tau*$
      Erfz((I^(2.0*1/4.0)*(tBez*(beta - 2.0*I*!Pi*(beta - 1)*beta*BW*tau) - I*!Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta - I^(2.0*3.0/4.0)*Sqrt(2.0)*$
      exp((I*!Pi*(2.0*tBez*beta*BW + 1)^2.0)/(2.0*tBez*beta^2.0*BW) + (tBez*(beta + 1))/tau +$
      (I*tBez)/(2.0*!Pi*BW*tau^2.0))*!Pi*tBez*tau*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I + 2.0*!Pi*BW*tau)))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*$
      tau))*beta - I^(2.0*1/4.0)*Sqrt(2.0)*exp((tBez*(beta + 1))/tau +$
      I*!Pi*(2.0*tBez*BW + 4.0/beta + 1/(tBez*beta^2.0*BW)))*!Pi*tBez*tau*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*!Pi*BW*tau)))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*$
      beta + I^(2.0*3.0/4.0)*exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) +$
      1/(beta*BW*tau))*!Pi*tBez*tau*Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0)*beta +$
      I^(2.0*3.0/4.0)*exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) +$
      1/(beta*BW*tau))*!Pi*tBez*tau*Erfz((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I + 2.0*!Pi*(beta - 1)*BW*$
      tau)))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0)*beta +$
      I^(2.0*1/4.0)*exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) +$
      (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi*tBez*tau*$
      Erfz((I^(2.0*1/4.0)*(tBez*(2.0*I*!Pi*BW*tau*beta + beta) - I*!Pi*tau))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*$
      tau))*Sqrt(2.0)*beta + I^(2.0*3.0/4.0)*exp((I*!Pi*(2.0*tBez*beta*BW + 1)^2.0)/(2.0*tBez*beta^2.0*BW) +$
      (tBez*(beta + 1))/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0))*!Pi*tBez*tau*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0)*beta +$
      I^(2.0*1/4.0)*exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 4.0/beta + 1/(tBez*beta^2.0*BW)))*!Pi*tBez*tau*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0)*beta +$
      exp((2.0*tBez^2.0*(4.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 2.0*!Pi*(beta + 1)*BW*tau + I)*beta^2.0 +$
      2.0*!Pi*tBez*tau*(6*I*!Pi*BW*tau + 1)*beta + 3.0*I*!Pi^2.0*tau^2.0)/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*$
      (4.0 - 4.0*I)*(tBez^2.0*beta^2.0 + !Pi^2.0*tau^2.0)*Erfz(((1/2.0 + I/2.0)*Sqrt(tBez/BW))/(Sqrt(2.0*!Pi)*tau)) -$
      I^(2.0*3.0/4.0)*Sqrt(2.0)*exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) +$
      (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfz((I^(2.0*1/4.0)*(tBez*(beta - 2.0*I*!Pi*(beta - 1)*beta*BW*tau) - I*!Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) +$
      exp((4.0*tBez^2.0*BW*(beta + 2.0*I*!Pi*BW*tau + 1)*beta^2.0 + 3.0*I*!Pi*tau +$
      2.0*tBez*(6*I*!Pi*BW*tau*beta + beta))/(4.0*tBez*beta^2.0*BW*tau))*(4.0 - 4.0*I)*!Pi^2.0*tau^2.0*$
      Erfi(((1/2.0 + I/2.0)*Sqrt(tBez/BW))/(Sqrt(2.0*!Pi)*tau)) - 2.0*I^(2.0*1/4.0)*Sqrt(2.0)*$
      exp((2.0*tBez^2.0*(4.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 2.0*!Pi*(beta + 1)*BW*tau + I)*beta^2.0 +$
      2.0*!Pi*tBez*tau*(6*I*!Pi*BW*tau + 1)*beta + 3.0*I*!Pi^2.0*tau^2.0)/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*$
      tau^2.0*Erfi((I^(2.0*1/4.0)*tBez*(-I + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau)) -$
      2.0*I^(2.0*3.0/4.0)*Sqrt(2.0)*exp((4.0*tBez^2.0*BW*(beta + 2.0*I*!Pi*BW*tau + 1)*beta^2.0 + 3.0*I*!Pi*tau +$
      2.0*tBez*(6*I*!Pi*BW*tau*beta + beta))/(4.0*tBez*beta^2.0*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau)) -$
      8.0*exp((tBez^2.0*(8.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 4.0*!Pi*(beta + 2.0)*BW*tau + I)*beta^2.0 +$
      2.0*!Pi*tBez*tau*(6*I*!Pi*BW*tau + 1)*beta + 3.0*I*!Pi^2.0*tau^2.0)/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*$
      tau^2.0*FresnelS(Sqrt(2.0)*Sqrt(tBez*BW)) + I^(2.0*1/4.0)*$
      exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi^2.0*$
      tau^2.0*Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*!Pi*BW*tau) - !Pi*tau))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*$
      tau))*Sqrt(2.0) + I^(2.0*1/4.0)*exp((tBez*(beta + 1))/tau +$
      I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfz((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0) + I^(2.0*3.0/4.0)*$
      exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) + (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta +$
      1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfz((I^(2.0*1/4.0)*(tBez*(2.0*I*!Pi*BW*tau*beta + beta) - I*!Pi*tau))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*$
      tau))*Sqrt(2.0) + I^(2.0*1/4.0)*exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) +$
      (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*(tBez*beta*(-I + 2.0*!Pi*BW*tau) - !Pi*tau))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*$
      tau))*Sqrt(2.0) + I^(2.0*3.0/4.0)*exp((tBez*(beta + 1))/tau +$
      I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(tBez*beta*(I + 2.0*!Pi*BW*tau) - !Pi*tau))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*$
      Sqrt(2.0) + I^(2.0*3.0/4.0)*exp((tBez*(beta + 1))/tau + I*!Pi*(2.0*tBez*BW + 2.0/beta + 1/(tBez*beta^2.0*BW)) +$
      1/(beta*BW*tau))*!Pi^2.0*tau^2.0*Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(-I + 2.0*!Pi*(beta - 1)*BW*$
      tau)))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0) +$
      I^(2.0*1/4.0)*exp(((beta + 1)*tBez)/tau + (I*tBez)/(2.0*!Pi*BW*tau^2.0) +$
      (1/2.0)*I*!Pi*(4.0*tBez*BW + 8.0/beta + 1/(tBez*beta^2.0*BW)) + 1/(beta*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*Sqrt(2.0) + 2.0*I^(2.0*3.0/4.0)*$
      exp((4.0*tBez^2.0*BW*(beta + 2.0*I*!Pi*BW*tau + 1)*beta^2.0 + 3.0*I*!Pi*tau +$
      2.0*tBez*(6*I*!Pi*BW*tau*beta + beta))/(4.0*tBez*beta^2.0*BW*tau))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(-I + 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      Sqrt(2.0) + 2.0*I^(2.0*1/4.0)*exp((2.0*tBez^2.0*(4.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 2.0*!Pi*(beta + 1)*BW*tau + I)*$
      beta^2.0 + 2.0*!Pi*tBez*tau*(6*I*!Pi*BW*tau + 1)*beta + 3.0*I*!Pi^2.0*tau^2.0)/$
      (4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*tau^2.0*Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*!Pi*(beta - 1)*BW*tau))/$
      (2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*Sqrt(2.0)) -$
      8.0*exp((1/4.0)*((2.0*I*tBez)/(!Pi*BW*tau^2.0) + I*!Pi*(8.0*tBez*BW + 12.0/beta + 3.0/(tBez*beta^2.0*BW)) +$
      (8.0*tBez*beta*BW + 2.0)/(beta*BW*tau)))*!Pi^2.0*tBez*tau^2.0*$
      FresnelS(Sqrt(2.0)*(tBez - tBez*beta)*Sqrt(BW/tBez)))*gBez[ax])/$
      exp((2.0*tBez^2.0*(4.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 2.0*!Pi*(beta + 2.0)*BW*tau + I)*beta^2.0 +$
      2.0*!Pi*tBez*tau*(2.0*X*beta*BW + 6*I*!Pi*tau*BW + 1)*beta + 3.0*I*!Pi^2.0*tau^2.0)/$
      (4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))/(16.0*tBez*(tBez^2.0*beta^2.0 + !Pi^2.0*tau^2.0)*FresnelS(Sqrt(2.0)))))
  endfor

  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_FresnelT_fullmodel
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model for eddy current distortion and oscillation following a Fresnel test pulse
;
; Editing Information:
; Model expression commented does not account for windowing
function Bez_fit_time_FresnelT_oscillation,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0
  for j=0,n-1,4 do begin
    alpha = p[j]
    tau = p[j+1]
    f = p[j+2]
    phi0 = p[j+3]

    model+= imaginary((exp(-(X/tau) + 2.0*f*I*!Pi*X - (I*(2.0*tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 2.0*!Pi^2.0*BW^2.0*tau^2.0 - 2.0*I*!Pi*(beta + 1)*BW*tau + $
      4.0*f*!Pi*(-I + !Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*alpha*$
      (Sqrt(BW)*(I^(2.0*1/4.0)*Sqrt(2.0)*(exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) +$
      (I*tBez)/(4.0*!Pi*BW*tau^2.0))*I*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*I*exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*f^2.0*$
      !Pi^2.0*tBez^2.0*tau^2.0*Erfi((I^(2.0*3.0/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*!Pi*$
      tBez^2.0*tau*Erfi((I^(2.0*3.0/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*I*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*I*exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*I*exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*$
      tBez^2.0*Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*$
      f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*$
      BW*tau)))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*I*exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^$
      2.0 + 4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/$
      (4.0*!Pi*beta*BW*tau^2.0))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau - tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau - tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*I*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau - tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*I*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*I*exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*f^2.0*$
      !Pi^2.0*tBez^2.0*tau^2.0*Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*$
      tau)))/(2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*!Pi*$
      tBez^2.0*tau*Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*$
      tBez^2.0*Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*$
      f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*f*$
      I*!Pi*tBez^2.0*tau*Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      I*exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*f^2.0*I*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      4.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      16.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      16.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f*I*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      4.0*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*I*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      16.0*I*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 +$
      16.0*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      2.0*I*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      beta^2.0 + 8.0*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*$
      !Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f^2.0*I*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      beta^2.0 - 8.0*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*$
      !Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      beta^2.0 + 2.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*$
      !Pi*(-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*$
      tau))*beta^2.0 - 8.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*$
      tau + 4.0*f*!Pi*(-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*$
      (-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f^2.0*$
      !Pi^2.0*tBez^2.0*tau^2.0*Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau))/$
      (2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*beta^2.0 -$
      8.0*I*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*$
      tau))*beta^2.0 + 2.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*$
      tau + 4.0*f*!Pi*(-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*$
      (-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      beta^2.0 - 8.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*$
      !Pi*(-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f^2.0*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      beta^2.0 - 8.0*I*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*$
      f*!Pi*(-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*$
      BW*tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau))*$
      beta^2.0 - 2.0*I*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*$
      f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*tBez^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*$
      tau))*beta^2.0 + 8.0*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*$
      tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*$
      !Pi*BW*tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f^2.0*I*!Pi^2.0*tBez^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*$
      tau))*beta^2.0 - 8.0*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*$
      tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*$
      !Pi*BW*tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*f*!Pi*tBez^2.0*tau*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*$
      tau))*beta^2.0 + 2.0*exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) +$
      (I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*!Pi^2.0*tBez*tau^2.0*$
      Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta +$
      exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*I*!Pi*$
      tBez*tau*Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta -$
      2.0*I*exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^$
      2.0 + 4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/$
      (4.0*!Pi*beta*BW*tau^2.0))*f*!Pi^2.0*tBez*tau^2.0*$
      Erfz((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta +$
      exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*!Pi*$
      tBez*tau*Erfz((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta -$
      2.0*exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*$
      !Pi^2.0*tBez*tau^2.0*Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta -$
      I*exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*!Pi*tBez*$
      tau*Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta +$
      2.0*exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*f*I*!Pi^2.0*tBez*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta -$
      exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*!Pi*tBez*tau*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta -$
      2.0*exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*!Pi^2.0*tBez*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta -$
      I*exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*!Pi*tBez*tau*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta +$
      2.0*exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*f*!Pi^2.0*tBez*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau - tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta +$
      exp((I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (2.0*I*!Pi*f)/(beta*BW) + (3.0*tBez*f)/(BW*tau) +$
      I*!Pi*tBez*BW + (tBez*beta)/tau + (3.0*I*tBez)/(4.0*!Pi*BW*tau^2.0))*I*!Pi*tBez*tau*$
      Erfi((I^(2.0*1/4.0)*(!Pi*tau - tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta -$
      2.0*I*exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*f*!Pi^2.0*tBez*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta +$
      exp((I*((8.0*!Pi^2.0*(f + BW))/beta + (tBez*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 -$
      4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*(-I + 4.0*!Pi*BW*tau)*tau + 1))/tau^2.0 +$
      (2.0*!Pi^2.0)/(tBez*beta^2.0)))/(4.0*!Pi*BW))*!Pi*tBez*tau*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau))*beta +$
      exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f + (tBez*f)/(BW*tau) + I*!Pi*tBez*BW +$
      (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau + 1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*!Pi^2.0*$
      tau^2.0*Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) +$
      exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*I*$
      !Pi^2.0*tau^2.0*Erfz((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) +$
      exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*I*$
      !Pi*tau*(tBez*beta*(I + 2.0*f*!Pi*tau) - !Pi*tau)*$
      Erfz((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) - exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f +$
      (tBez*f)/(BW*tau) + I*!Pi*tBez*BW + (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau +$
      1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*!Pi^2.0*tau^2.0*$
      Erfz((I^(2.0*1/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) + exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f +$
      (tBez*f)/(BW*tau) + I*!Pi*tBez*BW + (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau +$
      1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*I*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau) - !Pi*tau))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) +$
      exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*$
      !Pi^2.0*tau^2.0*Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) + exp((3.0*I*!Pi*tBez*f^2.0)/BW + 4.0*I*!Pi*tBez*f +$
      (tBez*f)/(BW*tau) + I*!Pi*tBez*BW + (I*!Pi)/(2.0*tBez*beta^2.0*BW) + (tBez*beta)/tau +$
      1/(beta*BW*tau) + (I*tBez)/(4.0*!Pi*BW*tau^2.0))*I*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*(!Pi*tau + tBez*beta*(-I - 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) -$
      exp((4.0*!Pi*tau*(2.0*I*!Pi*BW*tau + 1) + tBez*beta*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 4.0*I*!Pi^2.0*BW^2.0*tau^2.0 +$
      4.0*!Pi*beta*BW*tau + 4.0*f*!Pi*(4.0*I*!Pi*BW*tau + 3.0)*tau + 3.0*I))/(4.0*!Pi*beta*BW*tau^2.0))*$
      !Pi^2.0*tau^2.0*Erfi((I^(2.0*1/4.0)*(!Pi*tau + tBez*beta*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau)))/$
      (2.0*Sqrt(!Pi)*beta*Sqrt(tBez*BW)*tau)) -$
      4.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau)) +$
      4.0*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*I*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau)) -$
      2.0*I*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau)) +$
      2.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*$
      tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau + 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*$
      tau)) + 2.0*exp((I*(tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*$
      f*!Pi*(-3.0*I + 4.0*!Pi*BW*tau)*tau + 3.0)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*$
      BW*tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*1/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*tau)) -$
      2.0*I*exp((I*(tBez^2.0*(12.0*f^2.0*!Pi^2.0*tau^2.0 + 4.0*!Pi^2.0*BW^2.0*tau^2.0 - 4.0*I*!Pi*beta*BW*tau + 4.0*f*!Pi*$
      (-I + 4.0*!Pi*BW*tau)*tau + 1)*beta^2.0 + 2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*$
      beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*tau^2.0*$
      Erfi((I^(2.0*3.0/4.0)*tBez*(I + 2.0*f*!Pi*tau - 2.0*!Pi*(beta - 1)*BW*tau))/(2.0*Sqrt(!Pi)*Sqrt(tBez*BW)*$
      tau)))*(Cos(phi0) + I*Sin(phi0)) -$
      8.0*exp((2.0*tBez^2.0*(4.0*f^2.0*I*!Pi^2.0*tau^2.0 + 2.0*I*!Pi^2.0*BW^2.0*tau^2.0 + 2.0*!Pi*(beta + 1)*BW*tau +$
      4.0*f*!Pi*(I*!Pi*BW*tau + 1)*tau + I)*beta^2.0 + 2.0*!Pi*tBez*tau*(2.0*I*!Pi*BW*tau +$
      2.0*f*I*!Pi*tau + 1)*beta + I*!Pi^2.0*tau^2.0)/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*!Pi^2.0*tau^2.0*$
      FresnelS(Sqrt(2.0)*Sqrt(tBez*BW))*(Cos(phi0) + I*Sin(phi0))) -$
      8.0*exp((I*(2.0*tBez^2.0*(4.0*f^2.0*!Pi^2.0*tau^2.0 + 2.0*!Pi^2.0*BW^2.0*tau^2.0 - 2.0*I*!Pi*BW*tau +$
      4.0*f*!Pi*(-I + !Pi*(beta + 1)*BW*tau)*tau + 1)*beta^2.0 +$
      2.0*!Pi*tBez*tau*(-I + 2.0*f*!Pi*tau + 2.0*!Pi*BW*tau)*beta + !Pi^2.0*tau^2.0))/(4.0*!Pi*tBez*beta^2.0*BW*tau^2.0))*$
      !Pi^2.0*Sqrt(BW)*tau^2.0*FresnelS(Sqrt(2.0)*(tBez - tBez*beta)*Sqrt(BW/tBez))*$
      (Cos(phi0) + I*Sin(phi0)))*gBez[ax])/$
      (16.0*Sqrt(BW)*(tBez^2.0*beta^2.0*(I + 2.0*f*!Pi*tau)^2.0 - !Pi^2.0*tau^2.0)*FresnelS(Sqrt(2.0))))
  endfor
  return, (y-model)/err

end


; Subroutine name: Bez_fit_time_SincT
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Time dependent eddy current distortion model for Sinc test pulse with Tukey window
;
; Editing Information:
; Model expression commented does not account for windowing

function Bez_fit_time_SincT_eddy,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0

  for j=0,n-1,2 do begin
    alpha = p[j]
    tau = p[j+1]

    model+= real_part((1.0/(8.0*(tBez^2.0*beta^2.0 + !DPI^2.0*tau^2.0)))*(exp(-((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta) -$
      (X + tBez + tBez*beta)/tau)*alpha*$
      (-2.0*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez^2.0*beta^2.0 -$
      2.0*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^3*tau^2.0 -$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      Ei((tBez*eta*(-1.0 - 2.0*I*!DPI*BW*tau))/tau) -$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei((tBez*eta*(-1.0 - 2.0*I*!DPI*BW*tau))/tau) +$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      Ei((tBez*eta*(-1.0 + 2.0*I*!DPI*BW*tau))/tau) +$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei((tBez*eta*(-1.0 + 2.0*I*!DPI*BW*tau))/tau) +$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      Ei(-((I*tBez*(-1.0 + beta + eta)*(-I + 2.0*!DPI*BW*tau))/tau)) +$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei(-((I*tBez*(-1.0 + beta + eta)*(-I + 2.0*!DPI*BW*tau))/tau)) -$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*!DPI*BW*tau))/tau) -$
      4.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*!DPI*BW*tau))/tau) +$
      exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei(-((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei(-((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei(-((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei(-((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei(-((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei(-((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei(-((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei(-((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      Ei((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*!DPI*BW*tau))/tau) +$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*!DPI*BW*tau))/tau) -$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*!DPI*BW*tau))/tau) -$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*!DPI*BW*tau))/tau) -$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU((tBez*(-1.0 + eta)*(1.0 + 2.0*I*!DPI*BW*tau))/tau) -$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU((tBez*(-1.0 + eta)*(1.0 + 2.0*I*!DPI*BW*tau))/tau) +$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU((tBez*(-1.0 + beta + eta)*(1.0 + 2.0*I*!DPI*BW*tau))/tau) +$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU((tBez*(-1.0 + beta + eta)*(1.0 + 2.0*I*!DPI*BW*tau))/tau) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta - I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta - I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta + I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta + I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta - I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta - I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(((-1.0 + beta + eta)*(tBez*beta - I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(((-1.0 + beta + eta)*(tBez*beta - I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta + I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(((-1.0 + eta)*(tBez*beta + I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau))/(beta*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      GammaU(-((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      GammaU(-((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog(-1.0 + beta + eta) - 2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta +$
      (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*alog(-1.0 + beta + eta) -$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog(-1.0 - 2.0*I*!DPI*BW*tau) - 2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta +$
      (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*alog(-1.0 - 2.0*I*!DPI*BW*tau) +$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog((-1.0 + beta + eta)*(1.0 + 2.0*I*!DPI*BW*tau)) +$
      2.0*I*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog((-1.0 + beta + eta)*(1.0 + 2.0*I*!DPI*BW*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta - I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta - I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau - 2.0*I*!DPI*tBez*beta*BW*tau)) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta - I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau)) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta - I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog((-1.0 + beta + eta)*(tBez*beta + I*!DPI*tau + 2.0*I*!DPI*tBez*beta*BW*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog(I*!DPI*tau + tBez*beta*(-1.0 - 2.0*I*!DPI*BW*tau)) +$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog(I*!DPI*tau + tBez*beta*(-1.0 - 2.0*I*!DPI*BW*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog((-I)*!DPI*tau + tBez*beta*(-1.0 + 2.0*I*!DPI*BW*tau)) -$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog((-I)*!DPI*tau + tBez*beta*(-1.0 + 2.0*I*!DPI*BW*tau)) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog((-I)*(!DPI*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau))) +$
      I*exp((2.0*I*!DPI*(tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog((-I)*(!DPI*tau + tBez*beta*(-I + 2.0*!DPI*BW*tau))) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*tBez^2.0*beta^2.0*$
      alog(I*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau))) -$
      I*exp((2.0*I*!DPI*(1.0 + tBez*beta*BW*(-1.0 + eta)))/beta + (tBez*(beta + eta))/tau)*!DPI^2.0*tau^2.0*$
      alog(I*(!DPI*tau + tBez*beta*(I + 2.0*!DPI*BW*tau))) +$
      4.0*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta)*(exp(tBez/tau) + exp((tBez*(1.0 + beta))/tau))*!DPI^2.0*$
      tau^2.0*Si(2.0*!DPI*tBez*BW*eta) +$
      4.0*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + (tBez*(1.0 + beta))/tau)*!DPI^2.0*tau^2.0*$
      Si(2.0*!DPI*tBez*BW - 2.0*!DPI*tBez*BW*eta) +$
      4.0*exp((I*!DPI*(1.0 + 2.0*tBez*beta*BW*(-1.0 + eta) + eta))/beta + tBez/tau)*!DPI^2.0*tau^2.0*$
      Si(2.0*!DPI*tBez*BW - 2.0*!DPI*tBez*beta*BW - 2.0*!DPI*tBez*BW*eta))*gBez(ax)))/scale_sinc
  endfor
  return, (y-model)/err

end

; Subroutine name: Bez_fit_time_SincT_fullmodel
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Model for eddy current distortion and oscillation following a Sinc test pulse
;
; Editing Information:
; Model expression commented does not account for windowing
function Bez_fit_time_SincT_oscillation,p,X=x,Y=y,ERR=err

  common Bez_params
  common scan_data
  ax = project.procPramArray[project.ci].Bez_axis

  I = dcomplex(0,1)
  n = n_elements(p)
  model = 0
  for j=0,n-1,4 do begin
    alpha = p[j]
    tau = p[j+1]
    f = p[j+2]
    phi0 = p[j+3]

    model+= imaginary(-((1.0/(8.0*((-!DPI^2.0)*tau^2.0 + tBez^2.0*beta^2.0*(I + 2.0*f*!DPI*tau)^2.0)))*$
      (exp(-((X*beta + 2.0*I*(!DPI*eta*tau + tBez*beta*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      I*phi0)*alpha*$
      (-4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))/tau) -$
      16.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))/tau) -$
      4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))/tau) +$
      16.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))/tau) +$
      4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))/tau) +$
      16.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))/$
      tau) + 4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI^2.0*tau^2.0*$
      Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))/tau) -$
      16.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau -$
      2.0*!DPI*BW*tau))/tau) +$
      4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/tau) +$
      16.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/tau) +$
      4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/tau) -$
      16.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*Ei((I*tBez*eta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/tau) -$
      4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/tau) -$
      16.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/$
      tau) - 4.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI^2.0*tau^2.0*$
      Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/tau) +$
      16.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*Ei((I*tBez*(-1.0 + beta + eta)*(I + 2.0*f*!DPI*tau +$
      2.0*!DPI*BW*tau))/tau) -$
      exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI*tBez*beta*tau*Ei($
      (I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*Ei($
      (I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      2.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*!DPI^2.0*tBez*beta*tau^2.0*$
      Ei((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/$
      (beta*tau)) + exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/$
      (beta*tau)) - I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI^2.0*tau^2.0*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/$
      (beta*tau)) - 2.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*!DPI^2.0*tBez*beta*tau^2.0*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/$
      (beta*tau)) + exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*Ei($
      (I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      2.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI^2.0*tBez*beta*tau^2.0*Ei($
      (I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI*tBez*beta*tau*Ei((I*(-1.0 + beta + eta)*(!DPI*tau +$
      tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*Ei((I*(-1.0 + beta + eta)*$
      (!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      2.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI^2.0*tBez*beta*tau^2.0*Ei($
      (I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI*tBez*beta*tau*Ei($
      (I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*Ei($
      (I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      2.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*!DPI^2.0*tBez*beta*tau^2.0*$
      Ei((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/$
      (beta*tau)) - exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/$
      (beta*tau)) + I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI^2.0*tau^2.0*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/$
      (beta*tau)) + 2.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*!DPI^2.0*tBez*beta*tau^2.0*$
      Ei((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/$
      (beta*tau)) - exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*!DPI*tBez*beta*tau*$
      Ei((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*Ei($
      (I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      2.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI^2.0*tBez*beta*tau^2.0*Ei($
      (I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI*tBez*beta*tau*Ei((I*(-1.0 + beta + eta)*(!DPI*tau +$
      tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*Ei((I*(-1.0 + beta + eta)*$
      (!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) -$
      2.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI^2.0*tBez*beta*tau^2.0*Ei($
      (I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau)) +$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/tau) +$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/tau) +$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/tau) -$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/$
      tau) - 2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*tBez^2.0*beta^2.0*$
      GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/tau) -$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/tau) -$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/tau) +$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))/$
      tau) - 2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*tBez^2.0*beta^2.0*$
      GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/tau) -$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/tau) -$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/tau) +$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU((tBez*(-1.0 + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/$
      tau) + 2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*tBez^2.0*beta^2.0*$
      GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/tau) +$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/tau) +$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/tau) -$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU((tBez*(-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))/$
      tau) + I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*tBez^2.0*beta^2.0*$
      GammaU(-((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      -((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*GammaU($
      ((-1.0 + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/(beta*tau)) +$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      ((-1.0 + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/(beta*tau)) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU($
      ((-1.0 + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/(beta*tau)) -$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      GammaU(((-1.0 + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/(beta*tau)) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*GammaU(((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/(beta*tau)) -$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      ((-1.0 + beta + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/(beta*tau)) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU(((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/(beta*tau)) +$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      GammaU(((-1.0 + beta + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)))/$
      (beta*tau)) - I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*tBez^2.0*beta^2.0*$
      GammaU(-((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      -((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      GammaU(-((I*(-1.0 + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      -((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      GammaU(-((I*(-1.0 + beta + eta)*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/$
      (beta*tau))) - I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*tBez^2.0*beta^2.0*$
      GammaU(-((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      -((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*GammaU($
      -((I*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/(beta*tau))) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog(tBez) -$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog(tBez) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog(tBez) +$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog(tBez) -$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau)) -$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau)) -$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau)) +$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau)) -$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*alog(-1.0 + 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau) -$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*alog(-1.0 + 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau) -$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*alog(-1.0 + 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau) +$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog(-1.0 + 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau) +$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)) +$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)) +$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)) -$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog((-1.0 + beta + eta)*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau)) +$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      tBez^2.0*beta^2.0*alog(-1.0 + 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau) +$
      8.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f*$
      !DPI*tBez^2.0*beta^2.0*tau*alog(-1.0 + 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau) +$
      2.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*alog(-1.0 + 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau) -$
      8.0*I*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*$
      f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog(-1.0 + 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) +$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) -$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog(I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) +$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog(I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog(I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) -$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      alog(I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau - 2.0*!DPI*BW*tau))) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))) +$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))) +$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))) -$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      alog((-1.0 + beta + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau - 2.0*I*!DPI*BW*tau))) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog((-1.0 + beta + eta)*((-I)*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) -$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog((-1.0 + beta + eta)*((-I)*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog((-1.0 + beta + eta)*((-I)*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) +$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog((-1.0 + beta + eta)*((-I)*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) -$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog((-1.0 + beta + eta)*(I*!DPI*tau +$
      tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) +$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      alog((-1.0 + beta + eta)*(I*!DPI*tau + tBez*beta*(1.0 - 2.0*I*f*!DPI*tau + 2.0*I*!DPI*BW*tau))) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) -$
      4.0*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) -$
      I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) +$
      4.0*I*exp(I*((!DPI*(-1.0 + 2.0*f*X*beta + 3.0*eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*$
      alog(I*((-!DPI)*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*tBez^2.0*beta^2.0*alog((I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/tBez) -$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog((I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/tBez) -$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog((I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/tBez) +$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog((I*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/$
      tBez) + I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta +$
      (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/tau))*tBez^2.0*beta^2.0*$
      alog((-I)*(-1.0 + beta + eta)*(!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) +$
      4.0*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f*!DPI*tBez^2.0*beta^2.0*tau*alog((-I)*(-1.0 + beta + eta)*$
      (!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) +$
      I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*!DPI^2.0*tau^2.0*alog((-I)*(-1.0 + beta + eta)*$
      (!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) -$
      4.0*I*exp(I*((!DPI*(1.0 + 2.0*f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 4.0*!DPI*BW*tau))/$
      tau))*f^2.0*!DPI^2.0*tBez^2.0*beta^2.0*tau^2.0*alog((-I)*(-1.0 + beta + eta)*$
      (!DPI*tau + tBez*beta*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))) +$
      4.0*exp(I*((2.0*!DPI*(f*X*beta + eta))/beta +$
      (tBez*(I*beta + 2.0*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau)))/tau))*$
      (exp(2.0*I*f*!DPI*tBez*beta) + exp((tBez*beta)/tau))*!DPI^2.0*tau^2.0*Si(2.0*!DPI*tBez*BW*eta) +$
      4.0*exp(2.0*I*((!DPI*(f*X*beta + eta))/beta + (tBez*(-1.0 + eta)*(I + 2.0*f*!DPI*tau + 2.0*!DPI*BW*tau))/tau))*$
      !DPI^2.0*tau^2.0*Si(2.0*!DPI*tBez*BW - 2.0*!DPI*tBez*BW*eta) +$
      4.0*exp(I*((2.0*!DPI*eta)/beta + 2.0*f*!DPI*(X + tBez*(-2.0 + beta + 2.0*eta)) +$
      (tBez*(I*beta + 2.0*(-1.0 + eta)*(I + 2.0*!DPI*BW*tau)))/tau))*!DPI^2.0*tau^2.0*$
      Si(2.0*!DPI*tBez*BW - 2.0*!DPI*tBez*beta*BW - 2.0*!DPI*tBez*BW*eta))*gBez[ax])))/scale_sinc

  endfor
  return, (y-model)/err

end

; Subroutine name: load_state2_data_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
; data - Pointer containing the processed state 1 eddy current data
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To prepare the eddy current field for display after rotation/flip/zoom/thresholding
;
; Editing Information:

function load_state2_data_Bez, data

  common scan_data

  if ~ptr_valid(data) then return, !NULL

  scheme = project.imndarray[project.ci].scheme
  xzoom = project.procpramArray[project.ci].x_zoom  ; x-zoom factor in MAS main widget
  yzoom = project.procpramArray[project.ci].y_zoom  ; y-zoom factor in MAS main widget

  N = size(*data,/Dimension)
  axis = project.procPramArray[project.ci].Bez_axis

  temp = *data
  if ptr_valid(Bez) then ptr_free,Bez
  Bez = ptr_new(make_array(xzoom*N[0],yzoom*N[1],N[2],N[3],2,value = 0,/double))

  for i = 0,N[3]-1 do begin
    (*Bez)[*,*,*,i,0] = congrid(temp[*,*,*,i,axis],xzoom*N[0],yzoom*N[1],N[2],interp = 0)   ; Zooming
    (*Bez)[*,*,*,i,1] = congrid(temp[*,*,*,i,3],xzoom*N[0],yzoom*N[1],N[2],interp = 0)   ; Zooming
  endfor
  mask_Bez,Bez
  rotate_Bez, Bez
  flip_Bez,Bez

  return, Bez

end

; Subroutine name: TimeSeries
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To fit the eddy current magnetic field or mechanical oscillations over time and space
;
; Editing Information:

function Smoothing, data, kernel_size

  kernel = make_Array(kernel_size,value = 1)/float(kernel_size)
  return, convol(data,kernel,/EDGE_ZERO)

end

; Subroutine name: orient_Bez_data
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To transform Bez from logical to physical gradient axes
;
; Editing Information:
pro orient_Bez_data, Bez

  M = log2physgradtransform()
  
  ; Rotating the vector
  N = size(*Bez,/Dimensions)
  for k = 0,N[3]-1 do begin
    Xraster = make_array(N[0]*N[1]*N[2],3)
    for i = 0,2 do Xraster[*,i] = reform((*Bez)[*,*,*,k,i],N[0]*N[1]*N[2],1)  ; Rastering the velocity field (3xRxPxS) into 2D (3xRPS) matrix
    Xrot = M##Xraster                                                         ; Rotating the rastered array
    for i = 0,2 do (*Bez)[*,*,*,k,i] = reform(Xrot[*,i],N[0],N[1],N[2])       ; Reforming the rotated array back to 3D  
  endfor
 
end

; Subroutine name: get_Bez_params
; Created by: Magdoom Kulam
; Calling Information:
;
; T       - Time between the test gradient pulse and echo (ms)
; axis    - Test gradient axis (X,Y,Z)
; shape   - Shape of the test gradient, trapezoid by default (Other : Sinusoid, Gauss, Shifted cosine)
; gBez    - Strength of the test gradient (mT/m)
; tBez    - Duration of the test gradient (ms)
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To retrieve the eddy current measurement related parameters for the selected scan
;
; Editing Information:

pro get_Bez_params, T, axis, shape, gBez, tBez, tasc = tasc, tplateau = tplateau, tdesc = tdesc, BW, beta = beta, scheme, eta = eta, scale_sinc = scale_sinc

  common scan_data
  CI = project.ci

  orient  = project.imndarray[CI].orientation                  ; Slice orientation
  ; Echo time (ms)
  if ptr_valid(project.imndarray[CI].echo_time_ptr) then TE = *project.imndarray[CI].echo_time_ptr $
    else TE = project.imndarray[CI].echo_time                   
  shape   = project.imndarray[CI].tgrshape                     ; Shape of test gradient pulse
  if shape eq 'Asymmetric Trapezoid' then tBez = tasc+tdesc+tplateau else  tBez = project.imndarray[CI].tgrtime
  T = []
  for i = 0,n_elements(TE)-1 do T = [T,TE[i]+*project.imndarray[CI].tgrdelay]
  T = T[sort(T)]
  gBez    = (*project.imndarray[CI].tgramp)*10                   ; Strength of test gradient pulse (mT/m)  
  tgraxis = project.imndarray[CI].tgraxis
  rotate_Direction = project.procpramArray[CI].rotate_Direction  ; Rotation direction in MAS main widget
  
  scheme = project.imndarray[CI].scheme
  case scheme of
    'f'  : scheme = 'Four Point'
            
    's'  : scheme = 'Six point'
            
    'h'  : scheme = 'Hadamard'
                
    else :  scheme = 'Test'
            
  endcase
    
  case shape of
    's'   : begin
            shape = 'Sinusoid'
            BW = 1/float(tBez)
            end
            
    'c'   : begin
            shape = 'ShiftedCosine'
            BW = 1/float(tBez)
            end
           
    'g'   : begin
            shape = 'Gauss'
            BW = 9.6*sqrt(2*!Pi)/float(tBez)
            end
    
    'i'   : begin
            shape = 'Sinc-Tukey'
            BW    = 1e-3*project.imndarray[CI].BWSinc                  ; Excitation bandwidth shifted Sine integral pulse
            beta  = project.imndarray[CI].betaSinc                     ; Tukey window factor
            eta  = project.imndarray[CI].etaSinc                       ; Shift fraction
            scale_sinc = project.imndarray[CI].scaleSinc
            end
            
    't'   : begin
            shape = 'Symmetric Trapezoid'
            max_slew_rate = 10*project.imndarray[CI].max_slew_rate      ; Maximum slew rate (mT/m.us)
            trise = float(gBez[project.procPramArray[CI].Bez_axis])/float(max_slew_rate)
            ; Corrected rise time to account for 4 us delay
            if ((trise mod 4) eq 0) then begin
                trise_corr = trise
            endif else begin
                trise_corr = 4*(1+floor(trise/4.0))
            endelse
            tasc = 1e-3*trise_corr                                   ; Ascending ramp time for the trapezoidal gradient (ms)
            tdesc = 1e-3*trise_corr                                  ; Descending ramp time for the trapezoidal gradient (ms)
            tplateau = tBez-tasc-tdesc                               ; Plateau time for the trapezoidal gradient (ms)
            BW = 1/float(tdesc)   
            end
    'a'   : begin
            shape = 'Asymmetric Trapezoid'
            tasc  = project.imndarray[CI].tasc                       ; Ascending ramp time for the trapezoidal gradient (ms)
            tdesc = project.imndarray[CI].tdesc                      ; Descending ramp time for the trapezoidal gradient (ms)
            tplateau = tBez-tasc-tdesc                               ; Plateau ramp time for the trapezoidal gradient (ms)
            BW = 1/float(tdesc)
            end
    'f'   : begin
            shape = 'Fresnel-Tukey'    
            BW    = 1e-3*project.imndarray[CI].BWFresnel             ; Excitation bandwidth of Fresnel pulse
            beta  = project.imndarray[CI].betaFresnel                ; Tukey window factor
            end
  endcase
  trans_mat = log2physgradtransform()
  
  case (strsplit(orient,/extract))[0] of

    'trans' : begin
      case tgraxis of
        'r'  : if rotate_Direction eq 0 or 2 then axis = 'Y (Left-Right)' else axis = 'Y (Up-Down)'
        'p'  : if rotate_Direction eq 0 or 2 then axis = 'X (Up-Down)' else axis = 'X (Left-Right)'
        's'  : axis = 'Z (In-Out)'
        'g'  : axis = 'XYZ'
        'b'  : axis = 'B0'
        else : axis = 'Test'
      endcase
    end

    'trans90': begin
      case tgraxis of
        'r'  : if rotate_Direction eq 0 or 2 then axis = 'X (Left-Right)' else axis = 'X (Up-Down)'
        'p'  : if rotate_Direction eq 0 or 2 then axis = 'Y (Up-Down)' else axis = 'Y (Left-Right)'
        's'  : axis = 'Z (In-Out)'
        'g'  : axis = 'XYZ'
        'b'  : axis = 'B0'
        else : axis = 'Test'
      endcase
    end

    'cor'    : begin
      case tgraxis of
        'r'  : if rotate_Direction eq 0 or 2 then axis = 'Z (Left-Right)' else axis = 'Z (Up-Down)'
        'p'  : if rotate_Direction eq 0 or 2 then axis = 'X (Up-Down)' else axis = 'X (Left-Right)'
        's'  : axis = 'Y (In-Out)'
        'g'  : axis = 'XYZ'
        'b'  : axis = 'B0'
        else : axis = 'Test'
      endcase
    end

    'cor90'  : begin
      case tgraxis of
        'r'  : if rotate_Direction eq 0 or 2 then axis = 'X (Left-Right)' else axis = 'X (Up-Down)'
        'p'  : if rotate_Direction eq 0 or 2 then axis = 'Z (Up-Down)' else axis = 'Z (Left-Right)'
        's'  : axis = 'Y (In-Out)'
        'g'  : axis = 'XYZ'
        'b'  : axis = 'B0'
        else : axis = 'Test'
      endcase
    end

    'sag'    : begin
      case tgraxis of
        'r'  : if rotate_Direction eq 0 or 2 then axis = 'Z (Left-Right)' else axis = 'Z (Up-Down)'
        'p'  : if rotate_Direction eq 0 or 2 then axis = 'Y (Up-Down)' else axis = 'Y (Left-Right)'
        's'  : axis = 'X (In-Out)'
        'g'  : axis = 'XYZ'
        'b'  : axis = 'B0'
        else : axis = 'Test'
      endcase
    end

    'sag90'  :  begin
      case tgraxis of
        'r'  : if rotate_Direction eq 0 or 2 then axis = 'Y (Left-Right)' else axis = 'Y (Up-Down)'
        'p'  : if rotate_Direction eq 0 or 2 then axis = 'Z (Up-Down)' else axis = 'Z (Left-Right)'
        's'  : axis = 'X (In-Out)'
        'g'  : axis = 'XYZ'
        'b'  : axis = 'B0'
        else : axis = 'Test'
      endcase
    end
    else    : begin
              axis = '---'
              end
  endcase
 gBez = abs(trans_mat##gBez) 
      
end




; Subroutine name: mask_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
; Bez  - Pointer containing Bez
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To apply image and field threshold mask
;
; Editing Information:

pro mask_Bez, Bez

common scan_data
common common_widgets

CI = project.ci
Ithreshold = project.procPramArray[CI].Bez_Ithreshold

N = size(*Bez,/Dimension)                          ; State 2 dimensions

for i = 0,N[3]-1 do begin  
  Imask = (*Bez)[*,*,*,i,1] ge max((*Bez)[*,*,*,i,1])*float(Ithreshold)/100.0
  (*Bez)[*,*,*,i,0] *= Imask
  (*Bez)[*,*,*,i,1] *= Imask  
endfor

end


; Subroutine name: flip_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
;  Bez - Pointer containing 3D Bez field along with the MR image.
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To flip Bez as specified in the main MAS GUI
;
; Editing Information:

pro flip_Bez,Bez

common scan_data

flip_Direction = project.procPramArray[project.ci].flip_Direction
if flip_Direction eq 0 then return
N = size(*Bez,/Dimension)

; Flipping the image
for i = 0,N[2]-1 do begin
    for j = 0,N[3]-1 do begin
      (*Bez)[*,*,i,j,0] = reverse((*Bez)[*,*,i,j,0],flip_Direction)
      (*Bez)[*,*,i,j,1] = reverse((*Bez)[*,*,i,j,1],flip_Direction)
  endfor
endfor

end

; Subroutine name: rotate_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
;  Bez - Pointer containing 3D Bez field along with the MR image.
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To rotate Bez to the angle specified in the main MAS GUI
;
; Editing Information:

pro rotate_Bez, Bez

common scan_data

rotate_Direction = project.procpramArray[project.ci].rotate_Direction  ; Rotation direction in MAS main widget
if rotate_Direction eq 0 then return

N = size(*Bez,/Dimension)

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

temp = *Bez
Bez = ptr_new(make_array(N[0],N[1],N[2],N[3],2,value = 0,/double))   ; Resizing the Bez pointer based on rotation

; Rotating the image
for i = 0,N[2]-1 do begin
  for j = 0,N[3]-1 do begin
      (*Bez)[*,*,i,j,0] = rotate(temp[*,*,i,j,0],rotate_Direction)
      (*Bez)[*,*,i,j,1] = rotate(temp[*,*,i,j,1],rotate_Direction)
  endfor
endfor

end



; Subroutine name: unwrap_phase_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
; Bez - Pointer containing 3D Bez field along with the MR image.
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To unwrap the eddy current field using Fourier based method
;
; Editing Information:

pro unwrap_phase_Bez, Bez, scheme, T

common scan_data

N = size(*Bez,/Dimension)

if scheme eq 'Four Point' or scheme eq 'Six Point' or scheme eq 'Hadamard' then Naxes = 3 else Naxes = 1

progressbar = Obj_New('progressbar', Color = 'red', Text = 'Phase Unwrapping', /NOCANCEL)
progressbar -> Start
counter = 0UL
for i = 0,N[3]-1 do begin
    for j = 0, Naxes-1 do begin
      phase_image = reform((*Bez)[*,*,*,i,j])
      phase_unwrap_3d, phase_image, tolerance = 0.1
      (*Bez)[*,*,*,i,j] = phase_image
      counter++
      progressBar -> Update, 100*float(counter)/float(Naxes)
    endfor   
endfor     
progressbar -> Destroy

end

; Subroutine name: Calculate_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To calculate the normalized eddy current magnetic field from phase via complex division
;
; Editing Information:
; 1) Fixed sign error in Bez (Magdoom, 1/30/17)

pro Calculate_Bez, CI

common scan_data
common Bez_params

gamma = 2*!Pi*42.58E-3                      ; Proton gyromagnetic ratio (rad/ms.uT)

;; Loading complex data
set_signal_type, 9, ci = CI
mas_load_state_1, ci = CI

orient  = project.imndarray[CI].orientation          ; Slice orientation
data = *project.dataarray[CI].state1
if ptr_valid(project.imndarray[CI].echo_time_ptr) then TE = *project.imndarray[CI].echo_time_ptr $
else TE = project.imndarray[CI].echo_time
Ndelays = n_elements(*project.imndarray[CI].tgrdelay)
N = size(data,/Dimension)

case scheme of
  'Four Point' : begin
                 nskip = 4
                 xi = 1
                 end
  
  'Six Point' : begin
                nskip = 6
                xi = 2
                end
                
   'Hadamard'  : begin
                 nskip = 4
                 xi = 4
                 end
                                     
   else        : begin
                 nskip = 2
                 xi = 2
                 end
endcase

Nt = N[3]/nskip
if ptr_valid(project.dataarray[CI].Bez) then ptr_free,project.dataarray[CI].Bez
project.dataarray[CI].Bez = ptr_new(make_array(N[0],N[1],N[2],Nt,4,value = 0,/double))
k = 0
for i = 0,Nt-1 do begin
  case scheme of       
  'Four Point' : begin
                 (*project.dataarray[CI].Bez)[*,*,*,i,0] = atan(data[*,*,*,k+3]/data[*,*,*,k],/Phase)
                 (*project.dataarray[CI].Bez)[*,*,*,i,1] = atan(data[*,*,*,k+3]/data[*,*,*,k+1],/Phase)
                 (*project.dataarray[CI].Bez)[*,*,*,i,2] = atan(data[*,*,*,k+3]/data[*,*,*,k+2],/Phase)                    
                 (*project.dataarray[CI].Bez)[*,*,*,i,3] = mean(abs(data[*,*,*,k:k+3]),dimension = 4)               
                 end
      
   'Six Point' : begin              
                 (*project.dataarray[CI].Bez)[*,*,*,i,0] = atan(data[*,*,*,k+1]/data[*,*,*,k],/Phase)
                 (*project.dataarray[CI].Bez)[*,*,*,i,1] = atan(data[*,*,*,k+3]/data[*,*,*,k+2],/Phase)
                 (*project.dataarray[CI].Bez)[*,*,*,i,2]=  atan(data[*,*,*,k+5]/data[*,*,*,k+4],/Phase)
                 (*project.dataarray[CI].Bez)[*,*,*,i,3] = mean(abs(data[*,*,*,k:k+5]),dimension = 4)                 
                 end             

   'Hadamard' : begin  
                (*project.dataarray[CI].Bez)[*,*,*,i,0] = atan(data[*,*,*,k]*data[*,*,*,k+3]/(data[*,*,*,k+2]*data[*,*,*,k+1]),/Phase)
                (*project.dataarray[CI].Bez)[*,*,*,i,1] = atan(data[*,*,*,k]*data[*,*,*,k+2]/(data[*,*,*,k+3]*data[*,*,*,k+1]),/Phase)
                (*project.dataarray[CI].Bez)[*,*,*,i,2] = atan(data[*,*,*,k]*data[*,*,*,k+1]/(data[*,*,*,k+2]*data[*,*,*,k+3]),/Phase)                 
                (*project.dataarray[CI].Bez)[*,*,*,i,3] = mean(abs(data[*,*,*,k:k+3]),dimension = 4)
                end
                
      else    : begin
                (*project.dataarray[CI].Bez)[*,*,*,i,0] = atan(data[*,*,*,k+1]/data[*,*,*,k],/Phase)
                (*project.dataarray[CI].Bez)[*,*,*,i,1] = mean(abs(data[*,*,*,k:k+1]),dimension = 4)
                end            
  endcase  
  k+=nskip
endfor

if n_elements(TE) gt 1 then begin
  progressbar = Obj_New('progressbar', Color='red', Text='Calculating phase derivative', /NOCANCEL)
  progressbar -> Start
  counter = 0UL
  ; First derivative matrix (forward difference)
  D = identity(n_elements(TE)) - shift(identity(n_elements(TE)),-n_elements(TE))
  D[0,n_elements(TE)-1] = 0
  D[n_elements(TE)-2,n_elements(TE)-1] = -1
  D/=(TE[0]-TE[1])
  
  for i = 0,N[0]-1 do begin
    for j = 0,N[1]-1 do begin
      for k = 0,N[2]-1 do begin 
          for l = 0,Ndelays-1 do begin
            for m = 0,2 do begin            
              echo_data = reform((*project.dataarray[CI].Bez)[i,j,k,l:Nt-1:Ndelays,m],n_elements(TE))
              echo_data_uw = phunwrap(echo_data)
              (*project.dataarray[CI].Bez)[i,j,k,l:Nt-1:Ndelays,m] = D##echo_data_uw/float(gamma*xi)  
              counter++   
              progressBar -> Update, 100.0*float(counter)/float(N[0]*N[1]*N[2]*Ndelays*3)        
            endfor           
          endfor
      endfor
    endfor
  endfor
  progressbar -> Destroy
endif else begin
  if scheme eq 'Hadmard' or scheme eq 'Four Point' or scheme eq 'Six Point' then (*project.dataarray[CI].Bez)[*,*,*,*,0:2]/=float(TE*gamma*xi) $
    else (*project.dataarray[CI].Bez)[*,*,*,*,0]/=float(TE*gamma*xi)
endelse

if project.procPramArray[CI].Bez_up_stat eq 1 then unwrap_phase_Bez,project.dataarray[CI].Bez, scheme

Time = []
for i = 0,n_elements(TE)-1 do Time = [Time,TE[i]+*project.imndarray[CI].tgrdelay]
*project.dataarray[CI].Bez = (*project.dataarray[CI].Bez)[*,*,*,sort(Time),*]

orient_Bez_data, project.dataarray[CI].Bez

project.procPramArray[CI].Bez_proccess_flag = 1

end

; Subroutine name: Modelfit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Using MPFIT to perform the exponential fitting of eddy currents
;
; Editing Information:

pro Modelfit, n, T, Bez, BW, params, pvalue, r2r, r2c, model_time

common scan_data
CI = project.ci

maxiter     = 50                                                     ; Maximum number of random exploratory steps
rand_matrix = indgen(maxiter)
ind         = rand_matrix[sort(randomu(!NULL,maxiter))]
Nt        = project.procPramArray[CI].Bez_Nfit
alpha_max = project.procPramArray[project.ci].Bez_alpha_max/100.0   ; Maximum eddy current amplitude
shape   = project.imndarray[project.ci].tgrshape                     ; Shape of test gradient pulse
case shape of
  's'  : if model_time eq 0 then model_fn = 'Bez_fit_time_sine_eddy' else model_fn = 'Bez_fit_time_sine_oscillation'         
  'c'  : if model_time eq 0 then model_fn = 'Bez_fit_time_scosine_eddy' else model_fn = 'Bez_fit_time_scosine_oscillation'
  'g'  : if model_time eq 0 then model_fn = 'Bez_fit_time_Gauss_eddy' else model_fn = 'Bez_fit_time_Gauss_oscillation'
  'f'  : if model_time eq 0 then model_fn = 'Bez_fit_time_FresnelT_eddy' else model_fn = 'Bez_fit_time_FresnelT_oscillation'
  'i'  : if model_time eq 0 then model_fn = 'Bez_fit_time_SincT_eddy' else model_fn = 'Bez_fit_time_SincT_oscillation'
  else : if model_time eq 0 then model_fn = 'Bez_fit_time_trapz_eddy' else model_fn = 'Bez_fit_time_trapz_oscillation'
endcase

if model_time eq 0 then pinfo = replicate({value:0.0,limited:[0,0],limits:[0.0,0.0]},2*n) else pinfo = replicate({value:0.0,limited:[0,0],limits:[0.0,0.0]},4*n)
Np = n_elements(pinfo)

; Initial estimates for time constants
tau_est  = -Bez/deriv(T[0:Nt-1],Bez)              ; Estimated time constants based on simple algebra of the model
tau_min = min(abs(tau_est),/nan) 
tau_max = T[Nt-1]
tau     = make_array(maxiter,/double)
tau     = tau_min + (tau_max-tau_min)*findgen(maxiter)/(maxiter-1)

; Constraint on time constant based on excitation bandwidth and sampling time
pinfo[1:Np-1:Np/n].limited[0] = 1
pinfo[1:Np-1:Np/n].limits[0]  = 1/float(2.0*!Pi*BW)
pinfo[1:Np-1:Np/n].limited[1] = 1
pinfo[1:Np-1:Np/n].limits[1]  = 2*T[Nt-1]
pinfo[1:Np-1:Np/n].value = tau[ind[0:n-1]]

; Constraint on amplitudes based on widget input
pinfo[0:Np-1:Np/n].limited[0] = 1
pinfo[0:Np-1:Np/n].limited[1] = 1
pinfo[0:Np-1:Np/n].value = (2*randomu(!NULL,n)-1)*alpha_max
pinfo[0:Np-1:Np/n].limits[0] = -alpha_max
pinfo[0:Np-1:Np/n].limits[1]  = alpha_max

; Constraint on frequency and phase for oscillatory models based on Nyquist sampling criterion and +/- Pi respectively
if model_time eq 1 then begin
  pinfo[2:Np-1:Np/n].limited[0] = 1
  pinfo[2:Np-1:Np/n].limited[1] = 1  
  pinfo[2:Np-1:Np/n].limits[0] = -1/(2.0*min(deriv(T[0:Nt-1])))    ; Nyquist criterion
  pinfo[2:Np-1:Np/n].limits[1] = 1/(2.0*min(deriv(T[0:Nt-1])))    ; Nyquist criterion  
  pinfo[2:Np-1:Np/n].value = (2*randomu(!NULL,n)-1)/(2.0*min(deriv(T[0:Nt-1]))) 
  
  pinfo[3:Np-1:Np/n].limited[0] = 1
  pinfo[3:Np-1:Np/n].limited[1] = 1
  pinfo[3:Np-1:Np/n].limits[0] = -!Pi
  pinfo[3:Np-1:Np/n].limits[1]  = !Pi
  pinfo[3:Np-1:Np/n].value = (2*randomu(!NULL,n)-1)*!Pi
endif

fargs   = {x:T[0:Nt-1],y:Bez,err:1}
params  = MPFIT(model_fn, functargs = fargs, parinfo=pinfo,bestnorm = chisquared,status=status)
SS_tot  = (Nt-1)*variance(Bez,/double)
r2c     = 1 - chisquared/SS_tot
iter    = 0
while r2c lt 0.8 do begin
  iter++
  if iter gt maxiter then break
  pinfo[0:Np-1:Np/n].value = (2*randomu(!NULL,n)-1)*alpha_max  
  pinfo[1:Np-1:Np/n].value = tau[ind[0:n-1]] 
  if model_time eq 1 then begin
    pinfo[2:Np-1:Np/n].value =  (2*randomu(!NULL,n)-1)/(2.0*min(deriv(T[0:Nt-1])))
    pinfo[3:Np-1:Np/n].value = (2*randomu(!NULL,n)-1)*!Pi
  endif
  paramsfit  = MPFIT(model_fn, functargs = fargs, parinfo=pinfo,bestnorm = chisquared,status=status)
  if status le 0 then continue else begin
    SS_tot  = (Nt-1)*variance(Bez,/double)
    r2cfit  = 1 - chisquared/SS_tot 
    if r2cfit gt r2c then begin
      params = paramsfit
      r2c = r2cfit     
    endif 
  endelse   
endwhile

if model_time eq 0 then begin
  Fobs   = ((r2c-r2r)/2.0)/((1-r2c)/float(Nt-2*n-1))
  pvalue = 1-f_pdf(Fobs,2.0,Nt-2*n-1)  
endif else if model_time eq 1 then begin
  Fobs   = ((r2c-r2r)/4.0)/((1-r2c)/float(Nt-4*n-1))
  pvalue = 1-f_pdf(Fobs,4.0,Nt-4*n-1)
endif

end

; Subroutine name: SpatialFit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To perform linear fit of the eddy current magnetic field in space
;
; Editing Information:

pro SpatialFit, CI, B0, G, w = w, ax = ax, flag, model_space = model_space, it = it

common scan_data
common Bez_params

Lr        = 10*project.imndarray[CI].f_fov             ; Readout FOV (mm)
Lp1       = 10*project.imndarray[CI].p_fov             ; PE1 FOV (mm)
Lp2       = 10*project.imndarray[CI].s_fov             ; PE2 FOV (mm)
Loff_ro   = 10*project.imndarray[CI].f_fov_offset      ; Readout FOV offset (mm)
Loff_p1   = 10*project.imndarray[CI].p_fov_offset      ; PE FOV offset (mm)
Loff_p2   = 10*project.imndarray[CI].s_fov_offset      ; PE2 FOV offset (mm)
Nr        = project.imndarray[CI].fdim
Np1       = project.imndarray[CI].pdim
Np2       = project.imndarray[CI].sdim
Nt        = project.procPramArray[CI].Bez_Nfit

; Zeropad effect
if project.procPramArray[CI].zpad_flag eq 1 then begin
  Nr*=2           
  Np1*=2                 
  if project.imndArray[CI].dimensions eq 3 then Np2*=2              
endif 
orient       = project.imndarray[CI].orientation          ; Slice orientation
Ithreshold   = project.procPramArray[CI].Bez_Ithreshold   ; Image threshold  

if ~ptr_valid(project.dataarray[CI].Bez) then begin
  mbox = dialog_message(['Please calculate the field map for scan #'+strtrim(CI+1,2)],/error,/center)
  flag = 0
  return
endif else flag = 1

if n_elements(ax) le 0 then ax = project.procPramArray[CI].Bez_axis
Bez = (*project.dataarray[CI].Bez)[*,*,*,*,ax]       ; Calculated field map
; Magnitude image
if scheme eq 'Hadamard' or scheme eq 'Four Point' or scheme eq 'Six Point' then Imag = reform((*project.dataarray[CI].Bez)[*,*,*,*,3]) $
else  Imag = reform((*project.dataarray[CI].Bez)[*,*,*,*,1])        

Ntotal = Nr*Np1*Np2
L      = [[Lr],[Lp1],[Lp2]]
Loff   = [[Loff_ro],[Loff_p1],[Loff_p2]]
r      = findgen(Nr)/float(Nr-1)-0.5
p1     = findgen(Np1)/float(Np1-1)-0.5
p2     = findgen(Np2)/float(Np2-1)-0.5
meshgrid,r,p1, z = p2
 
; Transforming from image (r,p1,p2) to gradient coordinates (x,y,z)
trans_mat = log2physgradtransform()
L = trans_mat##L
Loff = trans_mat##Loff
pos = [[reform(r,Ntotal,1)],[reform(p1,Ntotal,1)],[reform(p2,Ntotal,1)]]
pos_trans = trans_mat##pos                   
x = L[0]*pos_trans[*,0]-Loff[0]
y = L[1]*pos_trans[*,1]-Loff[1]
z = L[2]*pos_trans[*,2]-Loff[2]

progressbar = Obj_New('progressbar', Color='red', Text='Spatial Fitting', /NOCANCEL)
progressbar -> Start
for i = 0,Nt-1 do begin
    if n_elements(it) gt 0 then if i ne it then continue
    Imask = Imag[*,*,*,i] ge max(Imag[*,*,*,i])*float(Ithreshold)/100.0
    Bez[*,*,*,i]*=Imask
    Bez_sample = reform(Bez[*,*,*,i],Ntotal)
    ind = where(Imask)
    xmask = x[ind]
    ymask = y[ind]
    zmask = z[ind]
    Bezmask = Bez_sample[ind]
    XYZ = [[xmask],[ymask],[zmask]]
    Ns = n_elements(ind)
        
    pinfo = replicate({value:1.0},4)
    fargs = {x:XYZ,y:Bezmask,err:1}
    SS_tot = (Ns-1)*variance(Bezmask)
      
    r2_0 = 1-total(Bezmask^2)/SS_tot
    
    params_c = MPFIT('Bez_fit_space_C', functargs = fargs, parinfo=pinfo[0],bestnorm = chisquared,status=status)
    if status gt 0 then r2_c = 1 - chisquared/SS_tot else r2_c = !values.D_INFINITY 
    
    params_x = MPFIT('Bez_fit_space_X', functargs = fargs, parinfo=pinfo[1],bestnorm = chisquared,status=status)
    if status gt 0 then r2_x = 1 - chisquared/SS_tot else r2_x = !values.D_INFINITY
    
    params_y = MPFIT('Bez_fit_space_Y', functargs = fargs, parinfo=pinfo[2],bestnorm = chisquared,status=status)
    if status gt 0 then r2_y = 1 - chisquared/SS_tot else r2_y = !values.D_INFINITY
    
    params_z = MPFIT('Bez_fit_space_Z', functargs = fargs, parinfo=pinfo[3],bestnorm = chisquared,status=status)
    if status gt 0 then r2_z = 1 - chisquared/SS_tot else r2_z = !values.D_INFINITY
    
    params_cx = MPFIT('Bez_fit_space_CX', functargs = fargs, parinfo=[pinfo[0],pinfo[1]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_cx = 1 - chisquared/SS_tot else r2_cx = !values.D_INFINITY
    
    params_cy = MPFIT('Bez_fit_space_CY', functargs = fargs, parinfo=[pinfo[0],pinfo[2]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_cy = 1 - chisquared/SS_tot else r2_cy = !values.D_INFINITY
    
    params_cz = MPFIT('Bez_fit_space_CZ', functargs = fargs, parinfo=[pinfo[0],pinfo[3]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_cz = 1 - chisquared/SS_tot else r2_cz = !values.D_INFINITY
    
    params_xy = MPFIT('Bez_fit_space_XY', functargs = fargs, parinfo=[pinfo[1],pinfo[2]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_xy = 1 - chisquared/SS_tot else r2_xy = !values.D_INFINITY

    params_xz = MPFIT('Bez_fit_space_XZ', functargs = fargs, parinfo=[pinfo[1],pinfo[3]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_xz = 1 - chisquared/SS_tot else r2_xz = !values.D_INFINITY

    params_yz = MPFIT('Bez_fit_space_YZ', functargs = fargs, parinfo=[pinfo[2],pinfo[3]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_yz = 1 - chisquared/SS_tot else r2_yz = !values.D_INFINITY
    
    params_cxy = MPFIT('Bez_fit_space_CXY', functargs = fargs, parinfo=[pinfo[0],pinfo[1],pinfo[2]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_cxy = 1 - chisquared/SS_tot else r2_cxy = !values.D_INFINITY
    
    params_cxz = MPFIT('Bez_fit_space_CXZ', functargs = fargs, parinfo=[pinfo[0],pinfo[1],pinfo[3]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_cxz = 1 - chisquared/SS_tot else r2_cxz = !values.D_INFINITY
    
    params_cyz = MPFIT('Bez_fit_space_CYZ', functargs = fargs, parinfo=[pinfo[0],pinfo[2],pinfo[3]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_cyz = 1 - chisquared/SS_tot else r2_cyz = !values.D_INFINITY
    
    params_xyz = MPFIT('Bez_fit_space_XYZ', functargs = fargs, parinfo=[pinfo[1],pinfo[2],pinfo[3]],bestnorm = chisquared,status=status)
    if status gt 0 then r2_xyz = 1 - chisquared/SS_tot else r2_xyz = !values.D_INFINITY
    
    params_cxyz = MPFIT('Bez_fit_space_CXYZ', functargs = fargs, parinfo=pinfo,bestnorm = chisquared,status=status)
    if status gt 0 then r2_cxyz = 1 - chisquared/SS_tot else r2_cxyz = !values.D_INFINITY
    
   rsquared = [r2_0, r2_c, r2_x, r2_y, r2_z, r2_cx, r2_cy, r2_cz, r2_xy, r2_xz, r2_yz, r2_cxy, r2_cxz, r2_cyz, r2_xyz, r2_cxyz]
   model = where(rsquared eq max(rsquared) and rsquared gt 0)
   if model eq -1 then model = 0 
   if n_elements(w) gt 0 then w[i] = rsquared[model]
   if n_elements(model_space) gt 0 then model_space[i] = model
    
    case model of
      0 : begin
          B0[i] = 0.0
          G[i,*] = dblarr(3)   
          end
      
      1 : begin
          B0[i] = params_c
          G[i,*] = dblarr(3)
          end       
     
      2 : begin
          B0[i] = 0
          G[i,*] = [params_x,0.0,0.0]
          end             

      3 : begin
          B0[i] = 0
          G[i,*] = [0.0,params_y,0.0]
          end 
             
      4 : begin
          B0[i] = 0
          G[i,*] = [0.0,0.0,params_z]
          end        
      
      5 : begin
          B0[i] = params_cx[0]
          G[i,*] = [params_cx[1],0.0,0.0]
          end     
      
      6 : begin
          B0[i] = params_cy[0]
          G[i,*] = [0.0,params_cy[1],0.0]
          end    
      
      7 : begin
          B0[i] = params_cz[0]
          G[i,*] = [0.0,0.0,params_cz[1]]
          end                      
      
      8 : begin
          B0[i] = 0.0
          G[i,*] = [params_xy[0],params_xy[1],0.0]
          end
          
       9 : begin
           B0[i] = 0.0
           G[i,*] = [params_xz[0],0.0,params_xz[1]]
           end                    

       10 : begin
            B0[i] = 0.0
            G[i,*] = [0.0,params_yz[0],params_yz[1]]
            end
       
       11 : begin
            B0[i] = params_cxy[0]
            G[i,*] = [params_cxy[1],params_cxy[2],0.0]
            end
        
        12 : begin
             B0[i] = params_cxz[0]
             G[i,*] = [params_cxz[1],0.0,params_cxz[2]]
             end

        13 : begin
             B0[i] = params_cyz[0]
             G[i,*] = [0.0,params_cyz[1],params_cyz[2]]
             end   
        
        14 : begin
             B0[i] = 0.0
             G[i,*] = [params_xyz[0],params_xyz[1],params_xyz[2]]
             end    
             
        15 : begin
             B0[i] = params_cxyz[0]
             G[i,*] = [params_cxyz[1],params_cxyz[2],params_cxyz[3]]
             end                                         
    endcase    
    
    progressBar -> Update, 100.0*float(i+1)/float(Nt) 
endfor

progressbar -> Destroy
end

; Subroutine name: TemporalFit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To fit the eddy current magnetic field or mechanical oscillations over time using multi-exponential function
;
; Editing Information:

pro TemporalFit, data, Ts, params, model_time, r2_time

common Bez_params

case model_time of
  0 : model_txt = 'Eddy Current'
  1 : model_txt = 'Oscillations'
endcase
SLE  = 0.2     ; Significance level

progressbar = Obj_New('progressbar', Color='red', Text='Temporal Fitting - ' + model_txt, /NOCANCEL)
progressbar -> Start
for i=0,3 do begin  
      r2r = 1-total(data[i,*]^2)/((n_elements(Ts)-1)*variance(data[i,*],/double))
      if model_time eq 0 then begin
          params0 = [0.0,!values.d_infinity]
          Fobs   = (r2r/2.0)/((1-r2r)/float(n_elements(Ts)-3))
          pvalue0 = 1-f_pdf(Fobs,2.0,n_elements(Ts)-3.0)
       endif else if model_time eq 1 then begin
          params0 = [0.0,!values.d_infinity,0.0,0.0]
          Fobs   = (r2r/4.0)/((1-r2r)/float(n_elements(Ts)-5))
          pvalue0 = 1-f_pdf(Fobs,4.0,n_elements(Ts)-5)
       endif 

      Modelfit, 1, Ts, reform(data[i,*]), BW, params1, pvalue1, r2r, r2c, model_time
      n = 1
      while pvalue1 lt SLE and pvalue1 gt 0.05  do begin
        n++
        if model_time eq 0 and n ge n_elements(Ts)/2 then break
        if model_time eq 1 and n ge n_elements(Ts)/4 then break
        pvalue0 = pvalue1
        params0 = params1
        r2r = r2c
        Modelfit, n, Ts, reform(data[i,*]), BW, params1, pvalue1, r2r, r2c, model_time
      endwhile
      if pvalue1 ge pvalue0 then begin
        *params[i] = params0 
        r2_time[i] = r2r
      endif else begin
        *params[i] = params1
        r2_time[i] = r2c
      endelse
      progressBar -> Update, (float(i+1)/(float(4)))*100.0    
endfor      
progressbar -> Destroy
end

; Subroutine name: Spectrum
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To fit the eddy current magnetic field or mechanical oscillations over time and space
;
; Editing Information:

pro Spectrum

  common scan_data
  common Bez_params

  CI = project.ci
  Nt = project.procPramArray[CI].Bez_Nfit
  B0 = make_array(Nt)
  G  = make_array(Nt,3)
  SpatialFit, CI, B0, G, ax = ax
  if project.procPramArray[CI].Bez_smooth eq 1 then begin
    B0 = Smoothing(B0,2)
    for i = 0,2 do G[*,i] = Smoothing(G[*,i],2)
  endif
  G*=1e3
  
  if project.procPramArray[CI].Bez_apodize eq 1 then begin
    phase = !Pi/2
    B0*= sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2
    G*= rebin(sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2,Nt,3)  
  endif
    
  if project.procPramArray[CI].Bez_zf_stat eq 1 then begin
    B0 = [B0, make_array(Nt,/double)]
    G = [G, make_array(Nt,3,/double)]
    Fs = (findgen(2*Nt)/(2*Nt-1)-0.5)/(T[1]-T[0])
    spectrum = complexarr(2*Nt,4)
  endif else begin
    Fs = (findgen(Nt)/(Nt-1)-0.5)/(T[1]-T[0])
    spectrum = complexarr(Nt,4)
  endelse
  
  case project.procPramArray[CI].Bez_axis of
    0 : ax = 'X'
    1 : ax = 'Y'
    2 : ax = 'Z'
    else : ax = axis
  endcase
  
  spectrum[*,0] = fft(B0,/center)
  spectrum[*,1] = fft(G[*,0],/center)
  spectrum[*,2] = fft(G[*,1],/center)
  spectrum[*,3] = fft(G[*,2],/center)
  h1 = plot(Fs,abs(spectrum[*,0]),layout = [2,1,1],"b2",window_title = ax+'-gradient spectrum', title = '$B_0$')
  h2 = plot(Fs,abs(spectrum[*,1]),"r2",layout = [2,1,2],/current,name = '$G_x$',title = 'Gradient')
  h3 = plot(Fs,abs(spectrum[*,2]),"g2",layout = [2,1,2],/current,/overplot, name = '$G_y$')
  h4 = plot(Fs,abs(spectrum[*,3]),"k2",layout = [2,1,2],/current,/overplot,name = '$G_z$')
  
  leg = LEGEND(TARGET=[h2,h3,h4],/DATA, /AUTO_TEXT_COLOR, FONT_SIZE = 15)
    
  ax1 = h1.axes
  ax1[0].TITLE = 'Frequency (KHz)'
  h1.FONT_SIZE = 20
  
  ax2 = h2.axes
  ax2[0].TITLE = 'Frequency (KHz)'
  h2.FONT_SIZE = 20
  
end  




; Subroutine name: TimeSeries
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To fit the eddy current magnetic field or mechanical oscillations over time and space
;
; Editing Information:

pro TimeSeries

  common scan_data
  common Bez_params

  CI = project.ci
  Nt = n_elements(T)
  B0 = make_array(Nt)
  G  = make_array(Nt,3)
  Nfit = project.procPramArray[CI].Bez_Nfit
    
  SpatialFit, CI, B0, G, flag
  if project.procPramArray[CI].Bez_smooth eq 1 then begin
    B0 = Smoothing(B0,2)
    for i = 0,2 do G[*,i] = Smoothing(G[*,i],2)
  endif
  
  if project.procPramArray[CI].Bez_apodize eq 1 then begin
    phase = !Pi/2.0
    B0*= sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2
    G*= rebin(sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2,Nt,3)
  endif
  
  case project.procPramArray[CI].Bez_axis of
    0 : ax = 'X'
    1 : ax = 'Y'
    2 : ax = 'Z'
    else : ax = axis
  endcase
  h1 = plot(T[0:Nfit-1],B0,"b2+",layout = [2,2,1],window_title = ax+' gradient eddy current amplitudes over time', margin = 0.2,/xlog)
  h2 = plot(T[0:Nfit-1],1e3*G[*,0],"b2+",/current,layout = [2,2,2], margin = 0.2,/xlog)
  h3 = plot(T[0:Nfit-1],1e3*G[*,1],"b2+",/current,layout = [2,2,3], margin = 0.2,/xlog)
  h4 = plot(T[0:Nfit-1],1e3*G[*,2],"b2+",/current,layout = [2,2,4], margin = 0.2,/xlog)

  ax1 = h1.axes
  ax1[0].TITLE = 'Time (ms)'
  ax1[1].TITLE = '$B_0$ ($\mu$T) '
  h1.FONT_SIZE = 20

  ax2 = h2.axes
  ax2[0].TITLE = 'Time (ms)'
  ax2[1].TITLE = '$G_x$ ($\mu$T/m)'
  h2.FONT_SIZE = 20

  ax3 = h3.axes
  ax3[0].TITLE = 'Time (ms)'
  ax3[1].TITLE = '$G_y$ ($\mu$T/m)'
  h3.FONT_SIZE = 20

  ax4 = h4.axes
  ax4[0].TITLE = 'Time (ms)'
  ax4[1].TITLE = '$G_z$ ($\mu$T/m)'
  h4.FONT_SIZE = 20
    
end 

; Subroutine name: export_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Dummy event handler for flow animation (required for the movie to loop indefinitely)
;
; Editing Information:

pro export_Bez

  common scan_data
  common Bez_params

  CI = project.ci
  Nt = n_elements(T)
  B0 = make_array(Nt)
  G  = make_array(Nt,3)
  Nfit = project.procPramArray[CI].Bez_Nfit
    
  SpatialFit, CI, B0, G, flag
  if project.procPramArray[CI].Bez_smooth eq 1 then begin
    B0 = Smoothing(B0,2)
    for i = 0,2 do G[*,i] = Smoothing(G[*,i],2)
  endif
  
  if project.procPramArray[CI].Bez_apodize eq 1 then begin
    phase = !Pi/2.0
    B0*= sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2
    G*= rebin(sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2,Nt,3)
  endif
  
  case project.procPramArray[CI].Bez_axis of
    0 : ax = 'X'
    1 : ax = 'Y'
    2 : ax = 'Z'
    else : ax = axis
  endcase
  
  data = [[B0],[1e3*G]]
  fpath = DIALOG_PICKFILE(Title = 'Select Directory')
  OpenW, lun,fpath + '.raw', /get_lun
  WriteU,lun,data
  close, lun
  free_lun, lun

  topbase = widget_base(title = 'File export information' ,tlb_frame_attr=1,/grid_layout, frame = 1)
  mainbase = widget_base(topbase)
  rowtxt1 = 'B0 and gradient field from' + ax +'-gradient exported as a binary file to ' + fpath + '.raw'
  rowtxt2 = 'Data type       : ' + strtrim(typename(data),2)
  rowtxt3 = 'Data order      : B0 (uT) Gx (uT/m) Gy (uT/m) Gz (uT/m)'
  rowtxt4 = 'File dimensions :' + strcompress(strjoin(size(data,/DIMENSIONS)))
  txt = [[rowtxt1],[rowtxt2],[rowtxt3],[rowtxt4]]
  rowlabel = widget_text(mainbase, value = txt, ysize = 5)
  widget_control, topbase, /realize

end 
; Subroutine name: Fitting
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To fit the eddy current magnetic field or mechanical oscillations over time and space 
;
; Editing Information:

pro Fitting, B0, G, w, pre_emph, oscns = oscns, r2_time, model_time, model_space
  
  common scan_data
  common Bez_params
 
  CI = project.ci  
  Nt = project.procPramArray[CI].Bez_Nfit
  SpatialFit, CI, B0, G, w = w, flag, model_space = model_space
  if project.procPramArray[CI].Bez_smooth eq 1 then begin
    B0 = Smoothing(B0,2)
    for i = 0,2 do G[*,i] = Smoothing(G[*,i],2)
  endif

  if project.procPramArray[CI].Bez_apodize eq 1 then begin
    phase = !Pi/2.0
    B0*= sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2
    G*= rebin(sin((!Pi-phase)*T[0:Nt-1]/T[Nt-1]+phase)^2,Nt,3)  
  endif
  
  data = [transpose(B0),transpose(G)]  
  if model_time eq 2 then begin
      TemporalFit, data, T[0:Nt-1], pre_emph, 0, r2_time
        case project.imndarray[project.ci].tgrshape of
          's' : for i = 0,3 do data[i,*] += Bez_fit_time_sine_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1)                  
               
          'c' : for i = 0,3 do data[i,*] += Bez_fit_time_scosine_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1)  
                
          'g' : for i = 0,3 do data[i,*] += Bez_fit_time_Gauss_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1)

          'f' : for i = 0,3 do data[i,*] += Bez_fit_time_FresnelT_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1)

          'i' : for i = 0,3 do data[i,*] += Bez_fit_time_SincT_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1)

          else : for i = 0,3 do data[i,*] += Bez_fit_time_trapz_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1)
          
        endcase
        TemporalFit, data, T[0:Nt-1], oscns, 1, r2_time     
    endif else TemporalFit, data, T[0:Nt-1], pre_emph, model_time, r2_time
end

; Subroutine name: Plot_Bez_fit
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To plot the fitted eddy current slope and intercept over time
;
; Editing Information:

pro Plot_Bez_fit, model_time
  common scan_data
  common Bez_params
  
  CI = project.ci
  gcoil = project.imndArray[CI].gcoil
  GradCoilParam = gradient_specs(gcoil)
  Imax_B0 = GradCoilParam.Imax_B0
  C0 = GradCoilParam.C0
  Gmax = GradCoilParam.Gmax
  
  Nt = project.procPramArray[CI].Bez_Nfit
  B0 = make_array(Nt)
  G  = make_array(Nt,3)
  w = replicate(0.0,Nt)
  model_space = replicate(0,Nt)
  r2_time = replicate(0.0,4)
  if ptr_valid(pre_emph) then ptr_free, pre_emph else pre_emph = ptrarr(4,/ALLOCATE_HEAP)
  if model_time eq 2 then begin
    if ptr_valid(oscns) then ptr_free, oscns else oscns = ptrarr(4,/ALLOCATE_HEAP)
    Fitting, B0, G, w, pre_emph, oscns = oscns, r2_time, model_time, model_space
  endif else Fitting, B0, G, w, pre_emph, r2_time, model_time, model_space

  result = [transpose(B0),transpose(G)]
  
  tabch = string(09B)
  text_output_space = strarr(Nt+3)
  text_output_space[0] = 'Linear eddy current slope and intercept for '+project.imndarray[ci].display_name
  text_output_space[1] = string(project.procpramarray[ci].freq_interp, project.procpramarray[ci].phase_interp, project.procpramarray[ci].slice_interp,format='(%"Interpolation (RPS): %d,%d,%d")')
  text_output_space[2] = 'Time (ms)' + tabch + 'B0(uT)' + tabch + 'Gx(uT/m)' + tabch + 'Gy(uT/m)'  + tabch + 'Gz(uT/m)' + tabch + 'Model' + tabch + 'Rsquared' 
  for i = 0, Nt-1 do text_output_space[i+3] = strtrim(string(T[i],format = '(F8.2)'),2) + tabch + strtrim(string(B0[i],format = '(F15.8)'),2) + tabch + $
  strtrim(string(1e3*G[i,0],format = '(F15.8)'),2) + tabch +  strtrim(string(1e3*G[i,1],format = '(F15.8)'),2) + tabch +  strtrim(string(1e3*G[i,2],format = '(F15.8)'),2) $
  + tabch +  strtrim(string(model_space[i],format = '(I2)'),2) + tabch +  strtrim(string(w[i],format = '(F8.2)'),2)    
  
  Nc = 500*(T[Nt-1]-T[0])/T[0]
  Tc = T[0] + (T[Nt-1]-T[0])*findgen(Nc)/(Nc-1)
  fit_result = make_array(4,n_elements(Tc),/double)
  nmax = max([n_elements(*pre_emph[0]),n_elements(*pre_emph[1]),n_elements(*pre_emph[2]),n_elements(*pre_emph[3])])
  if model_time eq 1 then nmax/=4.0 else nmax/=2.0
  text_output_time = strarr(8)
  text_output_time[0] = 'Pre-emphasis constants for '+project.imndarray[ci].display_name 
  text_output_time[1] = 'Test gradient - Shape :'+ shape +',Strength : '+ strtrim(gBez,2) +' mT/m, Duration : ' + strtrim(tBez,2) + ' ms'
  text_output_time[2] = string(project.procpramarray[ci].freq_interp, project.procpramarray[ci].phase_interp, project.procpramarray[ci].slice_interp,format='(%"Interpolation (RPS): %d,%d,%d")')
  text_output_time[3] = 'Distortion' + tabch
  for j=1,nmax do begin
    text_output_time[3] += 'alpha' +strtrim(j,2) +'(%)' + tabch + 'tau'+strtrim(j,2)+'(ms)' + tabch
    if model_time eq 1 then text_output_time[3] += 'f' +strtrim(j,2) +'(KHz)' + tabch + 'phi'+strtrim(j,2)+'(radians)' + tabch
  endfor
  if model_time eq 2 then begin
    nosmax = max([n_elements(*oscns[0]),n_elements(*oscns[1]),n_elements(*oscns[2]),n_elements(*oscns[3])])/4
    for j=1,nosmax do text_output_time[3] += 'beta' +strtrim(j,2) +'(%)' + tabch + 'tau'+strtrim(j,2)+'(ms)' + tabch $
           + 'f' +strtrim(j,2) +'(KHz)' + tabch + 'phi'+strtrim(j,2)+'(radians)' + tabch
  endif
  text_output_time[3] += 'Rsquared'
  for i = 0,3 do begin
    case i of
    0 : text_output_time[i+4] += 'B0'+ tabch
    1 : text_output_time[i+4] += 'Gx'+ tabch
    2 : text_output_time[i+4] += 'Gy'+ tabch
    3 : text_output_time[i+4] += 'Gz'+ tabch
    endcase
       
    case project.imndarray[CI].tgrshape of
    's' : begin
          case model_time of
          0 : fit_result[i,*] = -Bez_fit_time_sine_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) 
          1 : fit_result[i,*] = -Bez_fit_time_sine_oscillation(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
          2 : begin
              fp = -Bez_fit_time_sine_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0) $
                -Bez_fit_time_sine_oscillation(*oscns[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0)
              r2_time[i] = 1-total((result[i,*]-fp)^2)/((Nt-1)*variance(result[i,*],/double))
              fit_result[i,*] = -Bez_fit_time_sine_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) $
                              - Bez_fit_time_sine_oscillation(*oscns[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
               end               
          endcase
          end
    'c' : begin
          case model_time of
           0 : fit_result[i,*] = -Bez_fit_time_scosine_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) 
           1 : fit_result[i,*] = -Bez_fit_time_scosine_oscillation(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
           2 : begin
               fp = -Bez_fit_time_scosine_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0) $
                -Bez_fit_time_scosine_eddy(*oscns[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0)
               r2_time[i] = 1-total((result[i,*]-fp)^2)/((Nt-1)*variance(result[i,*],/double))
               fit_result[i,*] = -Bez_fit_time_scosine_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) $
                               - Bez_fit_time_scosine_oscillation(*oscns[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
               end                
          endcase   
          end                   
    'g' : begin
          case model_time of 
            0 : fit_result[i,*] = -Bez_fit_time_Gauss_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) 
            1 : fit_result[i,*] = -Bez_fit_time_Gauss_oscillation(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
            2 : begin
                fp = -Bez_fit_time_Gauss_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0) $
                  -Bez_fit_time_Gauss_oscillation(*oscns[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0)
                r2_time[i] = 1-total((result[i,*]-fp)^2)/((Nt-1)*variance(result[i,*],/double))
                fit_result[i,*] = -Bez_fit_time_Gauss_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) $
                                - Bez_fit_time_Gauss_oscillation(*oscns[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
                end                
           endcase 
           end
    'f' : begin
          case model_time of 
            0 : fit_result[i,*] = -Bez_fit_time_FresnelT_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) 
            1 : fit_result[i,*] = -Bez_fit_time_FresnelT_oscillation(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
            2 : begin
                fp = -Bez_fit_time_FresnelT_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0) $
                  -Bez_fit_time_FresnelT_oscillation(*oscns[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0)
                r2_time[i] = 1-total((result[i,*]-fp)^2)/((Nt-1)*variance(result[i,*],/double))
                fit_result[i,*] = -Bez_fit_time_FresnelT_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) $
                                -Bez_fit_time_FresnelT_oscillation(*oscns[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
                end                 
           endcase
           end      
    'i' : begin
          case model_time of 
            0 : fit_result[i,*] = -Bez_fit_time_SincT_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) 
            1 : fit_result[i,*] = -Bez_fit_time_SincT_oscillation(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
            2 : begin
                fp = -Bez_fit_time_SincT_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0) $
                  -Bez_fit_time_SincT_oscillation(*oscns[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0)
                r2_time[i] = 1-total((result[i,*]-fp)^2)/((Nt-1)*variance(result[i,*],/double))
                fit_result[i,*] = -Bez_fit_time_SincT_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) $
                                -Bez_fit_time_SincT_oscillation(*oscns[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
                end
          endcase
          end                      
     else : begin
            case model_time of 
            0 : fit_result[i,*] = -Bez_fit_time_trapz_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) 
            1 : fit_result[i,*] = -Bez_fit_time_trapz_oscillation(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
            2 : begin
                fp = -Bez_fit_time_trapz_eddy(*pre_emph[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0) $
                  -Bez_fit_time_trapz_oscillation(*oscns[i], X = T[0:Nt-1], y = dblarr(Nt), err = 1.0)
                r2_time[i] = 1-total((result[i,*]-fp)^2)/((Nt-1)*variance(result[i,*],/double))
                fit_result[i,*] = -Bez_fit_time_trapz_eddy(*pre_emph[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0) $
                                -Bez_fit_time_trapz_oscillation(*oscns[i], X = Tc, y = dblarr(n_elements(Tc)), err = 1.0)
                end                
            endcase
            end                        
     endcase
     fit_result[i,where(finite(fit_result[i,*],/NAN),/NULL)] = 0.0     

     ; B0 scaling
     if model_time eq 0 or model_time eq 2 then if i eq 0 then (*pre_emph[i])[0:n_elements(*pre_emph[i])-1:2]*=100.0/(6.0*C0*Imax_B0) else (*pre_emph[i])[0:n_elements(*pre_emph[i])-1:2]*=100.0
     if model_time eq 1 then if i eq 0 then (*pre_emph[i])[0:n_elements(*pre_emph[i])-1:4]*=100.0/(6.0*C0*Imax_B0) else (*pre_emph[i])[0:n_elements(*pre_emph[i])-1:4]*=100.0
     if model_time eq 2 then (*oscns[i])[0:n_elements(*oscns[i])-1:4]*=100.0

     text_output_time[i+4]+= strtrim(string(*pre_emph[i],format='('+strtrim(n_elements(*pre_emph[i]))+'(F8.3,:))'),2) + tabch
     if model_time eq 2 then  text_output_time[i+4]+= strtrim(string(*oscns[i],format='('+strtrim(n_elements(*oscns[i]))+'(F8.3,:))'),2) + tabch + strtrim(string(r2_time[i],format='(F8.3)'),2) $
     else text_output_time[i+4]+=strtrim(string(r2_time[i],format='(F8.3)'),2)

  endfor
  
  case project.procPramArray[CI].Bez_axis of  
    0 : ax = 'X'
    1 : ax = 'Y'
    2 : ax = 'Z'
    else : ax = axis
  endcase
  
  h1 = plot(T[0:Nt-1],B0,"b2+",layout = [2,2,1],window_title = ax+' gradient eddy current amplitudes over time',/xlog)
  h1fit = plot(Tc,fit_result[0,*],"r2-", NAME = 'B0',/current,/overplot,layout = [2,2,1],/xlog, margin = 0.2)
    
  h2 = plot(T[0:Nt-1],1e3*G[*,0],"b2+",/current,layout = [2,2,2],/xlog)
  h2fit = plot(Tc,1e3*fit_result[1,*],"r2-", NAME = 'Gx',/current,/overplot,layout = [2,2,2],/xlog, margin = 0.2)
  
  h3 = plot(T[0:Nt-1],1e3*G[*,1],"b2+",/current,layout = [2,2,3],/xlog)
  h3fit = plot(Tc,1e3*fit_result[2,*],"r2-", NAME = 'Gy',/current,/overplot,layout = [2,2,3],/xlog, margin = 0.2)
  
  h4 = plot(T[0:Nt-1],1e3*G[*,2],"b2+",/current,layout = [2,2,4],/xlog)
  h4fit = plot(Tc,1e3*fit_result[3,*],"r2-", NAME = 'Gz',/current,/overplot,layout = [2,2,4],/xlog, margin = 0.2)
  
  ax1 = h1.axes
  ax1[0].TITLE = 'Time (ms)'
  ax1[1].TITLE = '$B_0$ ($\mu$T) '
  h1.FONT_SIZE = 20
  
  ax2 = h2.axes
  ax2[0].TITLE = 'Time (ms)'
  ax2[1].TITLE = '$G_x$ ($\mu$T/m)'
  h2.FONT_SIZE = 20
  
  ax3 = h3.axes
  ax3[0].TITLE = 'Time (ms)'
  ax3[1].TITLE = '$G_y$ ($\mu$T/m)'
  h3.FONT_SIZE = 20
   
   ax4 = h4.axes
   ax4[0].TITLE = 'Time (ms)'
   ax4[1].TITLE = '$G_z$ ($\mu$T/m)'
   h4.FONT_SIZE = 20
    
  topbase1 = widget_base(title = axis+' gradient linear fit parameters' ,tlb_frame_attr=1,/grid_layout, frame = 1)
  mainbase1 = widget_base(topbase1)
  rowlabel1 = widget_text(mainbase1, value = text_output_space, ysize = Nt+3)
  widget_control, topbase1, /realize
  
  topbase2 = widget_base(title = axis+' gradient pre-emphasis constants' ,tlb_frame_attr=1,/grid_layout, frame = 1)
  mainbase2 = widget_base(topbase2)
  rowlabel2 = widget_text(mainbase2, value = text_output_time, ysize = 7)
  widget_control, topbase2, /realize
  
  close,/all  
end

; Subroutine name: Tensor_Bez
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To generate the ellipsoid for the impulse response tensor
;
; Editing Information:

pro Tensor_Bez

common scan_data
common Bez_params
CI = project.ci

if scheme ne 'Six Point' and scheme ne 'Four Point' and scheme ne 'Hadamard' then return

it = project.procPramArray[CI].Bez_time_index
Nt = n_elements(T)
B0x = make_array(Nt)
Gx  = make_array(Nt,3)
SpatialFit, CI, B0x, Gx, ax = 0, flag, it = it
if flag eq 0 then return

B0y = make_array(Nt)
Gy  = make_array(Nt,3)
SpatialFit, CI, B0y, Gy, ax = 1, flag, it = it
if flag eq 0 then return

B0z = make_array(Nt)
Gz  = make_array(Nt,3)
SpatialFit, CI, B0z, Gz, ax = 2, flag, it = it
if flag eq 0 then return

G = -1e3*[GX[it,*],GY[it,*],GZ[it,*]]
det = determ(G,/CHECK)
if det eq 0 then begin
  dummy = DIALOG_MESSAGE('Singular matrix. Please choose a different time point', /ERROR)
  return
endif
eval =  HQR(ELMHES(G), /DOUBLE)
evec = EIGENVEC(G, eval)
FOR i=0,2 DO evec[*,i] *= ABS(evec[0,i])/evec[0,i]
PRINT, 'Eigen value residuals:'
FOR i=0,2 DO print, G ## evec[*,i] - eval[i]*evec[*,i], $
  FORMAT ='(4("(",g9.2,",",g9.2,") "))'
       
dj_ellipsoid_ng, eval, transpose(evec),  window_title_txt = 'Gradient distortion tensor at t = '+strtrim(string(T[it],format = '(F5.2)'),2) + ' ms'
  
  
close,/all    
end

; Subroutine name: Bez_single_display
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display the eddy current field for a single slice at a given time
;
; Editing Information:

pro Bez_single_display

common scan_data
common Bez_params

CI = project.ci
gammabar = 42.58                      ; Proton gyromagnetic ratio (MHz/T)

Bez = load_state2_data_Bez(project.dataarray[CI].Bez)
if ~ptr_valid(Bez) then return

FOVr       = project.imndarray[CI].f_fov              ; RO-FOV (cm)
FOVp       = project.imndarray[CI].p_fov              ; PE-FOV (cm)
FOVs       = project.imndarray[CI].s_fov              ; PE2-FOV (cm)
slice_axis = project.procPramArray[CI].slice_axis     ; Slice axis (Read/Phase,Read/Slice or Phase/Slice)
it         = project.procPramArray[CI].Bez_time_index
TE         = project.imndarray[CI].echo_time  
case scheme of
  'Four Point' : xi = 1

  'Hadamard'  : xi = 4

   else       : xi = 2
endcase
range = 1e3*[-1,1]/float(xi*gammabar*TE)
case slice_axis of
  0 : begin
      is = project.procPramArray[CI].sdim_start
      iimage, abs((*Bez)[*,*,is,it,1]), aspect_ratio = FOVp/FOVr 
      iimage, (*Bez)[*,*,is,it,0], aspect_ratio = FOVp/FOVr, rgb_table=33, /current, /overplot, transparency=50, min_value = range[0], max_value = range[1], current_zoom = 4
      end
  1 : begin
        ip = project.procPramArray[CI].pdim_start
        iimage, reform(abs((*Bez)[*,ip,*,it,1])), aspect_ratio = FOVs/FOVr
        iimage, reform((*Bez)[*,ip,*,it,0]), aspect_ratio = FOVs/FOVr, rgb_table=33, /current, /overplot, transparency=50, min_value = range[0], max_value = range[1], current_zoom = 4
      end 
  2 : begin
        ir = project.procPramArray[CI].fdim_start
        iimage, reform(abs((*Bez)[ir,*,*,it,1])), aspect_ratio = FOVs/FOVp
        iimage, reform((*Bez)[ir,*,*,it,0]), aspect_ratio = FOVs/FOVp, rgb_table=33, /current, /overplot, transparency=50, min_value = range[0], max_value = range[1], current_zoom = 4
      end           
endcase

if ptr_valid(Bez) then ptr_free,Bez
close,/all
end

; Subroutine name: Bez_multi_display
; Created by: Magdoom Kulam
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display eddy current field for all slices or all times for a given slice
;
; Editing Information:

pro Bez_multi_display

common scan_data
common Bez_params

CI = project.ci
gammabar = 42.58                      ; Proton gyromagnetic ratio (MHz/T)

Bez = load_state2_data_Bez(project.dataarray[CI].Bez)
if ~ptr_valid(Bez) then return
    
FOVr     = project.imndarray[CI].f_fov            ; RO-FOV (cm)
FOVp1    = project.imndarray[CI].p_fov            ; PE-FOV (cm)
FOVp2    = project.imndarray[CI].s_fov            ; PE2-FOV (cm)
Nr       = project.imndarray[CI].fdim
Np1      = project.imndarray[CI].pdim
Np2      = project.imndarray[CI].sdim

; Zeropad effect
if project.procPramArray[CI].zpad_flag eq 1 then begin
  Nr*=2           
  Np1*=2                 
  if project.imndArray[CI].dimensions eq 3 then Np2*=2              
endif 
Nt = n_elements(T)
TE         = project.imndarray[CI].echo_time
case scheme of
  'Four Point' : xi = 1

  'Hadamard'  : xi = 4

  else       : xi = 2
endcase
range = 0.5*1e3*[-1,1]/float(2.0*xi*gammabar*TE)

sort_dir = project.procPramArray[CI].sort_dir     ; Sort direction (Space or time)
slice_axis = project.procPramArray[CI].slice_axis ; Slice axis (Read/Phase,Read/Slice or Phase/Slice)
case sort_dir of
  0 : begin
      it = project.procPramArray[CI].Bez_time_index
      case slice_axis of
        0   : begin
              for i = 1,Np2 do begin
                if i eq 1 then iimage, abs((*Bez)[*,*,i-1,it,1]), aspect_ratio = FOVp1/FOVr, layout=[1+floor(Np2/floor(sqrt(Np2))),floor(sqrt(Np2)),i], window_title = 'Bez (uT) variation in space at t = ' + strtrim(string(T[it],format = '(F6.2)'),2) + 'ms'  $
                else iimage, abs((*Bez)[*,*,i-1,it,1]), aspect_ratio = FOVp1/FOVr, /current, layout=[1+floor(Np2/floor(sqrt(Np2))),floor(sqrt(Np2)),i]
                iimage, (*Bez)[*,*,i-1,it,0], aspect_ratio = FOVp1/FOVr, rgb_table=33, /current, /overplot, transparency=50, layout=[1+floor(Np2/floor(sqrt(Np2))),floor(sqrt(Np2)),i], min_value = range[0], max_value = range[1]
              endfor
              end
        1   : begin
                for i = 1,Np1 do begin
                  if i eq 1 then iimage, reform(abs((*Bez)[*,i-1,*,it,1])), aspect_ratio = FOVp2/FOVr, layout=[1+floor(Np1/floor(sqrt(Np1))),floor(sqrt(Np1)),i], window_title = 'Bez (uT) variation in space at t = ' + strtrim(string(T[it],format = '(F6.2)'),2) + 'ms'  $
                  else iimage, reform(abs((*Bez)[*,i-1,*,it,1])), aspect_ratio = FOVp2/FOVr, /current, layout=[1+floor(Np1/floor(sqrt(Np1))),floor(sqrt(Np1)),i]
                  iimage, reform((*Bez)[*,i-1,*,it,0]), aspect_ratio = FOVp2/FOVr, rgb_table=33, /current, /overplot, transparency=50, layout=[1+floor(Np1/floor(sqrt(Np1))),floor(sqrt(Np1)),i], min_value = range[0], max_value = range[1]
                endfor
              end
        2   : begin
                for i = 1,Nr do begin
                  if i eq 1 then iimage, reform(abs((*Bez)[i-1,*,*,it,1])), aspect_ratio = FOVp2/FOVp1, layout=[1+floor(Nr/floor(sqrt(Nr))),floor(sqrt(Nr)),i], window_title = 'Bez (uT) variation in space at t = ' + strtrim(string(T[it],format = '(F6.2)'),2) + 'ms' $
                  else iimage, reform(abs((*Bez)[i-1,*,*,it,1])), aspect_ratio = FOVp2/FOVp1, /current, layout=[1+floor(Nr/floor(sqrt(Nr))),floor(sqrt(Nr)),i]
                  iimage, reform((*Bez)[i-1,*,*,it,0]), aspect_ratio = FOVp2/FOVp1, rgb_table=33, /current, /overplot, transparency=50, layout=[1+floor(Nr/floor(sqrt(Nr))),floor(sqrt(Nr)),i], min_value = range[0], max_value = range[1]
                endfor
              end
        endcase                       
      end
  1  : begin
      for i = 1,Nt do begin
        case slice_axis of
          0: begin 
             is = project.procPramArray[CI].sdim_start
             z = 10*(is-Np2/2.0)*FOVp2/Np2
             if i eq 1 then iimage, abs((*Bez)[*,*,is,i-1,1]), aspect_ratio = FOVp1/FOVr, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i], window_title = 'Bez (uT) variation in time at PE2 position = ' + strtrim(string(z,format = '(F6.2)'),2) + ' mm'  $
             else iimage, abs((*Bez)[*,*,is,i-1,1]), aspect_ratio = FOVp1/FOVr, /current, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i]
             iimage, (*Bez)[*,*,is,i-1,0], aspect_ratio = FOVp1/FOVr, rgb_table=33, /current, /overplot, title = 't = '+strtrim(string(T[i-1],format='(F6.2)'),2)+' ms', font_size = 12, transparency=50, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i], min_value = range[0], max_value = range[1]
             end
          1: begin
             ip = project.procPramArray[CI].pdim_start
             y = 10*(ip-Np1/2.0)*FOVp1/Np1
             if i eq 1 then iimage, reform(abs((*Bez)[*,ip,*,i-1,1])), aspect_ratio = FOVp2/FOVr, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i], window_title = 'Bez (uT) variation in time at PE1 position = ' + strtrim(string(y,format = '(F6.2)'),2) + ' mm' $
             else iimage, reform(abs((*Bez)[*,ip,*,i-1,1])), aspect_ratio = FOVp2/FOVr, /current, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i]
             iimage, reform((*Bez)[*,ip,*,i-1,0]), aspect_ratio = FOVp2/FOVr, rgb_table=33, /current, /overplot, title = 't = '+strtrim(string(T[i-1],format='(F6.2)'),2)+' ms', font_size = 12, transparency=50, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i], min_value = range[0], max_value = range[1]
             end    
          2: begin
             ir = project.procPramArray[CI].fdim_start
             x = 10*(ir-Nr/2.0)*FOVr/Nr
             if i eq 1 then iimage, reform(abs((*Bez)[ir,*,*,i-1,1])), aspect_ratio = FOVp2/FOVp1, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i], window_title = 'Bez (uT) variation in time at RO position = ' + strtrim(string(x,format = '(F6.2)'),2) + ' mm'$
             else iimage, reform(abs((*Bez)[ir,*,*,i-1,1])), aspect_ratio = FOVp2/FOVp1, /current, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i]
             iimage, reform((*Bez)[ir,*,*,i-1,0]), aspect_ratio = FOVp2/FOVp1, rgb_table=33, /current, /overplot, title = 't = '+strtrim(string(T[i-1],format='(F6.2)'),2)+' ms', font_size = 12, transparency=50, layout=[1+floor(Nt/floor(sqrt(Nt))),floor(sqrt(Nt)),i], min_value = range[0], max_value = range[1]
             end  
          endcase              
      endfor
      end     
endcase

if ptr_valid(Bez) then ptr_free,Bez
close,/all
end

; Subroutine name: Bez_movie
; Created by: Magdoom Kulam
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To display the eddy current field maps as a movie
;
; Editing Information:

pro Bez_movie

common scan_data

CI = project.ci
Bez = load_state2_data_Bez(project.dataarray[CI].Bez)
if ~ptr_valid(Bez) then return

N = size(*Bez,/Dimensions)
sort_dir = project.procPramArray[CI].sort_dir ; Sort direction (Space or time)

top_base = widget_base(title = 'Eddy Current Field',/grid_layout)

case sort_dir of
  0 : begin
      progressbar = Obj_New('progressbar', Color='red', Text='Loading frames', /NOCANCEL)
      progressbar -> Start
      anim = cw_animate(top_base,N[0],N[1],N[2], /NO_KILL)
      DEVICE, DECOMPOSED = 0
      loadCT,13
      it = project.procPramArray[CI].Bez_time_index
      for i=0,N[2]-1 do begin
        cw_animate_load, anim, frame = i, image = bytscl((*Bez)[*,*,i,it,0], min = -1, max = 1)
        progressBar -> Update, (float(i+1)/float(N[2]))*100.0
      endfor
      progressbar -> Destroy
      loadCT,0  
      end 
   1 : begin
       progressbar = Obj_New('progressbar', Color='red', Text='Loading frames', /NOCANCEL)
       progressbar -> Start
       anim = cw_animate(top_base,N[0],N[1],N[3], /NO_KILL)
       DEVICE, DECOMPOSED = 0
       loadCT,13
       is = project.procPramArray[CI].sdim_start
       for i=0,N[3]-1 do begin
        cw_animate_load, anim, frame = i, image = bytscl((*Bez)[*,*,is,i,0], min = -1, max = 1)
        progressBar -> Update, (float(i+1)/float(N[3]))*100.0
       endfor
       progressbar -> Destroy
       loadCT,0
      end      
endcase 

widget_control, /realize, top_base
cw_animate_run, anim,10
XMANAGER, 'CW_ANIMATE Demo', top_base, EVENT_HANDLER = 'anim_ehandler'

if ptr_valid(Bez) then ptr_free,Bez
close,/all
end

; Subroutine name: Bez_gui_cleanup
; Created by: Magdoom Kulam
; Calling Information:
;
; top_base
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Cleanup procedure before exiting the eddy current GUI
;
; Editing Information:

pro Bez_gui_cleanup, topbase

  widget_control, topbase, get_uvalue=state
  if (ptr_valid(state)) then ptr_free, state
  return

end