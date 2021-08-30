;+
; :Description:
;    This procedure is for creating the Bz image if the MREIT ICNE approach 
;    was employed. It is called through the MREIT GUI (mag_field_calculate.pro) when the ICNE check box is selected.
;    
;    NOTE:
;    It is assumed that the two scans i.e positive and negative current injections are acquired with the same
;    bandwidth and acqtimes with same readout gradient amplitude. 
;    For Agilent files (.fid), the raw data is read directly by this code. For the other file types (nifti,dicom,bruker files)
;    TT IS BEST IF THE DATA LOADED IS THE RAW KSPACE DATA. 
;
; :Params:
;    infoptr
;    image = The final B image (3D) for all the slices in the dataset.
;
; :Keywords:
;    output1
;
; :Author: aditya
;-

PRO mas_mreit_icne, infoptr, output1 = image

  COMMON scan_data
  forward_function get_dir_slash
  CI = project.ci
  dir_sep = (get_dir_slash())[1]
  
  scan_list_array = [(*infoptr).scan1, (*infoptr).scan2, (*infoptr).scan3]
  print, scan_list_array
  dataptr = [ptr_new(), ptr_new(), ptr_new()]
  
  
  ;;;;;;;;;;;; READING FIDs OF POSITIVE AND NEGATIVE CURRENT INJECTION DATASETS.
  FOR i = 0,2 DO BEGIN
    procpar_path = project.imndarray[scan_list_array[i]].file_path
    opp = OBJ_NEW('mas_varian_procpar', procpar_path)
    
    seqfil = opp->lookup('seqfil')
    apptype = opp->lookup('apptype')
    
    OPENR, lun, FILE_DIRNAME(project.imndarray[scan_list_array[i]].file_path) + dir_sep+"fid",$
                   /swap_if_little_endian, /get_lun, error=iserror
            if (iserror) then begin
                void = DIALOG_MESSAGE(['Unable to open FID file.', $
                                       !ERROR_STATE.MSG], /error, /center)
            endif
            dataptr[i] = mas_varian_read_fid_data_2d(lun, opp)
            help, dataptr
    CLOSE, lun
    OBJ_DESTROY, opp
  ENDFOR
  
  
  ;;;;;;;;;;;;;;;;;;;;;;;; COMPLEX KSPACE IMAGES OF THE POSITIVE AND NEGATIVE CURRENT DATASETS.
  
  rdim = (*infoptr).rdim                                         ; frequency encode dimension size
  pdim = (*infoptr).pdim                                         ; phase encode dimension size
  sdim = (*infoptr).sdim                                         ; number of slices
  slicepos = project.procpramarray[(*infoptr).scan1].sdim_start  ; slice number from the main MAS window
  Gread = project.imndarray[(*infoptr).scan2].Gfe                ; Readout Gradient strength from the procpar file (G/cm)  
  Ts = project.imndarray[(*infoptr).scan2].acqtime               ; Sampling time extracted in seconds
  Tc = (*infoptr).Tc_value                                       ; Injection current duration in seconds
  Tc_new = Tc + Ts                                               ; Total duration of injection current in seconds
  gyro_ratio = (*infoptr).Gamma_value                            ; Gyromagnetic ratio (rad/s/T)
  f_fov = project.imndarray[(*infoptr).scan2].f_fov              ; Readout field of view to compute delta kx. (in cm)
  
  kspace_complex_imagepos = *dataptr[1]                         ; positive current complex image 
  kspace_complex_imageneg = *dataptr[2]                         ; negative current complex image
  kspace_complex_no_current = *dataptr[0]
  help, kspace_complex_no_current
  
  FOR i = 0, sdim-1 DO BEGIN
      kspace_complex_no_current[*,*,i] = shift(fft(shift(kspace_complex_no_current[*,*,i],rdim/2,rdim/2),/inverse),-rdim/2,-rdim/2)
  END  

  no_current_image = (kspace_complex_no_current)                  ; Image domain dataset. Type: complex
  help, no_current_image
  
  i_value = complex(0,1)                            ; imaginary number - i
  complex_image = fltarr(rdim,pdim)
  fft_image = fltarr(rdim,pdim)
  field_image = fltarr(rdim,pdim)
  final_image = fltarr(rdim,pdim,sdim)

  tmin = 0
  tmax = Ts
  delta_t = linspace(tmin,tmax,rdim)

  help, delta_t
  
  FOR slice = 0,sdim-1 DO BEGIN
    FOR i = 0,pdim-1 DO BEGIN
       FOR j = 0,rdim-1 DO BEGIN
           complex_image[j,i] = (kspace_complex_imagepos[j,i,slice] - kspace_complex_imageneg[j,i,slice])/(2*i_value*(gyro_ratio*(Tc_new - Ts/2 + delta_t(j))))
       ENDFOR
    ENDFOR

    fft_image = shift(fft(shift(complex_image,rdim/2,rdim/2),/inverse),-rdim/2,-rdim/2)
    field_image = (fft_image)/no_current_image[*,*,slice]
    final_image[*,*,slice] = (field_image)
  ENDFOR

image = final_image
help, image
END