

FUNCTION get_filename, maps_path
    
     res = maps_path
     res_parts = STRSPLIT(res, PATH_SEP(), /EXTRACT)
     sz = SIZE(res_parts)
     index = sz[3] - 1  ; since res_parts will always be a vector, get the last, $
                        ; for number of elements, and minus 1 for index.
     fname = res_parts[index]
     RETURN, fname
 
END



FUNCTION set_map_flags, fnames
   
    COMMON scan_data
       FLAGS = MAKE_ARRAY(5)
       map_list = ['AD_map.nii.gz','ProtonDensity_map.nii.gz',$
                    'S0_map.nii.gz','T1_map.nii.gz','T2_map.nii.gz']
         FOR i = 0,4 DO BEGIN     
            index = WHERE(STRMATCH(fnames,map_list(i)))
            IF(index eq -1) THEN BEGIN
                FLAGS(i) = 0
            ENDIF ELSE BEGIN
                FLAGS(i) = 1
            ENDELSE
            PRINT, fnames(i)
            PRINT, index
         END 
       RETURN, FLAGS     
       
END



PRO trial_open_qmri_maps_dir

    COMMON scan_data
    COMMON INDEX, flag_array, maplist_start, maplist_end
      dirname = DIALOG_PICKFILE(/directory)
      PRINT, dirname
        IF (dirname eq '') THEN BEGIN
          RETURN
        ENDIF     
        files = FILE_SEARCH(dirname, '*.nii.gz')
        sz = SIZE(files)
        fnames = MAKE_ARRAY(5, /STRING)
        
        IF (sz[1] lt 6) THEN BEGIN
           FOR i=0, sz[1]-1 DO BEGIN
                fnames[i]=get_filename(files[i])
           END
        ENDIF ELSE BEGIN
           mesg = DIALOG_MESSAGE(['Detected more Nifti files than max number of maps.',$
                                  'Please check the data in the folder and ensure only maps exist.'], /error, /center)
           RETURN
        ENDELSE
        flag_array = set_map_flags(fnames)     
              IF((flag_array[1] ne 1) or (flag_array[3] ne 1) or (flag_array[4] ne 1)) THEN BEGIN
                  res = DIALOG_MESSAGE(['Please make sure at least the T1, T2 and ProtonDensity NIFTI maps are present in the folder.',$ 
                                    'Also please ensure they are sorted alphabetically.'], /error, /center)
                  RETURN   
              ENDIF         
        
        ni = project.ni             
              FOR i = 0 , sz[1]-1 DO BEGIN
                  nifti = mas_read_nifti(nifti_filename=files[i])
                  mnf_nifti_hdr_to_mas, nifti
                  mas_load_state_1
                  ni += 1
                  PRINT, 'The value of ni is : ', ni
                  fnames[i] = get_filename(files[i])
              END              
        maplist_end = ni-1  ; NOTE: ni is the index for the scan that follows.
        maplist_start = ni - sz[1]
        PRINT, flag_array
        PRINT, maplist_end , maplist_start
        
            
END
