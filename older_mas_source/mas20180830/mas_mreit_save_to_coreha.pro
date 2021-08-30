
PRO mas_mreit_save_to_coreha, dataptr, magdataptr

    COMMON scan_data
    Bzdata = *dataptr
    MRdata = *magdataptr
    
    sz = size(Bzdata) ; Usually 3D data.
    
    fdim = sz[1]
    pdim = sz[2]
    sdim = sz[3]
    
    MRfolder = DIALOG_PICKFILE(TITLE='Please select the folder where you want to save the magnitude images',/directory)
    FOR k = 0, sdim-1 DO BEGIN
      IF(k lt 9) THEN BEGIN
        fname = MRfolder + '00' + strtrim(string(k+1),1) + '.mri'
      ENDIF ELSE BEGIN
        fname = MRfolder + '0'+ strtrim(string(k+1),1) + '.mri'
      END
      OPENW, 1, fname
      WRITEU, 1, MRdata[*,*,k]
      CLOSE, 1
    END
    
    Bzfolder = DIALOG_PICKFILE(TITLE='Please select the folder where you want to save the Bz field maps',/directory)
    FOR k = 0, sdim-1 DO BEGIN
      IF(k lt 9) THEN BEGIN
        fname = Bzfolder + '00' + strtrim(string(k+1),1) + '.mri'
      ENDIF ELSE BEGIN
        fname = Bzfolder + '0' + strtrim(string(k+1),1) + '.mri'  
      END
      OPENW, 1, fname
      WRITEU, 1, Bzdata[*,*,k]
      CLOSE, 1
    END      
      
    mesg = DIALOG_MESSAGE('The images have been exported',/info)

    

END