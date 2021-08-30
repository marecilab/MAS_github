;+
; NAME:
;   mas_varian_procpar
;
; AUTHOR:
;   Bill Triplett (wtriplett@gmail.com)
;
; PURPOSE:
;   Provide an object-oriented interface to Vnmrj's 'procpar' file.
;   This file contains parameters associated with data acquired on
;   varian's Vnmrj console. This class provides access to the values.
;
; CREATION AND USE:
;       Create:
;           oProcpar = obj_new('mas_varian_procpar', '/path/to/procpar')
;
;       Lookup a value (method 1):
;           pvalue = oProcpar->lookup(paramname, is_array=ia, not_found=nf)
;           if (nf eq 1) then print, "Not Found (pvalue will be -1)."
;
;       Lookup a value (method 2):
;           boolean = oProcpar->paramExists('pslabel', value=pvalue)
;           if (boolean eq 1) then begin
;               print, 'The value of pslabel = "+pvalue
;           endif else begin
;               print, 'pslabel not found!!!'
;           endelse
;
;       Compute the data dimensions:
;           data_dims = oProcpar->computeDataDimensions()
;           print, "data dims = ", data_dims[0:3]
;
;       Dump a nice list of parameters and values in alphabetical
;       order.
;       Useful for finding differences (using unix diff(1)):
;           oProcpar->dumpParameters, toFile='param_dump.txt'
;
; EXAMPLE:
;       See CREATION AND USE
;
; CHANGES:
;       20100802 - Initial revision
;       20101015 - Updated to compute data dimensions for 1D data sets
;-

function mas_varian_procpar::init, filename

    ;; set up the check to make sure file exists and is specified
    ;; then parse it.

    if (n_elements(filename) ne 0 && file_test(filename, /read)) then begin
    
        self->parseFile, filename
        
    endif else begin
    
        message, "The name of an existent file is required."
        
    endelse
    
    message, "parsed: "+filename, /informational
    
    return, 1
    
end

pro mas_varian_procpar::cleanup

    ;; destroy the object's data members
    message, "cleaning up...", /informational
    ptr_free, self.procpar
    ptr_free, self.tags
    
end

;; Looks up the parameter named by 'param'. Optionally sets is_array to 1
;; if the parameter is an arrayed parameter. Also sets not_found to 1
;; if the parameter is not found, 0 otherwise. Note that if a parameter
;; is not found, the method will return -1.
function mas_varian_procpar::lookup, param, is_array=is_array, not_found=not_found

    ;; Lookup a parameter. function returns -1 and sets not_found=1
    ;; if the parameter is not found

    names = tag_names(*self.procpar)
    n_names = n_elements(names)
    if (names[n_names-1] ne '__ALL_TAGS') then begin
        not_found = 1
        return, -1
    endif
    
    if (n_elements(param) eq 1) then begin
        for n = 0, n_elements(*self.tags)-1 do begin
            if (param eq (*self.tags)[n]) then begin
                not_found = 0
                is_array = n_elements((*self.procpar).(n)) gt 1 ? 1 : 0
                return, (*self.procpar).(n)
            endif
        endfor
    endif else if (n_elements(param) gt 1) then begin
    
        message, "Arrayed lookups are not supported."
        
    ;; I am still considering how to implement this.
    ;; the main problem is that we need to return a list
    ;; of items of mixed type, and right now either a
    ;; pointer array or struct is the way to go.
        
    ;        start = 0
    ;        not_found = intarr(n_elements(param)) + 1
    ;        for p = 0, n_elements(param)-1 do begin
    ;
    ;            for n = 0, n_elements(*self.tags)-1 do begin
    ;
    ;                if (param[p] eq (*self.tags)[n]) then begin
    ;
    ;                    if (start eq 0) then begin
    ;
    ;                        out = create_struct(names[n], (*self.procpar).(n))
    ;                        start = 1
    ;
    ;                    endif else begin
    ;
    ;                        out = create_struct(out, names[n], (*self.procpar).(n))
    ;
    ;                    endelse
    ;
    ;                    not_found[p] = 0
    ;
    ;                endif
    ;
    ;            endfor
    ;
    ;        endfor
    ;
    ;        if (n_elements(out) ne 0) then return, out
    ;        return, -1
        
    endif else begin
    
        message, "a parameter name must be specified."
        
    endelse
    
    not_found = 1
    return, -1
    
end

;; See Varian Vnmrj User Programming Manual, p266
;; This method is internal and should not be called
;; outside of the class.
pro mas_varian_procpar::parseFile, filename

    ;; Note, parameters can be case sensitive (for example:
    ;; tDELTA and tdelta, and te and TE). Also there
    ;; is at least one parameter name that clashes with an
    ;; idl reserved work (ne). Therefore, for each parameter
    ;; in the procpar file, we will append a '_' to the parameter
    ;; name. If we find case-conflict, we will append a second '_'
    ;; to that parameter. So, 'te' becomes 'te_', then later if
    ;; 'TE' is found it will become 'te__'. The true, case-sensitive
    ;; paramter name is stored in a separate array. The indices of
    ;; the array are in 1:1 correspondence to the order of the
    ;; structure tags.

    if (filename eq '') then return
    
    openr, lun, filename, /get_lun
    
    paramln = ""
    dataln  = ""
    temp    = ""
    
    procpar = { procpar_file: filename, $
        parent_dir: file_dirname(filename) }
    procpar_tags = [ 'procpar_file', 'parent_dir' ]
    
    while (not eof(lun)) do begin
    
        readf, lun, paramln
        ;;;;; check if we have a param line: time_submitted 2 2 8 0 0 2 1 0 1 64
        newvar = stregex(paramln, "^([a-zA-Z0-9_]+) ([0-9]+) ([0-9]+) ([-0-9\.e\-\+]+) " + $
            "([-0-9\.e\-\+]+) ([-0-9\.e\-\+]+) ([0-9]+) ([0-9]+) " + $
            "([0-9]+) ([0-9]+) ([0-9]+)", /boolean)
        if (newvar eq 0) then continue
        
        vardata = strsplit(paramln, " ", /extract)
        
        if (vardata[1] eq '2' or vardata[2] eq '2') then begin
        
            readf, lun, paramln
            
            ;; matches a procpar string array, which can often be a multi-line
            ;; data structure. The first line contains the number of elements
            ;; and the first element. The subsequent lines contain the rest
            ;; of the elements.
            temp = stregex(paramln, '^([0-9]+) "([^"]+)"', /subexpr, /extract)
            num_items = long(temp[1])
            data = temp[2]
            if (num_items gt 1) then begin
                for i = 1, num_items-1 do begin
                    readf, lun, paramln
                    ;; remove quotation marks
                    data = [ data, strmid(paramln, 1, strlen(paramln)-2) ]
                endfor
            endif
            
        endif else if (vardata[1] eq '3') then begin
        
            readf, lun, paramln
            temp = strsplit(paramln, ' ', /extract)
            nitems = long(temp[0])
            data = (nitems gt 1) ? float(temp[1:*]) : float(temp[1])
            
        endif else if (vardata[1] eq '7') then begin
        
            readf, lun, paramln
            temp = strsplit(paramln, ' ', /extract)
            nitems = long(temp[0])
            data = (nitems gt 1) ? long(temp[1:*]) : long(temp[1])
            
        endif else if (vardata[1] eq '1' or vardata[2] eq '1') then begin
        
            readf, lun, paramln
            temp = strsplit(paramln, ' ', /extract)
            nitems = long(temp[0])
            data = (nitems gt 1) ? float(temp[1:*]) : float(temp[1])
            
        endif else begin
        
            print, "Unknown Data Type."
            continue
            
        endelse
        
        names = tag_names(procpar)
        new_name = strupcase(vardata[0]+'_')
        
        for n = 0, n_elements(names)-1 do begin
            if (new_name eq names[n]) then begin
                new_name = new_name+'_'
                break
            endif
        endfor
        
        procpar_tags = [ procpar_tags, vardata[0] ]
        procpar = create_struct(procpar, new_name, data)
        
    endwhile
    
    procpar_tags = [ procpar_tags, '__all_tags' ]
    
    procpar = create_struct(procpar, '__all_tags', procpar_tags)
    
    close, lun & free_lun, lun
    
    self.procpar = ptr_new(procpar, /no_copy)
    self.tags    = ptr_new(procpar_tags, /no_copy)
    
end

;; Checks if a parameter exists and optionall returns its value.
;; returns 1 if the parameter exists, 0 otherwise
;; arg1: param - string containing the paramer to be checked
function mas_varian_procpar::paramExists, param, value=value

    ;; check to see if a parameter exists, and if so optionally
    ;; return its value in a keyword.

    for i = 0, n_elements(*self.tags)-1 do begin
        if (param eq (*self.tags)[i]) then begin
            if (arg_present(value)) then begin
                value = (*self.procpar).(i)
            endif
            return, 1
        endif
    endfor
    
    return, 0
end

;; returns a 4 elements array with the following meaning;
;; dim[0] = dimension of readout
;; dim[1] = number of phaseencode steps
;; dim[2] = number of slices or number of second PE steps
;; dim[3] = number of volumes for arrayed experiments
function mas_varian_procpar::computeDataDimensions

    ;; Computes the FID data dimensions based on the
    ;; values stored in the procpar file.

    seqcon     = self->lookup('seqcon')
    ni         = self->lookup('ni')
    ni2         = self->lookup('ni2')
    np         = self->lookup('np')
    nf         = self->lookup('nf')
    arraydim   = self->lookup('arraydim')
    num_slices = self->lookup('ns')
    num_dims   = self->lookup('nD')
    num_echos  = self->lookup('ne')
    nv         = self->lookup('nv')
    nv2        = self->lookup('nv2')
    layout     = self->lookup('layout')
    seqfil     = self->lookup('seqfil')
    
    if (self->lookup('apptype') eq 'im2Depi') then begin
        fdim = self->lookup('nread')/2
        pdim = self->lookup('nphase')
        oversample = self->lookup('oversample', not_found=ntfd)
        if (ntfd) then oversample = 1
        
        sdim = num_slices
        adim = arraydim
        return, [fdim, pdim, sdim, adim]
    endif
        
    if (seqfil eq 'spuls') then begin
        return, [np/2, 1, 1, 1]
    endif
   
    seqcon = string(transpose(byte(seqcon)))
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; SEQCON CHECK                                                                  ;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if (strjoin(seqcon, '') eq 'nnnnn') then begin
        ;; data is 1D spect
        ;; we don't read that at this time
        message, '1D Spect not supported.'
        
    endif else if (seqcon[3] eq 'n')  then begin
        message, "Data set is 2D.",/informational
        
        fdim = np/2
        if seqcon[2] eq 'c' then pdim = nv else pdim = ni
        sdim = num_slices
        
        if seqcon[2] eq 's' then arraydim/=ni
        if seqcon[0] eq 'n' then adim = arraydim else adim = num_echos
      
    endif else if (seqcon[4] eq 'n') then begin
        message, "Data set is 3D.",/informational
        
        fdim = np/2
        if seqcon[2] eq 'c' then pdim = nv else pdim = ni 
        if seqcon[3] eq 'c' then sdim = nv2 else sdim = ni2
        
        if seqcon[2] eq 's' then arraydim/=ni
        if seqcon[3] eq 's' then arraydim/=ni2
        if seqcon[0] eq 'n' then adim = arraydim else adim = num_echos
     
    endif else begin
        ;; data is 4D
        message, '4D data (or whatever this is) is not supported at this time.'
    endelse
    
   if seqfil eq 'mge3d_Bez' or seqfil eq 'mge3d_Bez_v2' and arraydim gt 1 then return, [fdim, pdim, sdim, num_echos, arraydim] else return, [ fdim, pdim, sdim, adim ]
    
end

;; if toFile is not specified, it will dump the parameter list
;; to the console.
pro mas_varian_procpar::dumpParameters, toFile=toFile

    ;; prints a nice-looking parameter list to the console
    ;; or if toFile is set it will save it to a file.

    if (n_elements(toFile) ne 0) then begin
        openw, lun, toFile, /get_lun, error=iserror
        if (iserror ne 0) then begin
            message, "Unable to open "+toFile+" for writing." + $
                !ERROR_STATE.MSG, /continue
            return
        endif
    endif else begin
        lun = -1 ;; stdout
    endelse
    
    lengths   = strlen(*self.tags)
    maxlength = max(lengths)
    tagsort   = sort(*self.tags)
    
    for i = 0, n_elements(*self.tags)-1 do begin
    
        index = tagsort[i]
        
        ;; this tag is for internal use only.
        if ((*self.tags)[index] eq '__all_tags') then continue
        
        ;; figure out max length for spacing
        str = (*self.tags)[index]
        if (lengths[index] lt maxlength) then begin
            str = strjoin(replicate(' ', maxlength-lengths[index]), '')+str
        endif
        
        str += " = "
        
        case ( size((*self.procpar).(index), /structure) ).(0) of
            'STRING': format = ''
            'INT':    format = '(I0)'
            'LONG':   format = '(I0)'
            'FLOAT':  format = '(G0)'
            'DOUBLE': format = '(G0)'
            'BYTE':   format = '(I0)'
            else:     format = ''
        endcase
        
        if (n_elements((*self.procpar).(index)) gt 1) then begin
            str += "ARRAY: [ "+strjoin( string((*self.procpar).(index), format=format), ', ')+" ]"
        endif else begin
            str += string((*self.procpar).(index), format=format)
        endelse
        
        printf, lun, str
        
    endfor
    
    if (lun ne -1) then begin
        close, lun
        free_lun, lun
    endif
    
end

;; structure/class member definition
pro mas_varian_procpar__define

    struct = {  MAS_VARIAN_PROCPAR, $
        procpar: ptr_new(), $
        tags : ptr_new() $
        }
        
end



