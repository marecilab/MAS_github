;
; NAME:
; 	DICOM_WRITER
;
; VERSION:
;	0.21
;
; PURPOSE:
;	Generate a dicom file from within RSI IDL
; 
; AUTHOR:
;	Bhautik Joshi
;
; EMAIL:
;	bjoshi@geocities.com
;
; HOMEPAGE:
;	http://cow.mooh.org
;
; USE:
;	DICOM_WRITER, filename, image, VOXELSIZE = voxelsize, SSAI = ssai, $
;		PATIENT = patient, PHYSICIAN = physician, PATID = patid
;
; INPUT:
;	filename - string containing name of dicom file to be written to
;	image - Integer (BYTE, FIX, UINT, LONG or ULONG) image - type and bpp
;	is now automagically set
;
; OPTIONAL PARAMETERS
;	voxelsize - Array of 3 floating point values representing voxel size
;	with the format [x,y,z], otherwise set to default of [1.0,1.0,1.0]
;	ssai - Array of 4 integer values representing studyID, seriesnum,
;       acqnum,imagenum, with the format [studyID,seriesnum,acqnum,imagenum],
;	otherwise set to default of [0,0,0,0]
;	patient - patient name, if not defined, set to dummy name
;	physician - physician name, if not defined, set to dummy name
;	patid - patiend ID, if not defined, set to dummy name
;
; NOTES ON USAGE (READ! IMPORTANT!):
;	* At the moment the program only writes to a single slice
;	* Extra dicom tags can be easily added (see body of program, especially
;	  generate_VRtag function)
;	* There is little to no error-checking at the moment, so be 
;	  careful!
;	* Analyse seems to need a minimum image size of somewhere around
;	  100x100
;	* IMPORTANT: The DICOM writer tries to write 'Implicit VR' type
;	  DICOM files - see DICOM standard PS 3.5-1999, part 7
;	* Can write most VR (Value Represenation) tags via new function,
;	  generate_VRtag. Currently supported:
;	  AE, AS, AT, CS, DA, DS, DT, FL, FD, FD, IS, LO, LT, OB, OW,
;	  SH, SL, SS, ST, TM, UI, UL, UN, US, UT
;	  and SQ, PN unsupported (I got away with using UI in place of PN)
;	* See comments near generate_VRtag function for notes on usage
;	  of the function for adding your own additional tags
;
; EXAMPLE:
;	Create a horrendously boring byte image and store it in a 
;	dicom file, test.dcm, with voxel dimensions of [2.0,3.0,4.0],
;	and studyid=1,series number=2,acquisiton number=3 and image
;	number=4:
;
;	> rows = 200
;	> cols = 200
;	> image = indgen(rows,cols)
;	> dicom_writer, 'test.dcm', image, voxelsize=[2.0,3.0,4.0], ssai=[1,2,3,4]
;
; HISTORY:
;	Based on Marc O'Briens (m.obrien@sghms.ac.uk) TIFF_to_DICOM.c
;	version 0.1 	08-01-2002 - first working version produced
;	version 0.11	09-01-2002 - fixed endian-ness issue & added get_lun
;				     functionality
;	version 0.2	14-01-2002 - many fixes and additions, including:
;		* replaced most generate_* functions with generate_VRtag
;		* support for many VR (Value representation) types
;		* Autodetection of little/big endian architecture and
;		  automagic byte ordering as necessary (for tags/image)
;		* automagically detect image type and set bpp as necessary
;		* more data in the header can be set manually
;	version 0.21	15-01-2002 - uploaded all over the place & fixed bug
;				     that didn't update patient, patid etc.
;
; TODO:
;	* Allow for more robust dicom writing
;	* Part 10 compliance (!!!!!!!!!!!)
;	* Decent error checking
;
; DISCLAIMER:
; 
; Permission to use, copy, modify, and distribute this software and its
; documentation for any purpose and without fee is hereby granted,
; provided that the above copyright notice appear in all copies and that
; both that copyright notice and this permission notice appear in
; supporting documentation.
;
; This file is provided AS IS with no warranties of any kind.  The author
; shall have no liability with respect to the infringement of copyrights,
; trade secrets or any patents by this file or any part thereof.  In no
; event will the author be liable for any lost revenue or profits or
; other special, indirect and consequential damages.
; 
; The author accepts no responsibility for any action arising from use of 
; this package. The software is not guaranteed to write compliant DICOM
; files. If it causes damage to you or your system, you have been warned -
; this is a work in progress. If it bites your dog, its not my fault. If
; it causes you to curl up on the floor in the foetal position muttering
; about pixies and mushrooms, its not my fault. If it causes you or someone
; else to spontaneously burst into song and dance, its not my fault but
; I'd like to hear about it. You have been warned.
;

;convert a value, val, that is num bytes long, into 
;a series of ordered bytes
function mas_dicom_writer_getbytes, val, num

    ret = BYTARR(num)
    offset = 0

    ;;work in big endian ONLY
    ;;val=swap_endian(val)

    little_endian = (BYTE(1, 0, 1))[0]

    if (little_endian) then begin
       byteorder,val,/SWAP_IF_BIG_ENDIAN 
    endif else begin
       byteorder,val,/SWAP_IF_LITTLE_ENDIAN 
    endelse

    for i = 0,(num-1) do begin
       tmpres = BYTE(ISHFT(val, offset) AND 255)
       ret[i] = tmpres
       offset = offset-8
    endfor
    
    return, ret
end

;;generate any tag
function mas_dicom_writer_gen_anytag, group, element, data, STR=str

    pad = BYTE(0)
    
    ;;check to see if string type is set - if it is, change 
    ;;padding byte to a space
    if KEYWORD_SET(STR) then begin 
       pad = BYTE(STRING(' '))
    endif
    
    rs = [ mas_dicom_writer_getbytes(group,2), $
           mas_dicom_writer_getbytes(element,2) ]
 
    ;;correct to even length if necessary
    bs = BYTE(data)
    nl = n_elements(bs)
    if ((nl mod 2) ne 0) then begin
       bs = [bs,pad]
       nl = nl+1
    end

    ;;size of field
    rs = [rs,mas_dicom_writer_getbytes(nl,2)]

    ;;padding
    rs = [rs,[0,0]]

    ;;string itself
    rs = [rs,bs]

    return, rs

end

; generate a tag based on its data and VR (value representation)
; based on DICOM specs 3.5-1999 (table 6.2-1)
;
; usage: mas_dicom_writer_gen_VRtag(group, element, 'XX', data)
;        where XX is one of the supported VR types below
;
; This is a list of the current VR types supported/not supported
; and the expected data type for the 'data' variable
;
; AE:
; * Application Entity - normal string tag
; * 16 bytes max
; * STRING
;
; AS:
; * Age String - should be nnnX, where X={D,W,M,Y} (days, weeks, months, years),
; * 4 bytes fixed
; * STRING
;
; AT:
; * Attribute tag - should be a pair of unsigned integers representing a data
;   element tag eg. ['0018'x,'00FF'x]
; * 8 bytes fixed
; * [UINT,UINT]
;
; CS:
; * Code string
; * 32 bytes maximum
; * STRING
;
; DA:
; * Date string - 8 bytes fixed, formay yyyymmdd, or 10 bytes fixed
;   yyyy.mm.dd, which is compatible with versions prior dicom v3.0 -
;   so thats what will be used
; * 10 bytes fixed
; * STRING
;
; DS
; * Decimal string - convert an float into a string, and store
; * 16 bytes maximum
; * FLOAT/DOUBLE
;
; DT
; * Date/time string - 26 byte maximum string of format:
;   YYYMMDDGGMMSS.FFFFFF
; * 26 bytes max
; * STRING
; 
; FL:
; * Floating point single - 4 byte fp single val
; * storing as LITTLE ENDIAN - needs to be checked!!!!
; * 4 bytes fixed
; * FLOAT
; 
; FD:
; * Floating point double - 8 byte fp double val
; * storing as LITTLE ENDIAN - needs to be checked!!!!
; * 8 bytes fixed
; * DOUBLE 
;
; IS:
; * Decimal string - convert an int into a string, and store
; * 12 bytes maximum
; * FIX
; 
; LO:
; * long string - IDL doesn't care about this one
; * 64 bytes maximum
; * LONG
;
; LT:
; * long text - IDL doesn't care about this one too much
; * 10240 bytes maximum
; * STRING
; 
; OB
; * other byte string - padded by 00H
; * length variable
; * STRING/BYTE
; 
; OW
; * other word string - padded by 00H. not sure if this is working
; * length variable
; * STRING/BYTE
;
; PN:
; * person name - not supported! (yet?)
; 
; SH
; * short string 
; * 16 bytes maximum
; * STRING
; 
; SL:
; * signed long int
; * 4 bytes fixed
; * LONG
;
; SQ:
; * sequence of items - not supported!
; 
; SS
; * signed short
; * 2 bytes fixed
; * FIX
; 
; ST:
; * short text 
; * 1024 bytes maximum
; * STRING
; 
; TM:
; * time - of format hhmmss.frac
; * 16 bytes maximum
; * STRING
;
; UI:
; * unique identifier
; * 64 bytes maximum
; * STRING
; 
; UL:
; * unsigned long
; * 4 bytes fixed
; * ULONG
; 
; UN:
; * unknown - do whatever you please with this one
; * variable length
; * STRING/BYTE
;
; US:
; * unsigned short
; * 2 bytes fixed
; * UINT
; 
; UT:
; * unlimited text; could be huge!
; * variable length
; * STRING
;
function mas_dicom_writer_gen_VRtag, group, element, VR, data
	
    CASE VR of
       
       'AE': begin
          ;;Application Entity - normal string tag, truncated to 16 bytes
          dval = STRMID(data,0,16)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end

       'AS': begin
          ;;Age String - should be nnnX, where X={D,W,M,Y} (days, weeks, months, years),
          ;;truncated to 4 bytes
          dval = STRMID(data,0,4)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end

       'AT': begin
          ;;Attribute tag - should be a pair of unsigned integers representing a data
          ;;element tag eg. ['0018'x,'00FF'x]
          dval = [mas_dicom_writer_getbytes(UINT(data[0]),2), $
                  mas_dicom_writer_getbytes(UINT(data[1]),2)]
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end

       'CS': begin
          ;;Code string - 32 byte string
          dval = STRMID(data,0,32)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end

       'DA': begin
          ;;Date string - 8 bytes fixed, formay yyyymmdd, or 10 bytes fixed
          ;;yyyy.mm.dd, which is compatible with versions prior dicom v3.0 -
          ;;so thats what will be used
          dval = STRMID(data,0,10)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end

       'DS': begin
          ;;Decimal string - convert an float into a string, and store
          ;;16 bytes maximum
          dval = STRTRIM(STRING(data),1)
          split = strsplit(dval, '\', /extract)
          nitems = n_elements(split)
          for i = 0, nitems-1 do begin
            split[i] = strtrim(split[i], 2)
            split[i] = strmid(split[i], 0, 16)
          endfor
          dval = strjoin(split, '\')
          ;;dval = STRMID(dval,0,16)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end
       
       'DT': begin
          ;;Date/time string - 26 byte maximum string of format:
          ;;YYYMMDDGGMMSS.FFFFFF
          dval = STRMID(data,0,26)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end

       'FL': begin
          ;;Floating point single - 4 byte fp single val
          ;;storing as LITTLE ENDIAN - needs to be checked!!!!
          dval = FLOAT(data)
          ;;explicily cast to bytes
          dvaltmp = BYTE(dval,0,4)
          ;;fix byteorder (little-endian) if necessary; word length is 16 bits (?is this correct?)
          dval = [ mas_dicom_writer_getbytes(UINT(dvaltmp[0:1],0,1),2), $
                   mas_dicom_writer_getbytes(UINT(dvaltmp[2:3],0,1),2) ]
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'FD': begin
          ;;Floating point double - 8 byte fp double val
          ;;storing as LITTLE ENDIAN - needs to be checked!!!!
          dval = DOUBLE(data)
          ;;explicily cast to bytes
          dvaltmp = BYTE(dval,0,8)
          ;;fix byteorder (little-endian) if necessary; word length is 16 bits (?is this correct?)
          dval = [ mas_dicom_writer_getbytes(UINT(dvaltmp[0:1],0,1),2), $
                   mas_dicom_writer_getbytes(UINT(dvaltmp[2:3],0,1),2), $
                   mas_dicom_writer_getbytes(UINT(dvaltmp[4:5],0,1),2), $
                   mas_dicom_writer_getbytes(UINT(dvaltmp[6:7],0,1),2) ]
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end

       'IS': begin
          ;;Decimal string - convert an int into a string, and store
          ;;12 bytes maximum
          dval = STRTRIM(STRING(data),1)
          dval = STRMID(dval,0,12)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end
       
       'LO': begin
          ;;long string - IDL doesn't care about this one, 64 bytes max
          dval = STRMID(data,0,64)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)			
       end
       
       'LT': begin
          ;;long text - IDL doesn't care about this one too much, 10240 bytes max
          dval = STRMID(data,0,10240)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)			
       end
       
       'OB': begin
          ;;other byte string - padded by 00H
          dval = data
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'OW': begin
          ;;other word string - padded by 00H. not sure if this is working
          dval = data
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'PN': begin
          ;;person name - not supported! (yet?)
          print, 'PN currently unsupported!'
          rs = BYTE(0)
       end
       
       'SH': begin
          ;;short string - 16 bytes max
          dval = STRMID(data,0,16)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end
       
       'SL': begin
          ;;signed long
          dval = mas_dicom_writer_getbytes(LONG(data),4)
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'SQ': begin
          ;;sequence of items - not supported!
          print, 'SQ currently unsupported!'
          rs = BYTE(0)
       end
       
       'SS': begin
          ;;signed short
          dval = mas_dicom_writer_getbytes(FIX(data),2)
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'ST': begin
          ;;short text - 1024 bytes max
          dval = STRMID(data,0,1024)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end
       
       'TM': begin
          ;;time - of format hhmmss.frac, 16 bytes maximum
          dval = STRMID(data,0,16)
          rs = mas_dicom_writer_gen_anytag(group,element,dval,/STR)
       end
       
       'UI': begin
          ;;unique identifier, 64 bytes maximum
          dval = STRMID(data,0,64)
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'UL': begin
          ;;unsigned long
          dval = mas_dicom_writer_getbytes(ULONG(data),4)
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'UN': begin
          ;;unknown - do whatever you please with this one
          dval = data
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'US': begin
          ;;unsigned short
          dval = mas_dicom_writer_getbytes(UINT(data),2)
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
       'UT': begin
          ;;unlimited text; could be huge!
          dval = data
          rs = mas_dicom_writer_gen_anytag(group,element,dval)
       end
       
    ENDCASE
    return, BYTE(rs)
end

;;generate pixel tag
function mas_dicom_writer_gen_pixeltag, group, element, val
    return, BYTE([ mas_dicom_writer_getbytes(group,2), $
                   mas_dicom_writer_getbytes(element,2), $
                   mas_dicom_writer_getbytes(val,4) ])
end

; Subroutine name: mas_export_nifti_event
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Event handler for the nifti export choice window.
;
; Editing Information:

pro mas_export_dicom_event, event

    common scan_data

    name = widget_info(event.id, /uname)
    if (name eq 'nif_cancel') then begin
       widget_control, event.top, /destroy
       return
    endif else if (name ne 'nif_export') then begin
       return
    endif

    ci = project.ci
    uid = string(long(systime(1)), format='(I0)')
    dest = dialog_pickfile(path=project.current_path, /directory, /write)
    if (dest eq '') then return
    dest_dir = dest ;; file_dirname(dest)
    SL = (get_dir_slash())[0]

    btn_include_rescale = widget_info(event.top, find_by_uname='btn_include_rescale')
    include_rescale = widget_info(btn_include_rescale, /button_set)
    if (include_rescale eq 0) then begin
        noscale=1
    endif else begin
        noscale=0
    endelse
    
    id = widget_info(event.top, find_by_uname='nif_data')
    if (widget_info(id, /button_set)) then begin
       data_ptr = project.dataarray[ci].state1
       if (ptr_valid(data_ptr)) then begin

          output_ok = 'Yes'

          if (file_test(dest_dir+SL+'SCAN_DATA', /directory)) then begin
              void = dialog_message(['A "SCAN_DATA" directory already exists.', $
                                        'Please remove it and its contents and ', $
                                        'try this export option again.'], /error, /center)
              output_ok = 'No'
          endif

          if (output_ok eq 'Yes') then begin          
              print, "Exporting DICOM: Scan Data"
              file_mkdir, dest_dir+SL+'SCAN_DATA'
              mas_export_dicom, data_ptr=data_ptr, destdir=dest_dir+SL+'SCAN_DATA', uid=uid, noscale=noscale
          endif
          
       endif
    endif

    id = widget_info(event.top, find_by_uname='nif_ad')
    if (widget_info(id, /button_set)) then begin
       data_ptr = project.dataarray[ci].avg_dif
       if (ptr_valid(data_ptr)) then begin

          output_ok = 'Yes'

          if (file_test(dest_dir+SL+'AD', /directory)) then begin
              void = dialog_message(['An "AD" directory already exists.', $
                                        'Please remove it and its contents and ', $
                                        'try this export option again.'], /error, /center)
              output_ok = 'No'
          endif

          if (output_ok eq 'Yes') then begin          
              print, "Exporting DICOM: AD"
              file_mkdir, dest_dir+SL+'AD'
              mas_export_dicom, data_ptr=data_ptr, destdir=dest_dir+SL+'AD', uid=uid, noscale=noscale
          endif
          
       endif
    endif

    id = widget_info(event.top, find_by_uname='nif_fa')
    if (widget_info(id, /button_set)) then begin
       data_ptr = project.dataarray[ci].frac_ani
       if (ptr_valid(data_ptr)) then begin

          output_ok = 'Yes'

          if (file_test(dest_dir+SL+'FA', /directory)) then begin
              void = dialog_message(['An "FA" directory already exists.', $
                                        'Please remove it and its contents and ', $
                                        'try this export option again.'], /error, /center)
              output_ok = 'No'
          endif

          if (output_ok eq 'Yes') then begin          
              print, "Exporting DICOM: FA"
              file_mkdir, dest_dir+SL+'FA'
              mas_export_dicom, data_ptr=data_ptr, destdir=dest_dir+SL+'FA', uid=uid, noscale=noscale
          endif
          
       endif
    endif

    id = widget_info(event.top, find_by_uname='nif_s0')
    if (widget_info(id, /button_set)) then begin
       data_ptr = ptr_new(reform((*project.dataarray[ci].adt)[*,*,*,0]))
       
       if (ptr_valid(data_ptr)) then begin

          output_ok = 'Yes'

          if (file_test(dest_dir+SL+'S0', /directory)) then begin
              void = dialog_message(['An "S0" directory already exists.', $
                                        'Please remove it and its contents and ', $
                                        'try this export option again.'], /error, /center)
              output_ok = 'No'
          endif

          if (output_ok eq 'Yes') then begin          
              print, "Exporting DICOM: S0"
              file_mkdir, dest_dir+SL+'S0'
              mas_export_dicom, data_ptr=data_ptr, destdir=dest_dir+SL+'S0', uid=uid, noscale=noscale
          endif
          
       endif
       
       ptr_free, data_ptr
       
    endif

    widget_control, event.top, /destroy

end

pro mas_dicom_writer, filename, image, VOXELSIZE=voxelsize, SSAI=ssai, $
                      PATIENT=patient, PHYSICIAN=physician, PATID=patid, $
                      sliceloc=sliceloc, TOTAL_IMAGES=TOTAL_IMAGES, $
                      BVEC=bvec, BVAL=bval, uid=uid, $
                      rsslope=rsslope, rsintercept=rsintercept 

    common scan_data
    ci = project.ci
    ;;print, '-------------------------------------------'
    ;;print, 'IDL DICOM writer    by   Bhautik Joshi 2002'
    ;;print, 'http://cow.mooh.org   bjoshi@geocities.com '
    ;;print, '-----------------v 0.2---------------------'

    ;;determine type of image and set bpp
    
    bpp = 0
    
    imtype = size(image,/type)
    
    case imtype of
       1: bpp = 1               ;;byte
       2: bpp = 2               ;int
       3: bpp = 4               ;long int
       12: bpp = 2              ;unsigned int
       13: bpp = 4              ;unsigned long int
       else: bpp = 0
    endcase
    
    if (bpp eq 0) then begin
       print, 'Only integer type images (byte/int/long) supported at this time, sorry!'
       return
    endif
    
 ;;   if (bpp eq 1) then print, 'Byte type (bpp=1) image'
 ;;   if (bpp eq 2) then print, 'Integer type (bpp=2) image'
 ;;   if (bpp eq 4) then print, 'Long type (bpp=4) image'
    
    byteswapped_image = image
    
    if (bpp ge 2) then begin
       little_endian  =  (BYTE(1, 0, 1))[0]
       if (little_endian) then begin
;;          print, 'Little endian architecture detected'
;;          byteorder,byteswapped_image,/SWAP_IF_LITTLE_ENDIAN 
       endif else begin
           print, 'Big endian architecture detected'
          byteorder,byteswapped_image,/SWAP_IF_BIG_ENDIAN 
       endelse
    endif
    
    ;;dummy fill-in variables. 
    if keyword_set(uid) then begin
        study_uid = string(uid)
    endif else begin
        study_uid = '123456'
    endelse
        
    scan_datetime = stregex(project.imndarray[0].scan_date, $
                             '^([a-zA-Z]{3})[ ]+([0-9]+)[ ]+([0-9]{4})[ ]([0-9]{2}):([0-9]{2}):([0-9]{2})$', $
                             /extract, /subexpr)
    case scan_datetime[1] of
        'Jan': mon = '01'
        'Feb': mon = '02'
        'Mar': mon = '03'
        'Apr': mon = '04'
        'May': mon = '05'
        'Jun': mon = '06'
        'Jul': mon = '07'
        'Aug': mon = '08'
        'Sep': mon = '09'
        'Oct': mon = '10'
        'Nov': mon = '11'
        'Dec': mon = '12'
        else: begin
            print, "mas_dicom_writer: unknown month: "+scan_datetime[1]
            mon = '01'
        end
    endcase
    
    scan_date = scan_datetime[3]+mon+scan_datetime[2]
    scan_time = scan_datetime[4]+scan_datetime[5]+scan_datetime[6]
    
    ;;image variables
    sz = size(image)
    cols = sz[1]
    rows = sz[2]  
    ;;parameter set variables
    
    ;;voxelsize=[x,y,z]
    if (KEYWORD_SET(voxelsize)) then begin
       thickness = string(voxelsize[2], format='(f0.4)')
       spacing = string(voxelsize[0], format='(f0.4)')+'\'+string(voxelsize[1],format='(f0.4)')
    endif else begin
       thickness = '1.0'
       spacing = '1.0\1.0'
    endelse
    
    if (KEYWORD_SET(total_images)) then begin
       totalNumImages = total_images
    endif else begin
       totalNumImages = 1
    endelse

    if (KEYWORD_SET(ssai)) then begin
       StudyID = ssai[0]
       Seriesnum = ssai[1]
       Acqnum = ssai[2]
       Imagenum = ssai[3]
    endif else begin
       StudyID = 0
       Seriesnum = 0
       Acqnum = 0
       Imagenum = 0
    endelse

    if (n_elements(rsintercept) eq 0) then begin
        rsintercept = 0.0
    endif    
    if (n_elements(rsslope) eq 0) then begin
        rsslope = 1.0
    endif
    
    window_width = 32767L * rsslope + rsintercept
    window_center = window_width/2.0
    
    if (KEYWORD_SET(patient)) then patient = STRING(patient) else patient = 'Unnamed'
    if (KEYWORD_SET(physician)) then physician = STRING(physician) else physician = 'Unnamed'
    if (KEYWORD_SET(patid)) then patid = STRING(patid) else patid = 'PATXXX'
       
    ;;SOP class set to MR - see Annex A in PS 3.6-2001
    SOPClass = '1.2.840.10008.5.1.4.1.1.4' ;; 0008/0016/UI
    SOPInstancePrefix = '1.2.9999.10101.733.23.3' ;; not dicom
    
    StudyInstanceUID       = SOPInstancePrefix + study_uid ;; 0020/000D/UI
    SeriesInstanceUID      = StudyInstanceUID + '.'+string(Acqnum, format='(G0)') ;; 0020/000E/UI
    RelFrameOfReferenceUID = StudyInstanceUID + '.2' ;; 0020/0052/UI
    
    ;; this is unique per image!
    instance_id = string(imagenum, format='(I0)')
    SOPInstance = SeriesInstanceUID+'.'+instance_id;; 0008/0018/UI
    
    GET_LUN, U
    
    OPENW, U, filename
    
;;    print, 'Writing to file ', filename
    ;; Transfer syntax uid
    ;;WRITEU, U, mas_dicom_writer_gen_VRtag('0002'x,'0010'x,'UI','1.2.840.10008.1.2')
    ;;MR type	
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0008'x,'CS','ORIGINAL\PRIMARY\OTHER')
    ;;Instance date 
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0012'x,'DA',scan_date)
    ;;Instance time
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0013'x,'TM',scan_time)
    ;;SOP class
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0016'x,'UI',SOPClass)
    ;;SOP instance
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0018'x,'UI',SOPInstance)
    ;;Study Date
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0020'x,'DA',scan_date)
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0021'x,'DA',scan_date)
    ;;Acquisition Date
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0022'x,'DA',scan_date)
    ;;Study Time
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0030'x,'TM',scan_time)
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0031'x,'TM',scan_time)
    
    ;;Acquisition Time
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0032'x,'TM',scan_time)
    ;;Modality
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0060'x,'CS','MR')
    ;;Manufacturer
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0070'x,'LO','MAS DICOM EXPORTER')
    ;;Study Physicians Name
    WRITEU, U, mas_dicom_writer_gen_VRtag('0008'x,'0090'x,'UI',physician)
    
    ;; 0010 tags
    
    ;;Patient name
    WRITEU, U, mas_dicom_writer_gen_VRtag('0010'x,'0010'x,'UI',patient)
    ;;Patient ID
    WRITEU, U, mas_dicom_writer_gen_VRtag('0010'x,'0020'x,'LO',patid)
    ;;Patient birth date
    WRITEU, U, mas_dicom_writer_gen_VRtag('0010'x,'0030'x,'DA','20020114')
    ;;Patient sex
    WRITEU, U, mas_dicom_writer_gen_VRtag('0010'x,'0040'x,'CS','M')
    ;;Patient Weight
    WRITEU, U, mas_dicom_writer_gen_VRtag('0010'x,'1030'x,'DS','100')

    ;; 0018 tags
    ;;Scanning Sequence
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0020'x,'CS','RM')
    ;;Scanning Variant
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0021'x,'CS','NONE')
    ;;Scanning Options
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0022'x,'CS','')
    ;;Acquisition type
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0023'x,'CS','2D')
    ;;Sequence Name
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0024'x,'SH',project.imndarray[ci].scan_name)
    ;;Slice thickness
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0050'x,'DS',thickness)
    ;;Repetition Time
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0080'x,'DS',project.imndarray[ci].recov_time)
    ;;Echo Time
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0081'x,'DS',project.imndarray[ci].echo_time)
    ;;Number of Averages
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0083'x,'DS',fix(project.imndarray[ci].n_avg))
    ;;Slice spacing
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0088'x,'DS',thickness)
    ;;Echo Train Length
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'0091'x,'IS',1)
    ;;Protocol Name
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'1030'x,'LO',project.imndarray[ci].scan_name)
    ;;Patient Position
    WRITEU, U, mas_dicom_writer_gen_VRtag('0018'x,'5100'x,'CS','HFS')
    
    ;; 0020 tags
    
    ;;Study instance
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'000D'x,'UI',StudyInstanceUID)
    ;;Series instance UID
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'000E'x,'UI',SeriesInstanceUID)
    ;;StudyID
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'0010'x,'IS',StudyID)
    ;;Series number
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'0011'x,'IS',Acqnum)
    ;;Acquisition number
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'0012'x,'IS',acqnum)
    ;;Image number
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'0013'x,'IS',imagenum)
    ;;Image Position
    acq_mat = reform(float((*project.imndarray[ci].acq_matrix)[2,0:2])) * sliceloc
    impos = strcompress(string(acq_mat[0], format='(f0.5)')+'\'+ $
                          string(acq_mat[1], format='(f0.5)')+'\'+ $
                          string(acq_mat[2], format='(f0.5)'), /remove_all)
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'0032'x,'DS',impos)
    ;;Image Orientation
    acq_mat = fix((*project.imndarray[ci].acq_matrix)[0:2,0:2])
    img_orient = strcompress(string(acq_mat[0,0])+'\'+ $
                               string(acq_mat[0,1])+'\'+ $
                               string(acq_mat[0,2])+'\'+ $
                               string(acq_mat[1,0])+'\'+ $
                               string(acq_mat[1,1])+'\'+ $
                               string(acq_mat[1,2]), /remove_all)
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'0037'x,'DS',img_orient)
    ;;Frame of Reference UID
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'0052'x,'UI',RelFrameOfReferenceUID)
    ;;Total # of images in acquisition
    WRITEU, U, mas_dicom_writer_gen_VRtag('0020'x,'1002'x,'IS',totalNumImages)
    ;;Slice Location
    WRITEU, u, mas_dicom_writer_gen_VRtag('0020'x,'1041'x,'DS',string(sliceloc, format='(f0.5)'))
    
    ;; 0028 tags
    
    ;;samples per pixel
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0002'x,'US',1)
    ;;Photometric interpretation
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0004'x,'CS','MONOCHROME2')
    ;;Rows in image
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0010'x,'US',rows)
    ;;Columns in image
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0011'x,'US',cols)
    ;;pixel spacing
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0030'x,'DS',spacing)
    ;;bits allocated per sample
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0100'x,'US',bpp*8)
    ;;bits stored per sample
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0101'x,'US',bpp*8)
    ;;high bit
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0102'x,'US',(bpp*8)-1)
    ;;pixel representation
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0103'x,'US','0001'x)
    ;;min pixel value
;    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0106'x,'SS',FIX(0))
    ;;max pixel value
;    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'0107'x,'SS',FIX(32767+10))
    ;;window center
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'1050'x,'DS',window_center)
    ;;window width
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'1051'x,'DS',window_width)
    ;;rescale values
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'1052'x,'DS',rsintercept)
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'1053'x,'DS',rsslope)
    
    ;; window algorithm
    WRITEU, U, mas_dicom_writer_gen_VRtag('0028'x,'1055'x,'LO','MinMax')

    WRITEU, U, mas_dicom_writer_gen_VRtag('2001'x,'0010'x,'LO','MAS_DICOM_WRITER')
    if (n_elements(bval) ne 0) then begin
        WRITEU, U, mas_dicom_writer_gen_VRtag('2001'x,'1003'x,'DS',bval)
    endif 
    
    WRITEU, U, mas_dicom_writer_gen_VRtag('2005'x,'0010'x,'LO','MAS_DICOM_WRITER')
    if (n_elements(bvec) eq 3) then begin
        WRITEU, U, mas_dicom_writer_gen_VRtag('2005'x,'10B0'x,'DS',bvec[0])
        WRITEU, U, mas_dicom_writer_gen_VRtag('2005'x,'10B1'x,'DS',bvec[1])
        WRITEU, U, mas_dicom_writer_gen_VRtag('2005'x,'10B2'x,'DS',bvec[2])
    endif
    
    ;;write image data
    imsize = rows*cols*bpp
    
    WRITEU, U, mas_dicom_writer_gen_pixeltag('7FE0'x,'0010'x,imsize)
    
    WRITEU, U, byteswapped_image
    
    CLOSE, U
    FREE_LUN, U
    
;;    print, 'Successfully wrote ', filename
end

pro mas_export_dicom_gui

    mas_export_nifti_gui, /DICOM

end


pro mas_export_dicom, data_ptr=data_ptr, destdir=destdir, uid=uid, noscale=noscale

    common scan_data
    forward_function get_dir_slash
    
    ci = project.ci
    
    patient = 'Unnamed Patient'
    
    ;; get the destination dir if not provided
    if (not keyword_set(destdir)) then begin
        destdir = dialog_pickfile(title="Select a destination directory...", $
            path='~/Desktop', $ ;;project.current_path, $
            /directory, /write)
        if (destdir eq '') then return
        
    endif else if (destdir eq '') then return
    
    ;; use state1 unless provided with other data
    if (n_elements(data_ptr) eq 0) then begin
        data_ptr = project.dataarray[ci].state1
    endif
    
    if (n_elements(uid) eq 0) then begin
        uid = string(long(systime(1)), format='(I0)')
    endif
    
    ;; voxel dimensions, converted to mm
    vox_dims = [ project.imndarray[ci].f_voxsz, $
                 project.imndarray[ci].p_voxsz, $
                 project.imndarray[ci].s_voxsz ] * 10.

    sz_data = size(*data_ptr, /structure)
    
    n_slices     = sz_data.dimensions[2] eq 0 ? 1 : sz_data.dimensions[2]
    n_volumes    = sz_data.dimensions[3] eq 0 ? 1 : sz_data.dimensions[3]
    total_images = n_volumes * n_slices
    
    dir_slash = (get_dir_slash())[1]
    
    progressbar = OBJ_NEW('progressbar', Color='red', /nocancel, $
        Title='Writing DICOM sequence...')
    progressbar->Start
    
    for vol = 0, n_volumes-1 do begin
    
        ;; convert floating point to "short" and scale accordingly
        rsslope = 1.0
        if (sz_data.type eq 4 or sz_data.type eq 5) then begin
            volume_n = reform((*data_ptr)[*,*,*,vol])
            max_vol = float(max(volume_n))
            volume_n = fix(volume_n/max_vol * 32767., type=2)
            rsslope = max_vol/32767.0
        endif
        
        if (keyword_set(noscale)) then rsslope = 1.0
        
        if (ptr_valid(project.imndarray[ci].bval_array)) then begin
            
            bval  = (*project.imndarray[ci].bval_array)[vol]
            
            if (bval gt 0.0) then begin
                theta = (*project.imndarray[ci].angle_theta)[vol]
                phi   = (*project.imndarray[ci].angle_phi)[vol]            
                bvec = cv_coord(from_sphere=[phi, 90-theta, 1.], /to_rect, /degrees)
            endif else begin
                bvec = float([0,0,0])
            endelse
            
        endif
        
        ;; write out each slice as a single DICOM image
        for sl = 0, n_slices-1 do begin
            seq_num = vol*n_slices + sl
            ssai = [1,1,vol+1,sl+1]
            file_name = destdir + dir_slash +'IM_'+string(seq_num, $
                format='(I04)')+'.dcm'
                
            mas_dicom_writer, file_name, $
                reverse(reform(volume_n[*,*,sl]),2), ssai=ssai, $
                sliceloc=(sl+1)*vox_dims[2], $
                voxelsize=vox_dims, $
                total_images=total_images, $
                bvec=bvec, bval=bval, $
                patient=patient, uid=uid, rsslope=rsslope
                
            if (seq_num mod 100 eq 0) then begin
                progressbar->update, 100. * (float(seq_num)/float((total_images))), $
                    TEXT="Writing file: "+strcompress(string(seq_num)+' of '+string(total_images))
            endif
            
        endfor
        
    endfor
    
    progressbar->destroy
    
end

pro mas_dicom_element__define

    struct = { MAS_DICOM_ELEMENT, $
               group:     '0000'x, $
               element:   '0000'x, $
               attr_name: 'Attribute Name', $
               VR:        'VR', $
               data_in:   ptr_new(), $
               data_enc:  ptr_new(), $
               comment:   '' }

end

;pro mas_dicom_export_GUI
;
;    common common_widgets
;    common scan_data
;
;    base = widget_base(group_leader=WID_BASE_MAIN, /column, /modal)
;    
;    tab_base = widget_tab(base)
;
;    pat_tab = widget_base(tab_base, title="Patient", /column)
;    stu_tab = widget_base(tab_base, title="Study", /column)
;    ser_tab = widget_base(tab_base, title="Series", /column)
;    eqp_tab = widget_base(tab_base, title="Equipment")
;    img_tab = widget_base(tab_base, title="Image")
;    mri_tab = widget_base(tab_base, title="MR")
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    void_base = widget_base(pat_tab, /row)
;    lbl_pat_name = widget_label(void_base, value="Patient name:")
;    txt_pat_name = widget_text(void_base, xsize=30, value='Unnamed Patient', /editable)
;
;    void_base = widget_base(pat_tab, /row)
;    lbl_pat_id = widget_label(void_base, value="Patient ID:")
;    txt_pat_id = widget_text(void_base, xsize=10, value='0001', /editable)
;    
;    void_base = widget_base(pat_tab, /row)
;    lbl_pat_dob = widget_label(void_base, value="Patient Birth Date:")
;    txt_pat_dob = widget_texT(void_base, xsize=10, value='01/01/2000', /editable)
;    lbl_void = widget_label(void_base, value='(MM/DD/YYYY)')
;
;    void_base = widget_base(pat_tab, /row)
;    lbl_pat_sex = widget_label(void_base, value="Patient Sex:")
;    dl_pat_sex = widget_droplist(void_base, value=['Male','Female'])
;
;    void_base = widget_base(pat_tab, /row)
;    lbl_pat_weight = widget_label(void_base, value="Patient Weight:")
;    txt_pat_weight = widget_text(void_base, xsize=10, value='0', /editable)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    void_base = widget_base(stu_tab, /row)
;    lbl_stu_instanceUID = widget_label(void_base, value="Study Instance UID:")
;    txt_stu_instanceUID = widget_text(void_base, xsize=30, value='<<>>')
;    
;    void_base = widget_base(stu_tab, /row)
;    lbl_stu_id = widget_label(void_base, value="Study ID:")
;    txt_stu_id = widget_text(void_base, xsize=30, value='1', /editable)
;    
;    void_base = widget_base(stu_tab, /row)
;    lbl_stu_date = widget_label(void_base, value="Study Date:")
;    txt_stu_date = widget_text(void_base, xsize=10, value='01/01/2000', /editable)
;    lbl_void = widget_label(void_base, value='(MM/DD/YYYY)')
;    
;    void_base = widget_base(stu_tab, /row)
;    lbl_stu_time = widget_label(void_base, value="Study Time:")
;    txt_stu_time = widget_text(void_base, xsize=10, value='00:00:00', /editable)
;    lbl_void = widget_label(void_base, value='(HH:MM:SS)')
;    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    void_base = widget_base(_tab, /row)
;    lbl_stu_instanceUID = widget_label(void_base, value="Study Instance UID:")
;    txt_stu_instanceUID = widget_text(void_base, xsize=30, value='<<>>')
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;    widget_control, base, /realize
;   
;    
;
; end
