;; $Id$
;;
FUNCTION extract_echo_label, echotxt;

	;type      0     1     2     3     4     5     99
	types = ['STD','REJ','PHX','FRX','NOI','NAV','END OF DATA VECTOR INDEX']
	; -1 = no valid type

	;# typ mix   dyn   card  echo  loca  chan  extr1 extr2 ky    kz    aver
	;  sign  rf    grad  enc   rtop  rr    size   offset
	result = LONARR(20);

	echo_type = STRMID(echotxt,2,3)
	type = -1

	CASE echo_type OF
		types[0]: type = 0
		types[1]: type = 1
		types[2]: type = 2
		types[3]: type = 3
		types[4]: type = 4
		types[5]: type = 5
		ELSE    : BEGIN
			IF STRPOS(echotxt, types[6]) NE -1 THEN type = 99 $
			ELSE print,'Error: Echo type unknown.'
		END
	ENDCASE

	IF type EQ -1 THEN print, echotxt


    result(0) = type
    IF (type GE 0) AND (type LE 5) THEN BEGIN
        pos = STRPOS (echotxt,types(type))
        txt = STRTRIM(STRMID(echotxt,pos+3),1)+' '

		FOR i = 1, (SIZE(result))[1] - 1 DO BEGIN
           pos = STRPOS(txt,' ')
           num = STRMID(txt,0,pos)
           result(i) = LONG(num)
           txt = STRTRIM(STRMID(txt,pos+1),1)
        ENDFOR
    ENDIF

	RETURN, result
END


FUNCTION sequential_header_reader, fid, search_string, NUMREAD=numread
	numread = 0
	intxt = ''

	REPEAT BEGIN

    	READF, fid, intxt
    	numread = numread + 1
    	pos = STRPOS(intxt, search_string)

    	IF pos NE -1 THEN RETURN, intxt

	ENDREP UNTIL EOF(fid)
END


FUNCTION process_DATALIST_files, path_pattern

	list_file = path_pattern + '.LIST'
	data_file = path_pattern + '.DATA'


	OPENR, lfid, list_file, /GET_LUN
	OPENR, dfid, data_file, /GET_LUN

	labcnt = 0L;

	;**************************
	;* Get Header information *
	;**************************

	;str_kxo = 'kx_oversample_factor'
	;str_kyo = 'ky_oversample_factor'
	;str_xre = 'X-resolution'
	;str_yre = 'Y-resolution'
	;str_xra = 'X_range'
	;str_yra = 'Y_range'
	;str_end = 'END OF DATA VECTOR INDEX'

	scanname = sequential_header_reader(lfid,'Scan name',NUMREAD=numread)
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number_of_mixes',NUMREAD=numread)
	mix = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number_of_encoding_dimensions',NUMREAD=numread)
	enc = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number_of_dynamic_scans',NUMREAD=numread)
	dyn = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number_of_cardiac_phases',NUMREAD=numread)
	cap = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number_of_echoes',NUMREAD=numread)
	ech = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number_of_locations',NUMREAD=numread)
	loc = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number_of_signal_averages',NUMREAD=numread)
	nsa = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	txt = sequential_header_reader(lfid,'number of coil channels',NUMREAD=numread)
	coc = FIX(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)))
	labcnt = labcnt + numread

	kxr = INTARR(3)
	txt = sequential_header_reader(lfid,'kx_range',NUMREAD=numread)
	ktxt= STRTRIM(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)),1)
	pspa = STRPOS(ktxt,' ')
	num = STRMID(ktxt,0,pspa)
	kxr(0)= FIX(num)
	ktxt = STRTRIM(STRMID(ktxt,pspa+1),2)
	kxr(1) = FIX(ktxt)
	kxr(2)=kxr(1)-kxr(0)+1
	labcnt = labcnt + numread

	kyr=INTARR(3);
	txt = sequential_header_reader(lfid,'ky_range',NUMREAD=numread)
	ktxt= STRTRIM(STRMID(txt,STRTRIM(STRPOS(txt,':')+1)),1)
	pspa = STRPOS(ktxt,' ')
	num = STRMID(ktxt,0,pspa)
	kyr(0)= FIX(num)
	ktxt = STRTRIM(STRMID(ktxt,pspa+1),2)
	kyr(1) = FIX(ktxt)
	kyr(2)=kyr(1)-kyr(0)+1
	labcnt = labcnt + numread

	;******************************
	;* MOVE TO START OF DATA INFO *
	;******************************

	txt = sequential_header_reader(lfid,'START OF DATA VECTOR INDEX ',NUMREAD=numread)
	labcnt = labcnt + numread

	; skip four table header lines
	READF, lfid, txt
	READF, lfid, txt
	READF, lfid, txt
	READF, lfid, txt

	;*******************
	;* Prepare k-space *
	;*******************

	num_chn = FIX(coc)

	; channels * kx * ky * slice * echos
	kselv = COMPLEXARR(num_chn,kxr(2),kyr(2),FIX(loc),FIX(ech))
	noise_found = 0

	;**************
	;* Get Echoes *
	;**************
	endvector = 0;
	coil_no=0;

	REPEAT BEGIN

		READF, lfid, txt

		echo_header = extract_echo_label(txt)

		; check echo type
		CASE echo_header(0) OF
			0: BEGIN ;STD, image data
				echo_data = COMPLEXARR(echo_header(19)/8)
				READU, dfid, echo_data

				kselv[coil_no,*, echo_header(9)-kyr(0),echo_header(5),echo_header(4)] = echo_data
				coil_no = (coil_no + 1) MOD num_chn
			END

			1: BEGIN
				ignore_data = COMPLEXARR(echo_header(19)/8)
				READU, dfid, ignore_data
			END

			2: BEGIN
				ignore_data = COMPLEXARR(echo_header(19)/8)
				READU, dfid, ignore_data
			END

			3: BEGIN
				ignore_data = COMPLEXARR(echo_header(19)/8)
				READU, dfid, ignore_data
			END

			4: BEGIN ;NOI, channel noise
				IF noise_found EQ 0 THEN BEGIN
					noise_array = COMPLEXARR(num_chn, echo_header(19)/8)
					noise_found = 1
				ENDIF

				noise_data = COMPLEXARR(echo_header(19)/8)
				READU, dfid, noise_data

				noise_array[coil_no,*] = noise_data
				coil_no = (coil_no + 1) MOD num_chn
			END

			5: BEGIN
				ignore_data = COMPLEXARR(echo_header(19)/8)
				READU, dfid, ignore_data
			END

			99: endvector = 1

			ELSE:

		ENDCASE

		labcnt = labcnt + 1
	ENDREP UNTIL (endvector NE 0) OR EOF(lfid)

	CLOSE, lfid
	CLOSE, dfid

	data = { $
		scanname:scanname, $
		mix:mix, $
		enc:enc, $
		dyn:dyn, $
		cap:cap, $
		ech:ech, $
		loc:loc, $
		nsa:nsa, $
		coc:coc, $
		kxr:kxr, $
		kyr:kyr, $
		num_chn:num_chn, $
		kselv:PTR_NEW(kselv), $
		noise_array:PTR_NEW(noise_array) $
	}

	RETURN, data
END


FUNCTION create_receiver_noise_matrix, noise_array
	i = COMPLEX(0,1)
	sz = SIZE(noise_array)

	no_chnl = sz[1]
	no_samp = sz[2]

	noise_mat = COMPLEXARR(no_chnl,no_chnl)

	FOR idx1 = 0, no_chnl-1 DO BEGIN
		FOR idx2 = 0, no_chnl-1 DO BEGIN
			noise_mat[idx1,idx2] = (VARIANCE(noise_array[idx1,*] + noise_array[idx2,*]) + i * $
				VARIANCE(noise_array[idx1,*] - i * noise_array[idx2,*]) - (1 + i) * (VARIANCE(noise_array[idx1,*]) $
				+ VARIANCE(noise_array[idx2,*]))) / 2
		ENDFOR
	ENDFOR

	RETURN, noise_mat
END


PRO READ_PHILIPS_RAW_DATA

	test_dir ="C:\Users\Durgin\Desktop\DATA"
	slice_no = 3 ;for display

	;****************************************************************************************************************

	;*****************************
	;* READ SENSE ENCODED IMAGES *
	;*****************************

	path_pattern = DIALOG_PICKFILE( Title = 'Select either of the DATA/LIST files to open.',path=test_dir)
	IF (path_pattern EQ '') THEN RETURN

	pos = STRPOS(STRMID(path_pattern,0,STRLEN(path_pattern)),'.',/REVERSE_SEARCH)
	path_pattern = STRMID(path_pattern,0,pos)

	kdata = process_DATALIST_files(path_pattern)

	;*****************************
	;* FFT, MAGNITUDE, and SHIFT *
	;*****************************

	shift_amt = kdata.kxr(2)/2;
	chn_img = SHIFT(FFT(FFT(*kdata.kselv,1,DIMENSION=2,/OVERWRITE),1,DIMENSION=3,/OVERWRITE),0,shift_amt,0,0)
	sum_img = TOTAL(ABS(chn_img),1)

	;; adjust for over sampling from philips manual

	;************************
	;* CHANNEL NOISE MATRIX *
	;************************
	ch_noise_mat = create_receiver_noise_matrix(*kdata.noise_array)


	;****************************************************************************************************************

	;****************************
	;* READ SENSITIVITY PROFILE *
	;****************************

	path_pattern = DIALOG_PICKFILE( Title = 'Select either of the DATA/LIST files to open.',path=test_dir)
	IF (path_pattern EQ '') THEN RETURN

	pos = STRPOS(STRMID(path_pattern,0,STRLEN(path_pattern)),'.',/REVERSE_SEARCH)
	path_pattern = STRMID(path_pattern,0,pos)

	sdata = process_DATALIST_files(path_pattern)

;	;**************************
;	;* KSpace Low Pass Filter *
;	;**************************
;
;	fborder = LONG(1*sdata.kxr[2]/8)
;	order = 8
;
;	; low pass butterworth filter
;	filter = SHIFT(1.0 / ( 1.0d + (DIST(sdata.kxr[2],sdata.kyr[2])/fborder)^2*order),sdata.kxr[2]/2,sdata.kyr[2]/2)
;	;Display_DrawResS_Create, bytscl(filter) , 'LP', 0, 0
;
;	; apply filter to all channels and all slices
;	FOR i = 0, sdata.num_chn - 1 DO FOR j = 0, sdata.loc - 1 DO BEGIN ; for each slice
;		(*sdata.kselv)[i,*,*,j] = REFORM((*sdata.kselv)[i,*,*,j]) * filter
;	ENDFOR

	;*****************************
	;* FFT, MAGNITUDE, and SHIFT *
	;*****************************

	shift_amt = sdata.kxr(2)/2;
	sen_chn_img = SHIFT(FFT(FFT(*sdata.kselv,1,DIMENSION=2,/OVERWRITE),1,DIMENSION=3,/OVERWRITE),0,shift_amt,0,0)
	sen_sum_img = TOTAL(ABS(sen_chn_img),1)

	;; adjust for over sampling from philips manual

	; consider thresholding


;****************************************************************************************************************

	;******************
	;* DISPLAY IMAGES *
	;******************

	; display reconstructed SENSE images

	WINDOW, 0, XSIZE=kdata.kxr[2]*2, YSIZE=kdata.kyr[2]*5, XPOS=0, YPOS=0, RETAIN=2

	TVSCL, REFORM(sum_img[*,*,slice_no]), 0

	FOR i = 0, kdata.num_chn-1 DO BEGIN
		scaled_img = ABS(REFORM(chn_img[i,*,*,slice_no])) / REFORM(sum_img[*,*,slice_no])
		TVSCL, BYTSCL(scaled_img), i+1
	ENDFOR

	; display reconstructed sensitivity profile

	WINDOW, 1, XSIZE=sdata.kxr[2]*2, YSIZE=sdata.kyr[2]*5, XPOS=0, YPOS=0, RETAIN=2
	TVSCL, REFORM(sen_sum_img[*,*,slice_no]), 0

	FOR i = 0, sdata.num_chn-1 DO BEGIN
		scaled_img = ABS(REFORM(sen_chn_img[i,*,*,slice_no])) / REFORM(sen_sum_img[*,*,slice_no])
		TVSCL, BYTSCL(scaled_img), i+1
	ENDFOR

;****************************************************************************************************************

	;************************
	;* SENSE RECONSTRUCTION *
	;************************

	;U =





;****************************************************************************************************************


	;********************************
	;* PHILIPS COIL SENSITIVITY MAP *
	;********************************

;	dir = 'C:\Users\Durgin\Desktop\DATA\SENSE TEST\11.13.07 HARDI TEST\SENSE PROFILES\'
;
;	coil1 = dir + 'sri2007111316253268293696i'
;	coil2 = dir + 'sri2007111316255193383696i'
;	coil3 = dir + 'sri2007111316261125473696i'
;	coil4 = dir + 'sri2007111316263070563696i'
;	coil5 = dir + 'sri2007111316265015653696i'
;	coil6 = dir + 'sri2007111316270935743696i'
;	coil7 = dir + 'sri2007111316272870833696i'
;	coil8 = dir + 'sri2007111316274823923696i'
;
;
;
;	rawd = READ_BINARY(coil8,DATA_TYPE=6)
;
;
;	;imgdata = TOTAL(REFORM(TEMPORARY(rawd),32,64,64,15),4)
;	;imgdata = ABS(FFT(FFT(imgdata,1,DIMENSION=1),1,DIMENSION=2))
;
;
;
;	imgdata = TOTAL(REFORM(TEMPORARY(rawd),64,64,32,15),4)
;	imgdata = ABS(FFT(FFT(imgdata,1,DIMENSION=1),1,DIMENSION=2))
;
;
;	Display_DrawResS_Create, BYTSCL(REFORM(imgdata[*,*,15])), '', 0, 0
;
;	;TV, imgdata[*,*,]


END




