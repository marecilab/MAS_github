;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Subroutine name: InvRecovFunction_MP 
; Created by: Bill Triplett
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; parameters:
;   TI         [in] this is the data that needs a curve fit to it. must be same dimensions as indp_var.
;   A          [in] paramters being fit
;
; Returns: Value of inversion recovery expression for TE and A
;
; Editing Information:
;
; MODIFICATION HISTORY:


function InvRecovFunction_MP, TI, A

    common scan_data, project
    
    TR = project.imndarray[project.ci].recov_time * 1e-3
    return,  abs( A[0]*( 1.0 - 2.0*exp(-TI*A[1]) + exp(-TR*A[1]) ) ) 

end

;+
; :Description:
;    Inversion recovery fitting function for use with Powell optimizer.
;
; :Params:
;    P - parameters of the fit
;
; :Returns: 
;    size of difference between model and measurements
;
; :Author: wtriplett
;-
function InvRecovFunction_Powell, P

    common CURVEFIT_POWELL, fit_params
    
    d  = fit_params.measurements
    TI = fit_params.params
    return, norm(d - InvRecovFunction_MP(TI,P))
    
end

function InvRecovFunction_LM, TI, A

    return, [ [ A[0]*(1.0-2.0*exp(-TI*A[1])) ], $
              [ (1.0-2.0*exp(-TI*A[1]))      ], $
              [ TI*(2.0)*A[0]*exp(-TI*A[1])  ] ]

end

;+
; :Description:
;    WIP: Inversion recovery fitting routine. ** WORK IN PROGRESS **
;
; :Params:
;    dep_var  - measured values
;    indp_var - inversion times
;    A        - parameters of the fit
;    chisquare - chisq from fit
;    sigma    - sigma err from fit
;
; :Keywords:
;    noiseerr - not used
;
; :Author: wtriplett
;-
pro InvRecovSE, dep_var, indp_var , A , chisquare, sigma, noiseerr=noiseerr
    
    common CURVEFIT_POWELL, fit_params
    
    catch, cfit_error
    if (cfit_error ne 0) then begin
        catch, /cancel
        message, /info, 'Error detected: '+!ERROR_STATE.MSG
        if (!ERROR_STATE.code eq -637) then return
        return
    endif
    
    ;make an array for the values that we want to converge.
    A=dblarr(2)
    
    ;we need the last value of y - y at x where x is nearest to 0

    A[0] = max(dep_var)
    A[1] = 1.0/mean(indp_var)
    
    ;the parameters we would like to fit. (1=fit, 0=do not fit) that parameter.
    fita = [1,1]

    ;; we are going to use powell optimizer to compute using the abs() of the signal
    ;; this seems to perform better, since we are not guaranteed to have a
    ;; complex signal value to work with, and the abs() of recovery function has a undefined
    ;; derivative at the zero crossing.
    if (1) then begin
    
        fit_params = { measurements: dep_var, params: indp_var }
        Xi = diag_matrix([1.0,1.0])
        dof = n_elements(dep_var)-n_elements(a)
        powell, A, Xi, 1e-7, fmin, 'InvRecovFunction_POWELL'
        
        expected  = invrecovfunction_MP(indp_var, a)
        chisquare = total(((expected - dep_var)^2)/expected)
        sigma     = abs([sqrt(-1), sqrt(-1)]) ;; make these meaningless. 
        chisquare = 1.0
        
    endif else begin
        coefs = MPFITFUN('InvRecovFunction_MM', indp_var, dep_var, err, A, $
                        covar=covariance, perror=sigma, weights=1D, dof=dof, $
                        yfit=yfit, bestnorm=bestnorm, /quiet)
        
        chisquare = sqrt(bestnorm/dof) ;;total((dep_var - yfit)^2)
    endelse
    
 end


; Subroutine name: ProgSatFunct
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; parameters:
;   X     [in]   independent variables
;   A     [in]   values that are being fit
;   return    [out]  the fuction evaluated at all the entries and the partial derivatives WRT each A[]
; description:
;   this function is a support for the progresive saturation curve fitting routine.
;   this routine is the user supplied to fit the curve that describes progressive saturation
;   the function is f(X) = S0 * (1-exp(-t/T1))
;   substute
;     A[0] = S0 which is the intensity
;     A[1] = 1/T1 where T1 is the relaxiation parameter
;     A[2] = Y0 which represents the crossing of the y axis and the RMS of the noise.
;   the final formula is represented as
;     f(X) = A[0](1-exp(-X*A[1]))
;   to do the calculation for curve fitting the function lmfit requires the partial derivatives
;   for all of the A[] matrix.
;     df/dA[0] = 1-exp(-X*A[1])
;     df/dA[1] = A[0]* exp(-X*A[1]) * X
;     df/dA[2] = 1

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

function SatRecovGEFunct, X , A
    RETURN,[ [ A[0]*(1-exp(-X*A[1])) ] , [ 1-exp(-X*A[1]) ] , [ A[0]* exp(-X*A[1]) * X ] ]
END


; Subroutine name: SatRecovGE (formerly ProgSat)
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; parameters:
;   dep_var     [in] this is the data that needs a curve fit to it. must be same dimensions as indp_var.
;   indp_var    [in] independent variable that is increasing.
;   A          [out] array that will contain the fit parameters
;   chisquare  [out] the chisquare for this fit
; description:
;
; Editing Information:
;   Edited by HS 2006/10/09
;   Fix spelling mistakes and commenting.
; MODIFICATION HISTORY:
;   Edited by Alireza V 2006/10/16
;   Changing to the improved signal calculation.
;   Edited by GA 
;   Changed name of function 

pro SatRecovGE , dep_var, indp_var , A , chisquare, sigma
;   print, 'dep_var', dep_var
;   print, 'indp_var', indp_var


    ;make an array for the values that we want to converge.
    A=dblarr(2)
    ;we need the last value of y - y at x where x is nearest to 0
    temp = min(indp_var, Min_Subscript)
    A[0] = max(dep_var)
    ;a[1] is a representation of the slope at the crossing of the y axis
    A[1]= double(dep_var(0)) / double(indp_var(0))



    ;the parameters we would like to fit. (1=fit, 0=do not fit) that parameter.
    fita = [1,1]

    coefs = LMFIT( indp_var , dep_var, A ,FUNCTION_NAME = 'SatRecovGEFunct' $
             ,CONVERGENCE=converg, iter = iterations ,FITA = fita $
             ,chisq = chisquare, covar = covarance, ITMIN=6,SIGMA =sigma )

 end

; Subroutine name: satRecovFunct
; Created by: Garrett Astary
; Calling Information: Called from SatRecovSE which is called from roi_Curve_fit if fit for and ROI is selected, image_fit if fit for an image is selected
;                      In both cases image_type = 1, indentifier for T1 Saturation Recovery for SE to be used (image_fit = 0 for GE)

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; parameters:
;   X     [in]   independent variables
;   A     [in]   values that are being fit
;   return    [out]  the fuction evaluated at all the entries and the partial derivatives WRT each A[]
; description:
;   this function is a support for the saturation recovery for spin echo curve fitting routine.
;   this routine is the user supplied to fit the curve that describes saturation recovery for spin echo
;   Saturation recovery for spin echo allows for longitudinal relaxation between the pi and pi/2 pulses, this effect is more
;   prominent in samples that have a short T1 (samples that contain constrast agents)
;   the function is f(X) = So*(1-2*exp(-(TR-TE/2)*(1/T1))+exp(-TR*(1/T1))
;   substute
;     A[0] = S0 which is the intensity
;     A[1] = 1/T1 where T1 is the relaxiation parameter
;     TE = the echo time in seconds, comes from the T1_fit_data common block of variables 
;     X  is the independent variable (TR) 
;   the final formula is represented as
;     f(X) = A[0]*(1-2*exp(-(X-TE/2)*A[1])+exp(-X*A[1]))
;   to do the calculation for curve fitting the function lmfit requires the partial derivatives
;   for all of the A[] matrix.
;     df/dA[0] = 1-2*exp(-(X-TE/2)*A[1])+exp(-X*A[1]))
;     df/dA[1] = A[0]*(2*(X-TE/2)*exp(-(X-TE/2)*A[1])-X*exp(-X*A[1]))


; Editing Information:
    ; adapted from SatRecovGEFunct

function SatRecovSEFunct, X , A
    COMMON T1_fit_data, TE
    RETURN,[ [ A[0]*(1-2*exp(-(X-TE/2)*A[1])+exp(-X*A[1])) ], $
             [ (1-2*exp(-(X-TE/2)*A[1])+exp(-X*A[1])) ], $
             [ A[0]*(2*(X-TE/2)*exp(-(X-TE/2)*A[1])-X*exp(-X*A[1])) ] $
           ]
END
; Subroutine name: satRecov
; Created by: Garrett Astary
; Calling Information: Called from roi_Curve_Fit or image_fit if image_type = 1

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; parameters:
;   dep_var     [in] this is the data that needs a curve fit to it. must be same dimensions as indp_var.
;   indp_var    [in] independent variable that is increasing.
;   A           [out] array that will contain the fit parameters
;   chisquare   [out] the chisquare for this fit
; description:
;

pro SatRecovSE , dep_var, indp_var , A , chisquare, sigma

    common CURVEFIT_POWELL, fit_params
    
    ;make an array for the values that we want to converge.
    A=dblarr(2)
    ;we need the last value of y - y at x where x is nearest to 0
    temp = min(indp_var, Min_Subscript)
    A[0] = max(dep_var)
    ;a[1] is a representation of the slope at the crossing of the y axis
    A[1]= double(dep_var(0)) / double(indp_var(0))

    ;the parameters we would like to fit. (1=fit, 0=do not fit) that parameter.
    fita = [1,1]

    coefs = LMFIT( indp_var , dep_var, A ,FUNCTION_NAME = 'SatRecovSEFunct' $
                      ,CONVERGENCE=converg, iter = iterations ,FITA = fita $
                      ,chisq = chisquare, covar = covarance, ITMIN=6,SIGMA =sigma )
 end

; Subroutine name: curvefit_adcFunction
; Created by: BT - 10/20/2008
; Calling Information:
;
; X    - B-value of measurement
; A[0] - S0 to be fit
; A[1] - D Diffusion Coefficient to be fit
;
; Bugs or Important Comments to Developers:
;
; Currently used with MPFIT. NOT LMFIT due to 
; LMFIT being non-good.
;
; Purpose of subroutine:
; 
; Representa ADC model
;
; Editing Information:
;
; MODIFICATION HISTORY:
; Changed to use MPFIT - 11/2008 - BT

function curvefit_adcAlphaFunction, X, A
   foo = exp( -(X*A[1])^A[2] )
   return, A[0] * foo
end

; Subroutine name: curvefit_adcFunction
; Created by: BT - 10/20/2008
; Calling Information:
;
; X    - B-value of measurement
; A[0] - S0 to be fit
; A[1] - D Diffusion Coefficient to be fit
;
; Bugs or Important Comments to Developers:
;
; Currently used with MPFIT. NOT LMFIT due to 
; LMFIT being non-good.
;
; Purpose of subroutine:
; 
; Representa ADC model
;
; Editing Information:
;
; MODIFICATION HISTORY:
; Changed to use MPFIT - 11/2008 - BT

function curvefit_adcFunction, X, A
   foo = exp(-X*A[1])
   return, A[0] * foo
   ;; if this function is to be used with IDL's built in LMFIT
   ;; function, then the next line should be used.
   ;;return, [ [ A[0] * foo ], [ foo ], [ -X*A[0]*foo ] ]
end

; Subroutine name: curvefit_adc_alpha
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
;
; MODIFICATION HISTORY:
; Changed to use MPFIT package - 11/2008 - BT

pro curvefit_adc_alpha, dep_var, indp_var, A, chisquare, sigma, noiseerr=noiseerr

    ;; make an array for the values that we want to converge.
    A = dblarr(3)
    
    ;; we need the last value of y - y at x where x is nearist to 0
    A[0] = dep_var[0]
    
    ;;a[1] is a representation of the slope at the crossing of the y axis
    A[1] = 1./20
    A[2] = 1.0
    
    if (keyword_set(noiseerr)) then begin
       err = noiseerr
       coefs = MPFITFUN('curvefit_adcAlphaFunction', indp_var, dep_var, err, A, $
                        covar=covariance, perror=sigma, dof=dof, $
                        yfit=yfit, bestnorm=bestnorm, /quiet)

    endif else begin
       coefs = MPFITFUN('curvefit_adcAlphaFunction', indp_var, dep_var, err, A, $
                        covar=covariance, perror=sigma, weights=1D, dof=dof, $
                        yfit=yfit, bestnorm=bestnorm, /quiet)
       
    endelse

    chisquare = sqrt(bestnorm/dof) ;;total((dep_var - yfit)^2)

    A = coefs
end

; Subroutine name: curvefit_adc
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
;
; MODIFICATION HISTORY:
; Changed to use MPFIT package - 11/2008 - BT

pro curvefit_adc, dep_var, indp_var, A, chisquare, sigma, noiseerr=noiseerr

    ;; make an array for the values that we want to converge.
    A = dblarr(2)
    
    ;; we need the last value of y - y at x where x is nearist to 0
    A[0] = dep_var[0]
    
    ;;a[1] is a representation of the slope at the crossing of the y axis
    A[1] = 1./20
    
    if (keyword_set(noiseerr)) then begin
       err = noiseerr
       coefs = MPFITFUN('curvefit_adcFunction', indp_var, dep_var, err, A, $
                        covar=covariance, perror=sigma, dof=dof, $
                        yfit=yfit, bestnorm=bestnorm, /quiet)

    endif else begin
       coefs = MPFITFUN('curvefit_adcFunction', indp_var, dep_var, err, A, $
                        covar=covariance, perror=sigma, weights=1D, dof=dof, $
                        yfit=yfit, bestnorm=bestnorm, /quiet)
       
    endelse

    chisquare = sqrt(bestnorm/dof) ;;total((dep_var - yfit)^2)
    A = coefs

end

; Subroutine name: multiSpinEchoFunc
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
; MODIFICATION HISTORY:
;   Edited by Alireza V 2006/10/16
;   Changing to the improved signal calculation.

function multiSpinEchoFunc, X , A
    RETURN,[ [A[0] * exp(-X*A[1]) ] , [ exp(-X*A[1]) ], [ -X*A[0]*exp(-X*A[1]) ] ]
end

; Subroutine name: multiSpinEcho
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.
; MODIFICATION HISTORY:
;   Edited by Alireza V 2006/10/16
;   Changing to the improved signal calculation.

pro multiSpinEcho , dep_var, indp_var , A , chisquare, sigma

    ;make an array for the values that we want to converge.
    A=dblarr(2)
    ;we need the last value of y - y at x where x is nearist to 0

    A[0] = dep_var[0]
    ;a[1] is a representation of the slope at the crossing of the y axis
    A[1]= 1/20
;    ;this is our guess as to where it crosses the y axis. we assume no noise.
;    A[2] = 0.00

    ;the parameters we would like to fit. (1=fit, 0=do not fit) that parameter.
    fita = [1,1]

;    coefs = MPFITFUN('ADCFunction', indp_var, dep_var, err, A, $
;                     covar=covariance, weights=1D, /quiet)
;    A = coefs

    coefs = LMFIT( indp_var , dep_var, A  ,FUNCTION_NAME = 'multiSpinEchoFunc' $
             ,CONVERGENCE=converg, iter = iterations ,FITA = fita $
             ,chisq = chisquare, covar = covarance, ITMIN=6,SIGMA =sigma)

end


; Subroutine name: make_multi_graphs
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro make_multi_graphs, ROI_Data, Fit_Data, fit_type, name_array
common T1_fit_data
    ;fist plot the T* vs avg under ROI

    ;prep the data
    n_roi = (size(ROI_DATA))[2]-1
    ;print, 'n_roi', n_roi

    time = rotate(ROI_DATA[*,0],2)
    ;print, 'time',time

    ROI_Avg  = rotate(ROI_DATA[*,1:*],5)
    ;print, 'ROI_Avg',ROI_Avg

	;create the colors to loop thourgh when over plotting
    ;color_names = ['Red','Black','Blue','Green','Brown','Cyan','Purple','Maroon','Gray','Orange']
		color_names = ['Red','Blue','Forest Green','Orange','Gray','Maroon','Brown','Purple','Dark Gray','Hot Pink']
		length_color_names = (size(color_names))[1]


    if fit_type eq 0 then begin

       TR = float(indgen(max(time)*10+1))/10.0


       Y=fltArr((size(TR))[1],n_roi)

       for ii=0 , n_roi-1 do $
         Y[*,ii]=Fit_Data[ii,0]*(1-exp(-TR*Fit_Data[ii,1]))+Fit_Data[ii,2]

       ;name the window
       windowName = 'T1 Adim Vs mean least squares fit ROI'

       data = Y[*,0]

       ;print, size(data)
       ;print, size(TR)

       iPlot_identifier = 6

       IPLOT, TR, data , $
            title = windowName, $
            IDENTIFIER=iPlot_identifier, $
            NAME='Fit '+name_array[0], $
            XTitle = 'TR (sec)', $
            YTitle = 'S', $
            COLOR = FSC_Color('Black', /Triple, /Row), $
            thick = 2.0


       IPLOT, time , ROI_Avg[*,0], $
            OVERPLOT=iPlot_identifier, $
            /SCATTER, $
            name='Data '+name_array[0], $
            SYM_size = 1.0, $0.3, $
            COLOR = FSC_Color('Black', /Triple, /Row), $
            SYM_INDEX = 6 ;1

; This next counter is for changing the symbol for each data point sequentially
; starts at 2 because there is one plot already printed before the FOR loop starts
		counter = 2
		line_style_counter = 2

; this is the size of the color vector counter to cycle the markers
		marker_counter = 0

        for ii=1, n_roi-1 do begin

         IPLOT, TR, Y[*,ii], $
                OVERPLOT=iPlot_identifier, $
                   name='Fit '+name_array[ii], $
                   COLOR = FSC_Color(color_names[(ii-1) mod 10], /Triple, /Row), $
                   thick = 2.0, $
                   LINESTYLE = line_style_counter



;          IPLOT, time , ROI_Avg[*,ii], $
;                OVERPLOT=iPlot_identifier, $
;                /SCATTER, $
;                name='Data '+name_array[ii], $
;                SYM_size = 0.3, $
;                SYM_INDEX = 1, $
;                COLOR = FSC_Color(color_names[ii-1 mod 10], /Triple, /Row)

; Increasing the symbol used on every data series
 			IPLOT, time , ROI_Avg[*,ii], $
                OVERPLOT=iPlot_identifier, $
                /SCATTER, $
                name='Data '+name_array[ii], $
                SYM_size = 1.0, $; 0.3, $
                SYM_INDEX = counter, $
                COLOR = FSC_Color(color_names[(ii-1+marker_counter) mod 10], /Triple, /Row)

			; Lets increase the counter until it reaches 8 which is not a valid point marker
			; Also, lets add an offset everytime it reaches the 10 colors so the marker and
			; the plot lines are not completely the same in the next iteration.
			if (counter eq 8) then counter = 2 else counter = counter + 1
			;cycle the linestyle
			if (line_style_counter eq 5) then line_style_counter = 1 else line_style_counter = line_style_counter + 1

			if (ii eq length_color_names) then marker_counter = marker_counter + 2

        end


    endif else if fit_type eq 1 then begin
    TR = float(indgen(max(time)*10+1))/10.0


       Y=fltArr((size(TR))[1],n_roi)

       for ii=0 , n_roi-1 do $
         Y[*,ii]=Fit_Data[ii,0]*(1-2*exp(-(TR-TE/2)*Fit_Data[ii,1])+exp(-TR*Fit_Data[ii,1]))+Fit_Data[ii,2]

       ;name the window
       windowName = 'T1 Adim Vs mean least squares fit ROI'

       data = Y[*,0]

       ;print, size(data)
       ;print, size(TR)

       iPlot_identifier = 6

       IPLOT, TR, data , $
            title = windowName, $
            IDENTIFIER=iPlot_identifier, $
            NAME='Fit '+name_array[0], $
            XTitle = 'TR (sec)', $
            YTitle = 'S', $
            COLOR = FSC_Color('Black', /Triple, /Row), $
            thick = 2.0


       IPLOT, time , ROI_Avg[*,0], $
            OVERPLOT=iPlot_identifier, $
            /SCATTER, $
            name='Data '+name_array[0], $
            SYM_size = 1.0, $0.3, $
            COLOR = FSC_Color('Black', /Triple, /Row), $
            SYM_INDEX = 6 ;1

; This next counter is for changing the symbol for each data point sequentially
; starts at 2 because there is one plot already printed before the FOR loop starts
    counter = 2
    line_style_counter = 2

; this is the size of the color vector counter to cycle the markers
    marker_counter = 0

        for ii=1, n_roi-1 do begin

         IPLOT, TR, Y[*,ii], $
                OVERPLOT=iPlot_identifier, $
                   name='Fit '+name_array[ii], $
                   COLOR = FSC_Color(color_names[(ii-1) mod 10], /Triple, /Row), $
                   thick = 2.0, $
                   LINESTYLE = line_style_counter



;          IPLOT, time , ROI_Avg[*,ii], $
;                OVERPLOT=iPlot_identifier, $
;                /SCATTER, $
;                name='Data '+name_array[ii], $
;                SYM_size = 0.3, $
;                SYM_INDEX = 1, $
;                COLOR = FSC_Color(color_names[ii-1 mod 10], /Triple, /Row)

; Increasing the symbol used on every data series
      IPLOT, time , ROI_Avg[*,ii], $
                OVERPLOT=iPlot_identifier, $
                /SCATTER, $
                name='Data '+name_array[ii], $
                SYM_size = 1.0, $; 0.3, $
                SYM_INDEX = counter, $
                COLOR = FSC_Color(color_names[(ii-1+marker_counter) mod 10], /Triple, /Row)

      ; Lets increase the counter until it reaches 8 which is not a valid point marker
      ; Also, lets add an offset everytime it reaches the 10 colors so the marker and
      ; the plot lines are not completely the same in the next iteration.
      if (counter eq 8) then counter = 2 else counter = counter + 1
      ;cycle the linestyle
      if (line_style_counter eq 5) then line_style_counter = 1 else line_style_counter = line_style_counter + 1

      if (ii eq length_color_names) then marker_counter = marker_counter + 2

       end
    endif else if fit_type eq 2 then begin

       TE = float(indgen(max(time)*1000+1))/1000.0
       ;TE = time
       ;print, 'TE',TE
       ;print,'size(TE)',size(TE)
;     print,'max(TE)',max(TE)
;
;
;     print,'S0 roi 1'
;     print,Fit_Data[0,0]
;     print,'T1 roi 1'
;     print,Fit_Data[0,1]


       Y=fltArr((size(TE))[1],n_roi)

       for ii=0 , n_roi-1 do $
         Y[*,ii]=Fit_Data[ii,0]*(exp( -TE*Fit_Data[ii,1] ))+Fit_Data[ii,2]


       ;name the window
       windowName = 'T2 Adim Vs mean least squares fit ROI'



       iPlot_identifier = 6

       IPLOT, TE, Y[*,0] , $
            title = windowName, $
            IDENTIFIER=iPlot_identifier, $
            NAME='Fit '+name_array[0], $
            XTitle = 'TE (sec)', $
            YTitle = 'S', $
            COLOR = FSC_Color('Black', /Triple, /Row), $
            thick = 2.0


       IPLOT, time , ROI_Avg[*,0], $
            OVERPLOT=iPlot_identifier, $
            /SCATTER, $
            name='Data '+name_array[0], $
            SYM_size = 1.0, $; 0.3, $
            SYM_INDEX = 6, $; 1, $
            COLOR = FSC_Color('Black', /Triple, /Row)

; This next counter is for changing the symbol for each data point sequentially
; starts at 2 because there is one plot already printed before the FOR loop starts
		counter = 2
		marker_counter = 0
		line_style_counter = 2

        for ii=1, n_roi-1 do begin

         IPLOT, TE, Y[*,ii], $
                OVERPLOT=iPlot_identifier, $
                   name='Fitted '+name_array[ii], $
                   COLOR = FSC_Color(color_names[(ii-1) mod 10], /Triple, /Row), $
                   thick = 2.0,$
                   LINESTYLE = line_style_counter



          IPLOT, time , ROI_Avg[*,ii], $
                OVERPLOT=iPlot_identifier, $
                /SCATTER, $
                name='Data '+name_array[ii], $
                SYM_size = 1.0, $; 0.3, $
                SYM_INDEX = counter, $
                COLOR = FSC_Color(color_names[(ii-1+marker_counter) mod 10], /Triple, /Row)

			if (counter eq 8) then counter = 2 else counter = counter + 1
			;cycle the linestyle
			if (line_style_counter eq 5) then line_style_counter = 1 else line_style_counter = line_style_counter + 1

			if (ii eq length_color_names) then marker_counter = marker_counter + 2


        end


    endif else if fit_type eq 3 then begin

       time = rotate((ROI_DATA[*,0]),2)

       B_Value = time;(float(indgen(max(time))))

       Y=fltArr((size(B_Value))[1],n_roi)

       for ii=0 , n_roi-1 do begin
           Y[*,ii]=Fit_Data[ii,0]*(exp( -B_Value*Fit_Data[ii,1] ))+Fit_Data[ii,2]
           ;;Y[*,ii]=Fit_Data[ii,0]*(exp( -(B_Value*Fit_Data[ii,1])^Fit_Data[ii,2] )) ;; +Fit_Data[ii,2]
       endfor

       ;name the window
       windowName = 'B Value Vs mean least squares fit ROI'


       iPlot_identifier =  6


       IPLOT, B_Value, Y[*,0] , $
            title = windowName, $
            IDENTIFIER=iPlot_identifier, $
            NAME='Fit ' +name_array[0], $
            XTitle = 'b (sec/mm^2)', $
            YTitle = 'S',$
            COLOR = FSC_Color('Black', /Triple, /Row)


       IPLOT, time , ROI_Avg[*,0], $
            OVERPLOT=iPlot_identifier, $
            /SCATTER, $
            name='Data '+name_array[0], $
            SYM_size = 1.0, $;0.3, $
            SYM_INDEX = 6, $; 1,$
            COLOR = FSC_Color('Black', /Triple, /Row)

; This next counter is for changing the symbol for each data point sequentially
; starts at 2 because there is one plot already printed before the FOR loop starts
		counter = 2
		marker_counter = 0
		line_style_counter = 2

        for ii=1, n_roi-1 do begin

         IPLOT, B_Value, Y[*,ii], $
                OVERPLOT=iPlot_identifier, $
                   name='Fit '+name_array[ii], $
                   COLOR = FSC_Color(color_names[(ii-1) mod 10], /Triple, /Row), $
                   thick = 2.0, $
                   LINESTYLE = line_style_counter



          IPLOT, time , ROI_Avg[*,ii], $
                OVERPLOT=iPlot_identifier, $
                /SCATTER, $
                name='Data '+name_array[ii], $
                SYM_size = 1.0, $ ;0.3, $
                SYM_INDEX = 6, $ ;1, $
                COLOR = FSC_Color(color_names[(ii-1+marker_counter) mod 10], /Triple, /Row)

			if (counter eq 8) then counter = 2 else counter = counter + 1
			;cycle the linestyle
			if (line_style_counter eq 5) then line_style_counter = 1 else line_style_counter = line_style_counter + 1

			if (ii eq length_color_names) then marker_counter = marker_counter + 2

        end

    endif else if fit_type eq 4 then begin

       time = rotate((ROI_DATA[*,0]),2)

       B_Value = time;(float(indgen(max(time))))

       Y=fltArr((size(B_Value))[1],n_roi)

       for ii=0 , n_roi-1 do begin
           Y[*,ii]=Fit_Data[ii,0]*(exp( -(B_Value*Fit_Data[ii,1])^Fit_Data[ii,2] )) ;; +Fit_Data[ii,2]
       endfor

       ;name the window
       windowName = 'B Value Vs mean least squares fit ROI'


       iPlot_identifier =  6


       IPLOT, B_Value, Y[*,0] , $
            title = windowName, $
            IDENTIFIER=iPlot_identifier, $
            NAME='Fit ' +name_array[0], $
            XTitle = 'b (sec/mm^2)', $
            YTitle = 'S',$
            COLOR = FSC_Color('Black', /Triple, /Row)


       IPLOT, time , ROI_Avg[*,0], $
            OVERPLOT=iPlot_identifier, $
            /SCATTER, $
            name='Data '+name_array[0], $
            SYM_size = 1.0, $;0.3, $
            SYM_INDEX = 6, $; 1,$
            COLOR = FSC_Color('Black', /Triple, /Row)

; This next counter is for changing the symbol for each data point sequentially
; starts at 2 because there is one plot already printed before the FOR loop starts
        counter = 2
        marker_counter = 0
        line_style_counter = 2

        for ii=1, n_roi-1 do begin

         IPLOT, B_Value, Y[*,ii], $
                OVERPLOT=iPlot_identifier, $
                   name='Fit '+name_array[ii], $
                   COLOR = FSC_Color(color_names[(ii-1) mod 10], /Triple, /Row), $
                   thick = 2.0, $
                   LINESTYLE = line_style_counter



          IPLOT, time , ROI_Avg[*,ii], $
                OVERPLOT=iPlot_identifier, $
                /SCATTER, $
                name='Data '+name_array[ii], $
                SYM_size = 1.0, $ ;0.3, $
                SYM_INDEX = 6, $ ;1, $
                COLOR = FSC_Color(color_names[(ii-1+marker_counter) mod 10], /Triple, /Row)

            if (counter eq 8) then counter = 2 else counter = counter + 1
            ;cycle the linestyle
            if (line_style_counter eq 5) then line_style_counter = 1 else line_style_counter = line_style_counter + 1

            if (ii eq length_color_names) then marker_counter = marker_counter + 2

        endfor

    endif else  if fit_type eq 5 then begin

       TI = float(indgen(max(time)*10+1))/10.0

       Y=fltArr((size(TI))[1],n_roi)

       for ii=0 , n_roi-1 do $
         Y[*,ii]=abs(Fit_Data[ii,0]*(1.0-2.0*exp(-TI*Fit_Data[ii,1])))
       ;name the window
       windowName = 'T1 Adim Vs mean least squares fit ROI'

       data = Y[*,0]

       ;print, size(data)
       ;print, size(TR)

       iPlot_identifier = 6

       IPLOT, TI, data , $
            title = windowName, $
            IDENTIFIER=iPlot_identifier, $
            NAME='Fit '+name_array[0], $
            XTitle = 'TI (sec)', $
            YTitle = 'S', $
            COLOR = FSC_Color('Black', /Triple, /Row), $
            thick = 2.0


       IPLOT, time , ROI_Avg[*,0], $
            OVERPLOT=iPlot_identifier, $
            /SCATTER, $
            name='Data '+name_array[0], $
            SYM_size = 1.0, $0.3, $
            COLOR = FSC_Color('Black', /Triple, /Row), $
            SYM_INDEX = 6 ;1

; This next counter is for changing the symbol for each data point sequentially
; starts at 2 because there is one plot already printed before the FOR loop starts
        counter = 2
        line_style_counter = 2

; this is the size of the color vector counter to cycle the markers
        marker_counter = 0

        for ii=1, n_roi-1 do begin

         IPLOT, TI, Y[*,ii], $
                OVERPLOT=iPlot_identifier, $
                   name='Fit '+name_array[ii], $
                   COLOR = FSC_Color(color_names[(ii-1) mod 10], /Triple, /Row), $
                   thick = 2.0, $
                   LINESTYLE = line_style_counter



;          IPLOT, time , ROI_Avg[*,ii], $
;                OVERPLOT=iPlot_identifier, $
;                /SCATTER, $
;                name='Data '+name_array[ii], $
;                SYM_size = 0.3, $
;                SYM_INDEX = 1, $
;                COLOR = FSC_Color(color_names[ii-1 mod 10], /Triple, /Row)

; Increasing the symbol used on every data series
            IPLOT, time , ROI_Avg[*,ii], $
                OVERPLOT=iPlot_identifier, $
                /SCATTER, $
                name='Data '+name_array[ii], $
                SYM_size = 1.0, $; 0.3, $
                SYM_INDEX = counter, $
                COLOR = FSC_Color(color_names[(ii-1+marker_counter) mod 10], /Triple, /Row)

            ; Lets increase the counter until it reaches 8 which is not a valid point marker
            ; Also, lets add an offset everytime it reaches the 10 colors so the marker and
            ; the plot lines are not completely the same in the next iteration.
            if (counter eq 8) then counter = 2 else counter = counter + 1
            ;cycle the linestyle
            if (line_style_counter eq 5) then line_style_counter = 1 else line_style_counter = line_style_counter + 1

            if (ii eq length_color_names) then marker_counter = marker_counter + 2

        end


    endif

;    legend_titles = strarr(n_roi*2)
;       for ii=0 , n_roi*2-1 do legend_titles[ii] = string('ROI ',STRTRIM(floor(ii/2)+1,2) )
;
;       ;print, legend_titles

;    ;Set the legend name...
;    ;Get the objects associated with the Live_Tool
;    v = obj_valid()
;
;    ;Get the legend object.
;    a = where(obj_isa(v, 'idlgrlegend'))
;    olegend=v[a[0]]
;
;    olegend->idlgrlegend::setproperty, item_name = legend_titles
;    live_control, data, /update

end


;*****************************************************************************************************
;
; NAME:
;   roiCurveFit
;
; PURPOSE:
;
; Bugs and/or Comments to Developers:
; T* is arranged so that it is increasing in time. Presently the subroutine assumes it decrease
; in time and this probably needs to be updated because data may come in already increasing in time.
; e.g. The data may have been collected so that T* is increasing
;******************************
;
; ARGUMENTS:
;   scan_data [in]    main structure that holds all the state variables
;   type   [in]    specifies what type of curve to fit.
;               0: progressive saturation
;               1: inversion saturation
;               2: multi spin echo
;
; DESCRIPTION:
;   this procedure is called from the main pro for curve fitting. this pro sets up and
;   does the curve fitting.
;   1. collect the global varables that are neccessary
;   2. ask the user to specify an roi
;   3. extract the mask from that roi
;   4. calculate the image statistics namely the mean across the adim
;   5. set up the curve fit by estimating the inital parameters. A[]
;   6. arrange the data so that it is increasing
;   7. calculate the curvefit parameters
;
; MODIFICATION HISTORY:
;   Edited by HS 2006/10/09
;   Fix spelling mistakes and commenting.
;
;*****************************************************************************************************
pro roi_Curve_Fit , ev
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON T1_fit_data, TE

    CI = project.ci
    sdimStart = project.procPramArray[CI].sdim_Start
    adimStart = project.procPramArray[CI].adim_Start
    multi_slice = project.procPramArray[ci].single_multi_flag
    graph_flag = project.procPramArray[ci].curvefit_graph
    
    adim = project.imndArray[CI].adim
    sdim = project.imndarray[CI].sdim
    intensityCen = project.procPramArray[project.ci].intensity_Cen
    intensityMax = project.procPramArray[project.ci].intensity_Max
    intensityMin = project.procPramArray[project.ci].intensity_Min
    type = project.procPramArray[project.ci].ROI_fit_type

    display_format = '(E0.4)'
    
    if project.roi.ci eq 0 then begin
        void= dialog_message(['Please select a different item in the ROI list',$
                        'The first item in the ROI list is for creating a new ROI'],/info)
        return
    end

    ;chk to make sure that the there is an roi to draw on.
    if not ptr_valid(project.roi.pROIs[project.roi.ci]) then begin
        void =  dialog_message('Please either make an ROI or change selection',/error)
        return
    end

    ;get the data object from the roi program.
;    widget_control, ev.id, get_uvalue=og
;    regions = og->get(/all)
    regions = (*project.roi.pROIs[project.roi.ci])

    ;make sure the image is in state1
    mas_load_state_1

    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.
    case project.procpramarray[ci].slice_axis of
       0: begin ;; RP
          sdimStart = project.procpramarray[ci].sdim_start
          sdim      = project.imndarray[ci].sdim
          p_image = ptr_new(reform((*project.dataArray[ci].state1)[*,*,sdimStart,adimStart]))
       end
       1: begin ;; RS
          sdimStart = project.procpramarray[ci].pdim_start
          sdim      = project.imndarray[ci].pdim
          p_image = ptr_new(reform((*project.dataArray[ci].state1)[*,sdimStart,*,adimStart]))
       end
       2: begin ;; PS
          sdimStart = project.procpramarray[ci].fdim_start
          sdim      = project.imndarray[ci].fdim
          p_image = ptr_new(reform((*project.dataArray[ci].state1)[sdimStart,*,*,adimStart]))
       end
    endcase

    if (multi_slice) then begin
       slice_start = 0
       slice_end   = sdim - 1
    endif else begin
       slice_start = sdimStart
       slice_end   = sdimStart
    endelse

;; FIXME ROI
    if (project.procpramarray[ci].no_transform_roi eq 1) then begin
        mas_zoom , p_image
        mas_rotate_flip , p_image
    endif
    
    p_image_sz = size(*p_image, /dimensions)
        
    ;clear the status bar
    update_status_bar, ''

    ;chk to make sure the object regions is valid
    if not(total(obj_valid(regions))) then begin
       ;update_status_bar, 'No region returned from ROI tool'
       return
    end

    ;how many roi's do we have
    if  (size(regions))[0] eq 0 then numRoi = 1 $
    else numRoi =  (size(regions))[1]

    ;make an array to hold all the names entered by the user
    name_array = strarr(numRoi)

    ; make a pointer array to hold all the masks
    mask = ptrarr(numRoi)

    ;extract the masks
    for temp=0 , numRoi-1 do begin
       mask[temp] = ptr_new(regions[temp]->ComputeMask(dimensions =[p_image_sz[0],p_image_sz[1]], $
                                                       MASK_RULE=2))
       regions[temp] -> GETPROPERTY, name= name
       name_array[temp] = name
    end

    roi_names = *project.roi.pdisplay_names
    have_noise_mask = 0
    for n = 0, n_elements(roi_names)-1 do begin
       if (strupcase(roi_names[n]) eq 'BACKGROUND') then begin
          obj_ptr = project.roi.prois[n]
          for r = 0, n_elements(*obj_ptr)-1 do begin
             if (obj_valid((*obj_ptr)[r])) then begin
                (*obj_ptr)[r]->getProperty, name=mask_name
                if (strupcase(mask_name) eq 'BACKGROUND') then begin
                   noise_mask = (*obj_ptr)[r] -> computeMask(dimensions=[ p_image_sz[0], p_image_sz[1]], mask_rule=2)
                   print, 'roi_curvefit: found BACKGROUND mask.'
                   have_noise_mask = 1
                endif
             endif
          endfor
       endif
    endfor

    ;;after we know how many roi's there are we can now make an
    ;;array to hold all the info so that we can
    ;;display it later on after the data is calculated.
    ;;the first entry in the numRoi column will be the times, either
    ;;TR or TE then the average intensity for each time
    ROI_Data = fltArr(adim,numRoi+1)

    ;now we create an array to hold the fit parameters and chisquare for each ROI
    Fit_Data = fltArr(numRoi,4)

    ;in order for the virtual machine version of idl to display text we need a text box to display
    ;our info in. sence the process "display_stats, display_string_array, title" creates and displays
    ;text info all i have to do is create an array to pass into it that has what i want to say.
    ;this is what i am doing here
    ;the first entry into the string it the collum headers for the info that will eventually
    ;be written to the string array also.

    if (multi_slice eq 1) then begin
       display_string_array = strarr((numRoi * sdim) + 1)
    endif else begin
       display_string_array = strarr(numRoi + 1) 
    endelse

    ;print some pretty data that the user will want to the terminal
    case type of 
 
       0: begin
          display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])+'S0'+string([09B]) $
                                   +'Sigma'+string([09B])+'T1'+string([09B])+'Sigma'+string([09B])$
                                   +'Chisquare'
       end

       1: begin
          display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])+'S0'+string([09B]) $
                                   +'Sigma'+string([09B])+'T1'+string([09B])+'Sigma'+string([09B])$
                                   +'Chisquare'
       end

       2: begin
          display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])+'S0'+string([09B]) $
                                   +'Sigma'+string([09B])+'T2'+string([09B])+'Sigma'+string([09B])$
                                   +'Chisquare'
       end

       3: begin
 ;;         display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])+'S0'+string([09B]) $
 ;;                                  +'Sigma'+string([09B])+'ADC mm^2/sec'+string([09B])+'Sigma'+string([09B])$
 ;;                                  +'Chisquare'
          display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])$
                                   +'S0'+string([09B])+'Sigma'+string([09B]) $
                                   +'ADC mm^2/sec'+string([09B])+'Sigma'+string([09B])$
                                   +'Chisquare'
       end

       4: begin
 ;;         display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])+'S0'+string([09B]) $
 ;;                                  +'Sigma'+string([09B])+'ADC mm^2/sec'+string([09B])+'Sigma'+string([09B])$
 ;;                                  +'Chisquare'
          display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])$
                                   +'S0'+string([09B])+'Sigma'+string([09B]) $
                                   +'ADC mm^2/sec'+string([09B])+'Sigma'+string([09B])$
                                   +'alpha'+string([09B])+'Sigma'+string([09B])$
                                   +'Chisquare'
       end
       
       5: begin
          display_string_array[0] ='ROI#'+string([09B])+'Slice'+string([09B])+'S0'+string([09B]) $
                                   +'Sigma'+string([09B])+'T1'+string([09B])+'Sigma'+string([09B])$
                                   +'Chisquare'
          end

       else: begin
          print, "Unknown FIT type"
          return
       end

    endcase

    maskCounter = 0
    line = 0

    for maskCounter= 0, numRoi-1 do begin

;;;
       for jj = slice_start, slice_end do begin 
       ;make an array to hold the average intensity for each adim
       avgIntensityInROI = fltarr(adim)
       stddevInROI = fltarr(adim)

       ;loop through the adim and where the roi mask is
       ; average that data and save it into the avgIntensityInROI matrix
       ii = 0
       for ii=0, adim-1 do begin
         ;we have to bring state1 back in b/c we byte scalled it earlyer
         ; and we need to extract data from the ffted image and not the bytescalled
         ; also the mask was drawn on the zoomed image so we have to zoom the image also
;;;         p_image =
;;;         ptr_new((*project.dataArray[ci].state1)[*,*,jj,ii])
          sdimStart = jj

          case project.procpramarray[ci].slice_axis of
             0: begin ;; RP
                p_image = ptr_new((*project.dataArray[ci].state1)[*,*,sdimStart,ii])
             end
             1: begin ;; RS
                p_image = ptr_new((*project.dataArray[ci].state1)[*,sdimStart,*,ii])
             end
             2: begin ;; PS
                p_image = ptr_new((*project.dataArray[ci].state1)[sdimStart,*,*,ii])
             end
          endcase
          *p_image = reform(*p_image)

          ;;p_image = ptr_new((*project.dataArray[ci].state1)[*,*,sdimStart,ii])

         ;; FIXME ROI
         if (project.procpramarray[ci].no_transform_roi eq 1) then begin
            mas_zoom , p_image
            mas_rotate_flip , p_image
         endif
         ;;mas_zoom , p_image
         ;;mas_rotate_flip , p_image

         ;image_statistics will calclualte the mean automatically
         ; then just save that to the avg intensity roi
         
         IMAGE_STATISTICS , *p_image, mask=*(mask[maskCounter]) $
                            , SUM_OF_SQUARES=image_sum_of_squares $
                            , VARIANCE=image_var $
                            , STDDEV=std_dev $
                            , Count=image_count
         avgIntensityInROI[ii] = float(sqrt(image_sum_of_squares / image_count - 2 * image_var))

         if (have_noise_mask eq 1) then begin
            IMAGE_STATISTICS, *p_image, mask=noise_mask, mean=noise_mean
            stddevInROI[ii] = noise_mean
            ;;print, jj, ii, abs(noise_mean - std_dev)
         endif else begin
            stddevInROI[ii] = std_dev
         endelse

      end

       ROI_Data[*,maskCounter+1] = avgIntensityInROI

       ;call the procedure that will do progresive saturation curve fitting
       if (type eq 0) then   begin
         ;need to flip the data so that it is increasing in time
         ;
         ;note***************************
         ; this probably needs to be updated because data may come in already increasing in time.
         ; the data may have been collected so that T* is increasing
         ;******************************
         indp_var = float(ROTATE(*project.imndArray[CI].rep_Time_ptr, 2))

         ;dep_var is assumed to be given on the scale of milli seconds and it is converted to seconds
         indp_var = indp_var / 1000

         ;after the indp_var is in the right units then save it to the overriding data structure
         ROI_Data[*,0] = rotate(indp_var,2)

         ;since we are flipping the indp variable then we need to flip the dependent var also.
         Collected_Data = float(ROTATE(avgIntensityInROI, 2))
         sigma = fltarr(5)

         SatRecovGE , Collected_Data, indp_var , A , chisquare, sigma


         ;store the a and chisquare data
         Fit_Data[maskCounter,0:1] = A
         Fit_Data[maskCounter,3] = chisquare

         ;for the graph to look right there needs to alot of intermediate values between the
         ;    points we are trying to fit.
         TR = float(indgen((max(indp_var)*10)+1))/10
         S=A[0]*(1-exp(-TR*A[1]))



       ;this is the display string formatted so that there are 10 digits and 4 digits to the
       ;right of the decimal point. after each number there is a tab except for the end.
;;        display_string_array[maskCounter+1]= $
        display_string_array[line+1]= $
             name_array[maskCounter]+ string([09B])+$
             strcompress(jj, /remove_all)+string([09B])+$
             string(FORMAT=display_format,A[0])    +string([09B])+$
             string(FORMAT=display_format,sigma[0])+string([09B])+$
             string(FORMAT=display_format,1/A[1]  )+string([09B])+$
             string(FORMAT=display_format,(sigma[1]/(A[1]^2)))+string([09B])+$
;;             string(FORMAT='(%"%10.2e")',A[2]    )+string([09B])+$
;;             string(FORMAT='(%"%10.2e")',sigma[2])+string([09B])+$
             string(FORMAT=display_format,chisquare)

       endif else if (type eq 1 ) then begin
       ;Saturation Recovery - allows for longitudinal relaxation between pi/2 and pi pulses
       ;need to flip the data so that it is increasing in time
         ;
         ;note***************************
         ; this probably needs to be updated because data may come in already increasing in time.
         ; the data may have been collected so that T* is increasing
         ;******************************
         indp_var = float(ROTATE(*project.imndArray[CI].rep_Time_ptr, 2))
         TE = project.imndArray[CI].echo_time / 1000
         ;dep_var is assumed to be given on the scale of milli seconds and it is converted to seconds
         indp_var = indp_var / 1000

         ;after the indp_var is in the right units then save it to the overriding data structure
         ROI_Data[*,0] = rotate(indp_var,2)

         ;since we are flipping the indp variable then we need to flip the dependent var also.
         Collected_Data = float(ROTATE(avgIntensityInROI, 2))
         sigma = fltarr(5)

         satRecovSE , Collected_Data, indp_var , A , chisquare, sigma


         ;store the a and chisquare data
         Fit_Data[maskCounter,0:1] = A
         Fit_Data[maskCounter,3] = chisquare

         ;for the graph to look right there needs to alot of intermediate values between the
         ;    points we are trying to fit.
         TR = float(indgen((max(indp_var)*10)+1))/10
         S=A[0]*(1-exp(-TR*A[1]))



       ;this is the display string formatted so that there are 10 digits and 4 digits to the
       ;right of the decimal point. after each number there is a tab except for the end.
;;        display_string_array[maskCounter+1]= $
        display_string_array[line+1]= $
             name_array[maskCounter]+ string([09B])+$
             strcompress(jj, /remove_all)+string([09B])+$
             string(FORMAT=display_format,A[0])    +string([09B])+$
             string(FORMAT=display_format,sigma[0])+string([09B])+$
             string(FORMAT=display_format,1/A[1]  )+string([09B])+$
             string(FORMAT=display_format,(sigma[1]/(A[1]^2)))+string([09B])+$
;;             string(FORMAT='(%"%10.2e")',A[2]    )+string([09B])+$
;;             string(FORMAT='(%"%10.2e")',sigma[2])+string([09B])+$
             string(FORMAT=display_format,chisquare)
         

       endif else if (type eq 2 ) then begin
         ;this is the multi spin echo
         indp_var = (float(*project.imndArray[CI].echo_Time_ptr))

         ;dep_var is assumed to be given on the scale of milli seconds and it is converted to seconds
         indp_var = indp_var / 1000

         ;after the indp_var is in the right units then save it to the overriding data structure
         ROI_Data[*,0] = indp_var

         Collected_Data = float(avgIntensityInROI)
    ;       print, 'indp_var',indp_var
    ;       print, 'Collected_Data', Collected_Data

         sigma = fltarr(3)
         multiSpinEcho, Collected_Data, indp_var , A , chisquare, sigma


         ;store the a and chisquare data
         Fit_Data[maskCounter,0:1] = A
         Fit_Data[maskCounter,3] = chisquare


         indp_var2 = float(indgen(max(indp_var)*1000))/1000



       ;this is the display string formatted so that there are 10 digits and 4 digits to the
       ;right of the decimal point. after each number there is a tab except for the end.
;;       display_string_array[maskCounter+1]= $
       display_string_array[line+1]= $
            name_array[maskCounter]+ string([09B])+$
            strcompress(jj, /remove_all)+string([09B])+$
            string(FORMAT=display_format,A[0])    +string([09B])+$
            string(FORMAT=display_format,sigma[0])+string([09B])+$
            string(FORMAT=display_format,1/A[1]  )+string([09B])+$
            string(FORMAT=display_format,sigma[1])+string([09B])+$
            $;string(FORMAT=display_format,(sigma[1]/(A[1]^2)))+string([09B])+$            
;            string(FORMAT=display_format,A[2]    )+string([09B])+$
;            string(FORMAT=display_format,sigma[2])+string([09B])+$
            string(FORMAT=display_format,chisquare)



          ;PRINT, FORMAT='(%"%d\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e\t%10.2e")', STRTRIM(maskCounter+1,2),A[0],sigma[0],1/A[1],sigma[1],A[2],sigma[2],chisquare

       endif else if (type eq 4) then begin

          ;;this is the ADC + alpha parameter

          ;;note***************************
          ;; this probably needs to be updated because data may come in already increasing in time.
          ;; the data may have been collected so that T* is increasing
          ;;******************************
          indp_var = float(*project.imndArray[CI].bval_Array)
          ;;dep_var is assumed to be given on the scale of milli seconds and it is converted to seconds
          ;;indp_var = indp_var / 1000.0

          ;;after the indp_var is in the right units then save it to the overriding data structure
          ROI_Data[*,0] = indp_var

          Collected_Data = float(avgIntensityInROI)
          ;;Collected_Data = rotate(Collected_Data,2)

          sigma = fltarr(2)
         
          ;;the ADC forumla is the same as the multiSpin echo formula
          curvefit_adc_alpha, Collected_Data, indp_var, A, chisquare, sigma, noiseerr=stddevInROI;;, /weight_with_ones

          ;;store the a and chisquare data
          Fit_Data[maskCounter,0:2] = A
          Fit_Data[maskCounter,3] = chisquare

          ;indp_var2 = float(indgen(max(indp_var)))

          ;;this is the display string formatted so that there are 10 digits and 4 digits to the
          ;;right of the decimal point. after each number there is a tab except for the end.
;;;          display_string_array[maskCounter+1]= $
          display_string_array[line+1]= $
             name_array[maskCounter]+ string([09B])+$
             strcompress(jj, /remove_all)+string([09B])+$
             string(FORMAT=display_format,A[0])    +string([09B])+$
             string(FORMAT=display_format,sigma[0])+string([09B])+$
             string(FORMAT=display_format,A[1]  )+string([09B])+$
             $;;string(FORMAT=display_format,(sigma[1]/(A[1]^2)))+string([09B])+$
             string(FORMAT=display_format,sigma[1])+string([09B])+$
             string(FORMAT=display_format,A[2])+string([09B])+$      ;; alpha parameter
             string(FORMAT=display_format,sigma[2])+string([09B])+$  ;; alpha sigma
             string(FORMAT=display_format,chisquare)

          
          endif else if (type eq 3) then begin

          ;;this is the ADC

          ;;note***************************
          ;; this probably needs to be updated because data may come in already increasing in time.
          ;; the data may have been collected so that T* is increasing
          ;;******************************
          indp_var = float(*project.imndArray[CI].bval_Array)
          ;;dep_var is assumed to be given on the scale of milli seconds and it is converted to seconds
          ;;indp_var = indp_var / 1000.0

          ;;after the indp_var is in the right units then save it to the overriding data structure
          ROI_Data[*,0] = indp_var

          Collected_Data = float(avgIntensityInROI)
          ;;Collected_Data = rotate(Collected_Data,2)

          sigma = fltarr(2)
         
          ;;the ADC forumla is the same as the multiSpin echo formula
          curvefit_adc, Collected_Data, indp_var, A, chisquare, sigma, noiseerr=stddevInROI;;, /weight_with_ones

          ;;store the a and chisquare data
          Fit_Data[maskCounter,0:1] = A
          Fit_Data[maskCounter,3] = chisquare

          ;indp_var2 = float(indgen(max(indp_var)))

          ;;this is the display string formatted so that there are 10 digits and 4 digits to the
          ;;right of the decimal point. after each number there is a tab except for the end.
;;;          display_string_array[maskCounter+1]= $
          display_string_array[line+1]= $
             name_array[maskCounter]+ string([09B])+$
             strcompress(jj, /remove_all)+string([09B])+$
             string(FORMAT=display_format,A[0])    +string([09B])+$
             string(FORMAT=display_format,sigma[0])+string([09B])+$
             string(FORMAT=display_format,A[1]  )+string([09B])+$
             $;;string(FORMAT=display_format,(sigma[1]/(A[1]^2)))+string([09B])+$
             string(FORMAT=display_format,sigma[1])+string([09B])+$
             string(FORMAT=display_format,chisquare)
             
       end else if (type eq 5) then begin ;; inversion recovery
          ;need to flip the data so that it is increasing in time
         ;
         ;note***************************
         ; this probably needs to be updated because data may come in already increasing in time.
         ; the data may have been collected so that T* is increasing
         ;******************************
         indp_var = float(ROTATE(*project.imndArray[CI].inversion_Time_ptr, 2))
         
         ;dep_var is assumed to be given on the scale of milli seconds and it is converted to seconds
         indp_var = indp_var / 1000.0

         ;after the indp_var is in the right units then save it to the overriding data structure
         ROI_Data[*,0] = rotate(indp_var,2)

         ;since we are flipping the indp variable then we need to flip the dependent var also.
         Collected_Data = float(ROTATE(avgIntensityInROI, 2))
         
         sigma = fltarr(5)
         
         InvRecovSE , Collected_Data, indp_var , A , chisquare, sigma, noiseerr=stddevInROI

         ;store the a and chisquare data
         Fit_Data[maskCounter,0:1] = A
         Fit_Data[maskCounter,3] = chisquare

         ;for the graph to look right there needs to alot of intermediate values between the
         ;    points we are trying to fit.
         TI = float(indgen((max(indp_var)*10)+1))/10
         S = abs(A[0]*(1.0-2.0*exp(-TI*A[1])))
       ;this is the display string formatted so that there are 10 digits and 4 digits to the
       ;right of the decimal point. after each number there is a tab except for the end.
;;        display_string_array[maskCounter+1]= $
        display_string_array[line+1]= $
             name_array[maskCounter]+ string([09B])+$
             strcompress(jj, /remove_all)+string([09B])+$
             string(FORMAT=display_format,A[0])    +string([09B])+$
             string(FORMAT=display_format,sigma[0])+string([09B])+$
             string(FORMAT=display_format,1/A[1]  )+string([09B])+$
             string(FORMAT=display_format,(sigma[1]/(A[1]^2)))+string([09B])+$
;;             string(FORMAT='(%"%10.2e")',A[2]    )+string([09B])+$
;;             string(FORMAT='(%"%10.2e")',sigma[2])+string([09B])+$
             string(FORMAT=display_format,chisquare)
       
       end


       line++

    endfor

;;;
 endfor

    display_stats, display_string_array, 'Fitting Parameters'

    print,''

    ;this the display string array for the data under each roi
    ;the first row is the title to the columns.

    data_string_array = strarr(adim+1)

    ;this will print the headers of each collum of the data that was fit.
    ;09B is the ascii form of tab
    if (type eq 0) then $
       printString = 'TR'
    if (type eq 1) then $
       printString = 'TR'
    if (type eq 2) then $
       printString ='TE'
    if (type eq 3 or type eq 4) then $
       printString = 'b val'
    if (type eq 5) then $
       printString = 'TR'
    for ii=0, numRoi-1 do $
        printString = printString + string([09B]) +name_array[ii]
    data_string_array[0] = printString

    ; now we print the collected data and format it correctly with tabs
    printROI_Data = rotate(ROI_Data,3)
    for ii=0, adim-1 do begin
       printString = strtrim(printROI_Data[0,ii],2)
       for jj=1, numroi do begin
         printString = printString + string([09B]) + $
                       string(FORMAT=display_format,printROI_Data[jj,ii])
         
       end
       data_string_array[ii+1] = printString
    endfor

    display_stats, data_string_array, 'Measured data for each ROI'

    if (multi_slice eq 0) and (graph_flag eq 1) then begin
       make_multi_graphs, ROI_Data, Fit_Data, type, name_array
    endif

end


; Subroutine name: multi_graphs_toggle_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro multi_graphs_toggle_event , event
    COMPILE_OPT IDL2
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci

    project.procPramArray[CI].multi_ROI_graphs = event.SELECT

end


;*****************************************************************************************************
;
; NAME:
;   roiCurveFit_init
;
; PURPOSE:
;
;
; ARGUMENTS:
;   scan_data [in]    main structure that holds all the state variables
;   type   [in]    specifies what type of curve to fit.
;               0: progresive saturation
;               1: invershions saturation
;               2: multi spin echo
;
; DESCRIPTION:
;   this procedure is called from the main pro for curve fitting. this pro sets up and
;   does the curve fitting.
;   1. collect the global variables that are neccessary
;   2. ask the user to specify an roi
;   3. extract the mask from that roi
;   4. calculate the image statistics namely the mean across the adim
;   5. set up the curve fit by estimating the inital parameters. A[]
;   6. arange the data so that it is increasing
;   7. calculate the curvefit parameters
;
; MODIFICATION HISTORY:
; Editing Information:
;   Edited by HS 2006/10/09
;   Fix spelling mistakes and commenting.
;*****************************************************************************************************

pro roi_Curve_Fit_init , type
    COMPILE_OPT IDL2


    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci
    sdimStart = project.procPramArray[CI].sdim_Start
    adimStart = project.procPramArray[CI].adim_Start
    adim = project.imndArray[CI].adim
    intensityCen = project.procPramArray[project.ci].intensity_Cen
    intensityMax = project.procPramArray[project.ci].intensity_Max
    intensityMin = project.procPramArray[project.ci].intensity_Min

    ;make sure the image is in state1
    mas_load_state_1


    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.
    p_image = ptr_new((*project.dataArray[ci].state1)[*,*,sdimStart,adimStart])
    mas_windowing , *p_image
    mas_zoom , p_image
    mas_rotate_flip , p_image

    ;tell the user to draw 1 ROI and then exit the ROI
    update_status_bar, 'Draw ROIs and then close the ROI program to begin calculations'


    tlb=widget_base(/column)

;    if project.roi_valid eq 1 then begin
;        ;mas_xroi, (*p_image), outgroup=og,apptlb=tlb, group=tlb, $
;            regions_in=*project.roi
;
;    end else begin
        ;mas_xroi, (*p_image), outgroup=og,apptlb=tlb, group=tlb
;    end


    b=widget_button(tlb,value='Calculate', event_pro='roi_Curve_Fit', uvalue=og)

    tlb2=widget_base(tlb,/row)
    GRAPH_BASE = Widget_Base(tlb2, UNAME='GRAPH_BASE  ' ,COLUMN=1 ,/NONEXCLUSIVE)

    b= Widget_Button(GRAPH_BASE, UNAME='multi_Graphs' ,/ALIGN_LEFT ,VALUE='Multiple Graphs' $
       , event_pro='multi_graphs_toggle_event'   )
    if project.procPramArray[CI].multi_ROI_graphs eq 1 then WIDGET_CONTROL, /SET_BUTTON, b

 end

; Subroutine name: batch_T2_image_fit
; Created by: Garrett Astary
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine: Performs T2 fit for a multi-slice variable TE (multi spin echo) data set. The pixel-by-pixel fit is stored in a
; fdim x pdim array and then exported to a text file as a column vector with a 2-line header giving the fdim and pdim values.
; Editing Information:


pro batch_T2_image_fit
	COMPILE_OPT IDL2

  COMMON scan_data
  CI = project.ci
    sdimStart = project.procPramArray[CI].sdim_Start
    S0_min = project.procPramArray[project.ci].S0_min
    S0_max = project.procPramArray[project.ci].S0_max
    T_min = project.procPramArray[project.ci].T_min
    T_max = project.procPramArray[project.ci].T_max

    intensity_min = project.procPramArray[project.ci].intensity_min
    intensity_cen = project.procPramArray[project.ci].intensity_cen
    intensity_max = project.procPramArray[project.ci].intensity_max

    fdim = project.ImndArray[CI].fdim
    pdim = project.ImndArray[CI].pdim
    ffov = project.imndArray[ci].f_fov
    pfov = project.imndArray[ci].p_fov
    sfov = project.imndArray[ci].s_fov
    rotate_Direction = project.procpramArray[ci].rotate_Direction
    slice_axis = project.procpramArray[ci].slice_axis

    image_fit_threshold_const = project.procPramArray[project.ci].image_fit_threshold

    image_fit_flag = project.procPramArray[project.ci].image_fit_flag
    type = project.procPramArray[project.ci].image_fit_type
    ;need to make sure the image is ffted rotated and zoomed if not do so
    mas_load_state_1




    T2_stack = fltarr(fdim, pdim, project.imndArray[ci].slices)
    ;for loop to do the following slice by slice
    for slice=0, project.imndArray[ci].slices -1 DO BEGIN

    p_image =  ptr_new(reform((*project.dataArray[ci].state1)[*,*,slice,*]))

    sz_image = size(*p_image)

    ;now we need to generate a mask for our data set.
    ;since we have multiple slices we have to generate it from the begining
    ;and loop through then freq then phase then slice. that way we dont
    ;have to recreate them unnecessarly.
    mask_background_flag =  project.procPramArray[CI].remove_background_active
    IF mask_background_flag EQ 1 THEN BEGIN
       key_adim = project.procPramArray[CI].key_adim
       FOR ii=0, sz_image[3]-1 DO BEGIN
         p_image_mask = PTR_NEW(REFORM((*p_image)[*,*,key_adim]))
         (*p_image)[*,*,ii] *= FLOAT(*(generate_background_mask( p_image_mask)))
       END

    END

    mas_zoom , p_image
    mas_rotate_flip , p_image

    ;use the roi mask
    FOR ii=0, sz_image[3]-1 DO BEGIN
       p_image_temp = ptr_new((*p_image)[*,*,ii])
       mas_roi_mask, p_image_temp
       (*p_image)[*,*,ii] = *p_image_temp
    END

    sz_image = size(*p_image)

    xx = 0
    yy = 0


;   print, 'image_fit'

    ;set the threshold
    image_fit_threshold = max((*p_image)[*,*,0]) * image_fit_threshold_const

    S0 = fltarr(sz_image[1],sz_image[2])

     T2 = fltarr((size((*p_image)))[1],(size((*p_image)))[2])



         ;prep the indp var for each different type
         print,*project.imndArray[CI].echo_Time_ptr
             indp_var = float(*project.imndArray[CI].echo_Time_ptr)
         indp_var = indp_var / 1000

         maximum = max( (*p_image)[*,*,0])

;          print, 'maximum',maximum
;          print, 'threshold',image_fit_threshold


         ;update_status_bar, 'processing image'
         progressbar = Obj_New('progressbar', Color='red', Text='T2 Image fit Slice '+strtrim(string(slice+1),2)+'/'+strtrim(string(project.imndArray[ci].slices),2) ,/NOCANCEL)
         progressbar -> Start

         for xx=0 , ((size((*p_image)))[1]-1) do begin
          for yy=0 , ((size((*p_image)))[2]-1) do begin

              if (*p_image)[xx,yy,0] gt image_fit_threshold then begin

                 dep_var =  reform((*p_image)[xx,yy,*])
                 multiSpinEcho, dep_var, indp_var , A , chisquare
                 S0[xx,yy] = A[0]
                 T2[xx,yy] = 1/A[1]
               endif
             endfor
             ;update_status_bar, 'processing image '+ STRTRIM(xx,2)+'/'+STRTRIM(((size((*p_image)))[1]),2)
             progressBar -> Update, (float(xx)/float(((size((*p_image)))[2]-1)))*100.0
         endfor
         progressbar -> Destroy
         T2_neg = where(T2 lt 0, count_T2)
         IF count_T2 gt 0 THEN BEGIN
          T2[T2_neg] = 0
         ENDIF

         T2_stack[*,*,slice] = T2
         endfor

    ;scaling from 0-255 for Amira
    max_T2 = max(max(max(T2_stack)))
    T2_stack = T2_stack/max_T2
    T2_stack = T2_stack*255



    file = (Dialog_PickFile( Title = 'Save Current Scan to flt', path=project.current_path))

      if file eq '' then return
      ;print, file
      image_type = project.imndArray[ci].image_type

      ;get the file extension that goes at the end of the file name
      ext = get_file_extension( image_type )

      ;update_status_bar,'Saving .flt file...'
      sz_data = size(T2_stack)
      OpenW, lun0, file+ext, /SWAP_IF_LITTLE_ENDIAN, /get_lun
      WriteU, lun0, sz_data, T2_stack
      Free_LUN, lun0

      write_geo_info_file, file, sz_data

;         ;The array resulting from the curve-fit, T2, is written as a column vector to a .txt file
;         snum = strtrim(string(slice+1),2)
;         OpenW, lun,batch_fpath+'T2_image_slice' + snum + '.txt', /get_lun
;         ;printing 2-line header
;         PrintF, lun, 'fdim = ' + strtrim(string(project.imndArray[ci].fdim),2)
;         PrintF, lun, 'pdim = ' + strtrim(string(project.imndArray[ci].pdim),2)
;         ;printing column vector of data
;         sz_T2 = size(T2)
;                  for ii=0,sz_T2[4]-1 do begin
;                     PrintF, lun, strtrim(string(T2[ii]),2)
;                  endfor
;
;
;         close, lun
;         free_lun, lun



  ;endfor   before I was writing the text file out as the T2 fit was completed, now I need to write an FLT file so we need a stack
  ; of the T2 curve-fit data
end

;+
; :Description:
;    Prepares a curvefit image for display and calls the right displayer
;
; :Params:
;    imgfit - raw fit data
;    title  - title of the window which will display the image
;
; :Keywords:
;    sclmin - optional, set the min value to be used when scaling for display,
;             for example, set this to 0 to have all negative fit values set to
;             zero before display
;    sclmax - optional, set the max value for display. For example, set this to
;             50 to have all values greater than 50 set to 255 on display.
;
; :Author: btt
;-
pro image_fit_display_image, imgfit, title, $
                             sclmin=sclmin, sclmax=sclmax

    common scan_data

    ci = project.ci
    
    if (n_elements(title) eq 0) then title = ''
    if (n_elements(sclmin) eq 0) then begin
        sclmin = min(imgfit)
    endif
    
    if (n_elements(sclmax) eq 0) then begin
        sclmax = max(imgfit)
    endif

    ffov = project.imndArray[ci].f_fov
    pfov = project.imndArray[ci].p_fov
    sfov = project.imndArray[ci].s_fov
    rotate_Direction = project.procpramArray[ci].rotate_Direction
    slice_axis = project.procpramArray[ci].slice_axis

    if ((rotate_direction+1) mod 2 eq 0) then begin
        xfov = pfov & yfov = ffov
    endif else begin
        xfov = ffov & yfov = pfov
    endelse
    
    if (project.iImage_flag eq 0) then begin
    
        p_imgfit = ptr_new(imgfit)

        if (project.procpramarray[ci].no_transform_roi eq 0) then begin
            mas_rotate_flip, p_imgfit
            mas_zoom, p_imgfit
        endif

        imgfit_scl = bytscl(*p_imgfit, min=sclmin, max=sclmax)
                
        imgdata = *p_imgfit & ptr_free, p_imgfit
        
        mas_display_multi, imgfit_scl, $
                           tab_title=title,  $
                           /standalone,$
                           fov_x=xfov, $
                           fov_y=yfov, $
                           fp_vals=imgdata, $
                           show_axis_labels=0

    endif else begin
         iImage,  imgfit $
           ,IMAGE_DIMENSIONS=[xfov,yfov] $
           ,xtitle='cm' $
           ,ytitle='cm' $
           ,title = title

    endelse

end

; Subroutine name: image_fit
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro image_fit, use_ci=use_ci
    COMPILE_OPT IDL2

    COMMON scan_data
    common T1_fit_data, TE
    
    if (n_elements(use_ci) ne 0) then begin
        if (not ptr_valid(project.dataarray[use_ci].state1)) then begin
            void = dialog_message(['No data found.', $
                                   'Please use the main Display button to use load data.'], $
                                  /center, /error)
            return
        endif
        CI = use_ci
    endif else begin
        CI = project.ci
    endelse
    
    sdimStart = project.procPramArray[CI].sdim_Start
    S0_min = project.procPramArray[CI].S0_min
    S0_max = project.procPramArray[CI].S0_max
    T_min = project.procPramArray[CI].T_min
    T_max = project.procPramArray[CI].T_max

    intensity_min = project.procPramArray[CI].intensity_min
    intensity_cen = project.procPramArray[CI].intensity_cen
    intensity_max = project.procPramArray[CI].intensity_max

    ffov = project.imndArray[ci].f_fov
    pfov = project.imndArray[ci].p_fov
    sfov = project.imndArray[ci].s_fov
    rotate_Direction = project.procpramArray[ci].rotate_Direction
    slice_axis = project.procpramArray[ci].slice_axis

    image_fit_threshold = project.procPramArray[CI].image_fit_threshold

    image_fit_flag = project.procPramArray[CI].image_fit_flag
    type = project.procPramArray[CI].image_fit_type
    ;need to make sure the image is ffted rotated and zoomed if not do so
    mas_load_state_1

    p_image =  ptr_new(reform((*project.dataArray[ci].state1)[*,*,sdimStart,*]))

    ;; FIXME ROI
    if (project.procpramarray[ci].no_transform_roi eq 1) then begin
        mas_zoom , p_image
        mas_rotate_flip , p_image
    endif

    sz_image = size(*p_image)

    ;now we need to generate a mask for our data set.
    ;since we have multiple slices we have to generate it from the begining
    ;and loop through then freq then phase then slice. that way we dont
    ;have to recreate them unnecessarly.
    mask_background_flag =  project.procPramArray[CI].remove_background_active
    IF mask_background_flag EQ 1 THEN BEGIN
       key_adim = project.procPramArray[CI].key_adim
       FOR ii=0, sz_image[3]-1 DO BEGIN
         p_image_mask = PTR_NEW(REFORM((*p_image)[*,*,key_adim]))
         (*p_image)[*,*,ii] *= FLOAT(*(generate_background_mask( p_image_mask)))
       END

    END

    ;use the roi mask
    FOR ii=0, sz_image[3]-1 DO BEGIN
       p_image_temp = ptr_new((*p_image)[*,*,ii])
       mas_roi_mask, p_image_temp
       (*p_image)[*,*,ii] = *p_image_temp
    END

    sz_image = size(*p_image)

    xx = 0
    yy = 0

    ;set the threshold
    image_fit_threshold = max((*p_image)[*,*,0]) * image_fit_threshold

    S0 = fltarr(sz_image[1],sz_image[2])

    if type eq 0 then begin
       ;progressive saturation

       if image_fit_flag eq 0 then begin

         T1 = fltarr(sz_image[1],sz_image[2])

         ;prep the indp var for each different type
         indp_var = float(ROTATE(*project.imndArray[CI].rep_Time_ptr, 2))


         indp_var = indp_var / 1000
         progressbar = Obj_New('progressbar', Color='red', Text='T1 Image fit',/NOCANCEL)
         progressbar -> Start


         for xx=0 , sz_image[1]-1 do begin
          for yy=0 , sz_image[2]-1 do begin

             if (*p_image)[xx,yy,0] gt image_fit_threshold then BEGIN

                     dep_var = ROTATE( reform((*p_image)[xx,yy,*]) , 2)
                     SatRecovGE , dep_var, indp_var , A , chisquare
                     S0[xx,yy] = A[0]
                     T1[xx,yy] = 1/A[1]
             end

             endfor

             progressBar -> Update, (float(xx)/float((sz_image[2]-1)))*100.0
           endfor
         progressbar -> Destroy

         ;store the images of S0 and T1
         p_image = ptr_new(fltarr(2,(size((*p_image)))[1],(size((*p_image)))[2]))
         (*p_image)[0,*,*] = S0
         (*p_image)[1,*,*] = T1
         project.dataArray[ci].imgFit1 = p_image

         project.procPramArray[CI].image_fit_flag = 1

       endif else begin

         S0 = reform((*project.dataArray[ci].imgFit1)[0,*,*])
         T1 = reform((*project.dataArray[ci].imgFit1)[1,*,*])

       end

       IMAGE_STATISTICS , S0, MEAN=S0mean, MINIMUM=S0imageMin , MAXIMUM=S0imageMax , STDDEV=S0stddev
       IMAGE_STATISTICS , T1, MEAN=T1mean, MINIMUM=T1imageMin , MAXIMUM=T1imageMax, STDDEV=T1stddev


       sData = [ 'S0image min'+string(S0imageMin),$
              'S0image max'+string(S0imageMax),$
                     'S0mean     '+string(S0mean),$
                 'S0stddev   '+string(S0stddev),$
              'T1image min'+string(T1imageMin),$
                 'T1image max'+string(T1imageMax),$
                 'T1mean     '+string(T1mean),$
                 'T1stddev   '+string(T1stddev)]

       display_stats, sData, 'Fitting parameters'
       
       image_fit_display_image, S0, 'T1-S0 (arbitary floating point values)', $
                                sclmin=S0_min, sclmax=S0_max
       image_fit_display_image, T1, 'T1 (floating point values in secs)', $
                                sclmin=T_min, sclmax=T_max

    end else if type eq 1 then begin
     
      if image_fit_flag eq 0 then begin
           
           T1 = fltarr(sz_image[1],sz_image[2])
  
           ;prep the indp var for each different type
           indp_var = float(ROTATE(*project.imndArray[CI].rep_Time_ptr, 2))
           TE = project.imndArray[CI].echo_time / 1000
  
           indp_var = indp_var / 1000
           ;update_status_bar, 'processing image'
           progressbar = Obj_New('progressbar', Color='red', Text='T1 Image fit',/NOCANCEL)
           progressbar -> Start
  
  
           for xx=0 , sz_image[1]-1 do begin
            for yy=0 , sz_image[2]-1 do begin
  
               if (*p_image)[xx,yy,0] gt image_fit_threshold then BEGIN
  
                       dep_var = ROTATE( reform((*p_image)[xx,yy,*]) , 2)
                       SatRecovSE , dep_var, indp_var , A , chisquare
                       S0[xx,yy] = A[0]
                       T1[xx,yy] = 1/A[1]
               end
  
               endfor
  
               progressBar -> Update, (float(xx)/float((sz_image[2]-1)))*100.0
             endfor
           progressbar -> Destroy
  
  
  
           ;store the images of S0 and T1
           p_image = ptr_new(fltarr(2,(size((*p_image)))[1],(size((*p_image)))[2]))
           (*p_image)[0,*,*] = S0
           (*p_image)[1,*,*] = T1
           project.dataArray[ci].imgFit1 = p_image
  
           project.procPramArray[CI].image_fit_flag = 1
  
         endif else begin
  
           S0 = reform((*project.dataArray[ci].imgFit1)[0,*,*])
           T1 = reform((*project.dataArray[ci].imgFit1)[1,*,*])
  
  
         end
  
         IMAGE_STATISTICS , S0, MEAN=S0mean, MINIMUM=S0imageMin, MAXIMUM=S0imageMax, STDDEV=S0stddev
         IMAGE_STATISTICS , T1, MEAN=T1mean, MINIMUM=T1imageMin, MAXIMUM=T1imageMax, STDDEV=T1stddev
  
  
         sData = [ 'S0image min'+string(S0imageMin),$
                'S0image max'+string(S0imageMax),$
                       'S0mean     '+string(S0mean),$
                   'S0stddev   '+string(S0stddev),$
                'T1image min'+string(T1imageMin),$
                   'T1image max'+string(T1imageMax),$
                   'T1mean     '+string(T1mean),$
                   'T1stddev   '+string(T1stddev)]
  
         display_stats, sData, 'Fitting parameters'
         image_fit_display_image, S0, 'T1-S0 (arbitary floating point values)', $
                                  sclmin=S0_min, sclmax=S0_max
         image_fit_display_image, T1, 'T1 (floating point values in secs)', $
                                  sclmin=T_min, sclmax=T_max
         
    end else if type eq 2 then begin
       ;multi spin echo
       T2 = fltarr((size((*p_image)))[1],(size((*p_image)))[2])


       if image_fit_flag eq 0 then begin
         ;prep the indp var for each different type
         print,*project.imndArray[CI].echo_Time_ptr
         indp_var = float(*project.imndArray[CI].echo_Time_ptr)
         indp_var = indp_var / 1000

         maximum = max( (*p_image)[*,*,0])

;          print, 'maximum',maximum
;          print, 'threshold',image_fit_threshold


         ;update_status_bar, 'processing image'
         progressbar = Obj_New('progressbar', Color='red', Text='T2 Image fit',/NOCANCEL)
         progressbar -> Start

         for xx=0 , ((size((*p_image)))[1]-1) do begin
          for yy=0 , ((size((*p_image)))[2]-1) do begin

              if (*p_image)[xx,yy,0] gt image_fit_threshold then begin

                 dep_var =  reform((*p_image)[xx,yy,*])
                 multiSpinEcho, dep_var, indp_var , A , chisquare
                 S0[xx,yy] = A[0]
                 T2[xx,yy] = 1/A[1]
               endif
             endfor
             ;update_status_bar, 'processing image '+ STRTRIM(xx,2)+'/'+STRTRIM(((size((*p_image)))[1]),2)
             progressBar -> Update, (float(xx)/float(((size((*p_image)))[2]-1)))*100.0
         endfor
         progressbar -> Destroy

         ;store the images of S0 and T1
         p_image = ptr_new(fltarr(2,(size((*p_image)))[1],(size((*p_image)))[2]))
         (*p_image)[0,*,*] = S0
         (*p_image)[1,*,*] = T2
         project.dataArray[ci].imgFit1 = p_image

         project.procPramArray[CI].image_fit_flag = 1

       endif else begin

         S0 = reform((*project.dataArray[ci].imgFit1)[0,*,*])
         T2 = reform((*project.dataArray[ci].imgFit1)[1,*,*])

       end

       IMAGE_STATISTICS , S0, MEAN=S0mean, MINIMUM=S0imageMin , MAXIMUM=S0imageMax , STDDEV=S0stddev
       IMAGE_STATISTICS , T2, MEAN=T2mean, MINIMUM=T2imageMin , MAXIMUM=T2imageMax, STDDEV=T2stddev
       
       sData = [ 'S0image min'+string(S0imageMin),$
                 'S0image max'+string(S0imageMax),$
                 'S0mean     '+string(S0mean),$
                 'S0stddev   '+string(S0stddev),$
                 'T2image min'+string(T2imageMin),$
                 'T2image max'+string(T2imageMax),$
                 'T2mean     '+string(T2mean),$
                 'T2stddev   '+string(T2stddev)]

       display_stats, sData, 'Fitting parameters'
       image_fit_display_image, S0, 'S0 for T2 (arbitary floating point values)', $
                                sclmin=S0_min, sclmax=S0_max
       image_fit_display_image, T2, 'T2 (floating point values in secs)', $
                                sclmin=T_min, sclmax=T_max

    end else if type eq 3 then begin
       
       ;;ADC
       img_size = size(*p_image)

       B_Val = fltarr(img_size[1],img_size[2])
       
       if image_fit_flag eq 0 then begin
          ;;prep the indp var for each different type
          indp_var = float(*project.imndArray[CI].bval_Array)
          ;;indp_var = indp_var / 1000.0
;          print, indp_var

          progressbar = Obj_New('progressbar', Color='red', Text='ADC Image fit', /NOCANCEL)
          progressbar -> Start
          
          for xx=0 , (img_size[1]-1) do begin
             for yy=0 , (img_size[2]-1) do begin

                if (*p_image)[xx,yy,0] gt image_fit_threshold then begin
                   dep_var =  reform((*p_image)[xx,yy,*])
                   curvefit_adc, dep_var, indp_var , A , chisquare
                   S0[xx,yy] = A[0]
                   B_Val[xx,yy] = A[1]
                endif

             endfor
             
            progressBar -> Update, (float(xx)/float((img_size[2]-1)))*100.0

         endfor

         progressbar -> Destroy

         ;;store the images of S0 and T1
         p_image = ptr_new(fltarr(2,img_size[1],img_size[2]))
         (*p_image)[0,*,*] = S0
         (*p_image)[1,*,*] = B_Val

         project.dataArray[ci].imgFit1 = p_image
         ;;project.procPramArray[CI].image_fit_flag = 1
         ;; commented out so that the fit is redone every time.
         ;; there is a conflict between the "Standard" fit and the "alpha-fit"
       endif else begin

         S0 = reform((*project.dataArray[ci].imgFit1)[0,*,*])
         B_Val = reform((*project.dataArray[ci].imgFit1)[1,*,*])

       end

       IMAGE_STATISTICS , S0, MEAN=S0mean, MINIMUM=S0_imageMin , MAXIMUM=S0_imageMax , STDDEV=S0stddev
       IMAGE_STATISTICS , B_Val, MEAN=B_Val_mean, MINIMUM=B_Val_imageMin , MAXIMUM=B_Val_imageMax, STDDEV=B_Val_stddev

       sData = [ 'S0image min'+string(S0_imageMin),$
                 'S0image max'+string(S0_imageMax),$
                 'S0mean     '+string(S0mean),$
                 'S0stddev   '+string(S0stddev),$
                 'ADCimage min'+string(B_Val_imageMin),$
                 'ADCimage max'+string(B_Val_imageMax),$
                 'ADC_mean     '+string(B_Val_mean),$
                 'ADC_stddev   '+string(B_Val_stddev)]

       display_stats, sData, 'Fitting parameters'
       
       image_fit_display_image, S0, 'S0 for dif fit (arbitary floating point values)', $
                                sclmin=S0_min, sclmax=S0_max
       image_fit_display_image, B_Val, 'ADC for dif fit (floating point values in um^2/msecs)', $
                                sclmin=T_min, sclmax=T_max
               
    end else if type eq 4 then begin
       
       ;;ADC
       img_size = size(*p_image)

       B_Val = fltarr(img_size[1],img_size[2])
       alpha = fltarr(img_size[1],img_size[2])
       
       if image_fit_flag eq 0 then begin
          ;;prep the indp var for each different type
          indp_var = float(*project.imndArray[CI].bval_Array)
          ;;indp_var = indp_var / 1000.0
;          print, indp_var

          progressbar = Obj_New('progressbar', Color='red', Text='ADC Image fit', /NOCANCEL)
          progressbar -> Start
          
          for xx=0 , (img_size[1]-1) do begin
             for yy=0 , (img_size[2]-1) do begin

                if (*p_image)[xx,yy,0] gt image_fit_threshold then begin
                   dep_var =  reform((*p_image)[xx,yy,*])
                   curvefit_adc_alpha, dep_var, indp_var , A , chisquare
                   S0[xx,yy]    = A[0]
                   B_Val[xx,yy] = A[1]
                   alpha[xx,yy] = A[2]
                endif

             endfor
             
            progressBar -> Update, (float(xx)/float((img_size[2]-1)))*100.0

         endfor

         progressbar -> Destroy

         ;;store the images of S0 and T1
         p_image = ptr_new(fltarr(3,img_size[1],img_size[2]))
         (*p_image)[0,*,*] = S0
         (*p_image)[1,*,*] = B_Val
         (*p_image)[2,*,*] = alpha
         
         project.dataArray[ci].imgFit1 = p_image

         ;;project.procPramArray[CI].image_fit_flag = 1
         ;; commented out... see above (same region in previous block)
       endif else begin

         S0 = reform((*project.dataArray[ci].imgFit1)[0,*,*])
         B_Val = reform((*project.dataArray[ci].imgFit1)[1,*,*])
         alpha = reform((*project.dataArray[ci].imgFit1)[2,*,*])
         
       end

       IMAGE_STATISTICS , S0, MEAN=S0mean, MINIMUM=S0_imageMin , MAXIMUM=S0_imageMax , STDDEV=S0stddev
       IMAGE_STATISTICS , B_Val, MEAN=B_Val_mean, MINIMUM=B_Val_imageMin , MAXIMUM=B_Val_imageMax, STDDEV=B_Val_stddev
       IMAGE_STATISTICS , alpha, MEAN=alpha_mean, MINIMUM=alpha_imageMin , MAXIMUM=alpha_imageMax, STDDEV=alpha_stddev

       sData = [ 'S0image min'+string(S0_imageMin),$
                 'S0image max'+string(S0_imageMax),$
                 'S0mean     '+string(S0mean),$
                 'S0stddev   '+string(S0stddev),$
                 'ADCimage min'+string(B_Val_imageMin),$
                 'ADCimage max'+string(B_Val_imageMax),$
                 'ADC_mean     '+string(B_Val_mean),$
                 'ADC_stddev   '+string(B_Val_stddev), $
                 'alpha min   '+string(alpha_imageMin), $
                 'alpha max   '+string(alpha_imageMax), $
                 'alpha_mean   '+string(alpha_mean), $
                 'alpha_stddev   '+string(alpha_stddev) $
                 ]

       display_stats, sData, 'Fitting parameters'
       image_fit_display_image, S0, 'S0 for dif fit (arbitary floating point values)', $
                                sclmin=S0_min, sclmax=S0_max
       image_fit_display_image, B_Val, 'ADC for dif fit (floating point values in um^2/msecs)', $
                                sclmin=T_min, sclmax=T_max
       image_fit_display_image, alpha, 'alpha for dif fit (arbitrary floating point values)'

    end else if type eq 5 then begin
       ;inversion recovery

       if 1 or image_fit_flag eq 0 then begin

         T1 = fltarr(sz_image[1],sz_image[2])

         ;prep the indp var for each different type
         indp_var = float(ROTATE(*project.imndArray[CI].inversion_Time_ptr, 2))

         indp_var = indp_var / 1000.0
         ;update_status_bar, 'processing image'
         progressbar = Obj_New('progressbar', Color='red', Text='T1 Image fit',/NOCANCEL)
         progressbar -> Start

         for xx=0 , sz_image[1]-1 do begin
          for yy=0 , sz_image[2]-1 do begin
             
             if  abs((*p_image)[xx,yy,0]) gt image_fit_threshold then BEGIN

                     dep_var = abs(ROTATE( reform((*p_image)[xx,yy,*]) , 2))
                     InvRecovSE, dep_var, indp_var , A , chisquare
                     S0[xx,yy] = A[0]
                     T1[xx,yy] = 1.0/A[1]
             end

             endfor

             progressBar -> Update, (float(xx)/float((sz_image[2]-1)))*100.0
           endfor
         progressbar -> Destroy

         ;store the images of S0 and T1
         p_image = ptr_new(fltarr(2,(size((*p_image)))[1],(size((*p_image)))[2]))
         (*p_image)[0,*,*] = S0
         (*p_image)[1,*,*] = T1
         project.dataArray[ci].imgFit1 = p_image

         project.procPramArray[CI].image_fit_flag = 1

       endif else begin

         S0 = reform((*project.dataArray[ci].imgFit1)[0,*,*])
         T1 = reform((*project.dataArray[ci].imgFit1)[1,*,*])

       end

       IMAGE_STATISTICS , S0, MEAN=S0mean, MINIMUM=S0imageMin , MAXIMUM=S0imageMax , STDDEV=S0stddev
       IMAGE_STATISTICS , T1, MEAN=T1mean, MINIMUM=T1imageMin , MAXIMUM=T1imageMax, STDDEV=T1stddev


       sData = [ 'S0image min'+string(S0imageMin),$
              'S0image max'+string(S0imageMax),$
                     'S0mean     '+string(S0mean),$
                 'S0stddev   '+string(S0stddev),$
              'T1image min'+string(T1imageMin),$
                 'T1image max'+string(T1imageMax),$
                 'T1mean     '+string(T1mean),$
                 'T1stddev   '+string(T1stddev)]

       display_stats, sData, 'Fitting parameters'
       
       image_fit_display_image, S0, 'T1-S0 (arbitary floating point values)', $
                                sclmin=S0_min, sclmax=S0_max
       image_fit_display_image, T1, 'T1 (floating point values in secs)', $
                                sclmin=T_min, sclmax=T_max
       
    end

    update_status_bar, ' '

end

; Subroutine name: image_fit_init_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro image_fit_init_event, event

    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets

    widget_control, event.top, get_uvalue=state
    CI = (*state).proj_id
    
    case Event.id of

        (*state).S0_min_slider :BEGIN
          if event.value ge  project.procPramArray[CI].S0_max then begin
              project.procPramArray[CI].S0_min = float(project.procPramArray[CI].S0_max - 0.0001)
              widget_control,(*state).S0_min_slider, SET_VALUE=project.procPramArray[CI].S0_min
          end else begin
              project.procPramArray[CI].S0_min = float(event.value)
          end

          if project.procPramArray[CI].image_fit_auto eq 1 then image_fit, use_ci=CI
          
        end
        
        (*state).S0_max_slider :BEGIN
           if event.value le project.procPramArray[CI].S0_min then begin
              project.procPramArray[CI].S0_max = float(project.procPramArray[CI].S0_min + 0.0001)
              widget_control, (*state).S0_max_slider, SET_VALUE=project.procPramArray[CI].S0_max
          end else begin
              project.procPramArray[CI].S0_max = float(event.value)
          end

          if project.procPramArray[CI].image_fit_auto eq 1 then image_fit, use_ci=CI

        end
        
        (*state).T_min_slider :BEGIN
          if event.value ge project.procPramArray[CI].T_max then begin
              project.procPramArray[CI].T_min = float(project.procPramArray[CI].T_max - 0.0001)
              widget_control, (*state).T_min_slider, SET_VALUE=project.procPramArray[CI].T_min
          end else begin
              project.procPramArray[CI].T_min = float(event.value)
          end

          if project.procPramArray[CI].image_fit_auto eq 1 then image_fit, use_ci=CI
        end

        (*state).T_max_slider :BEGIN
           if event.value le  project.procPramArray[CI].T_min then begin
              project.procPramArray[CI].T_max = float(project.procPramArray[CI].T_min + 0.0001)
              widget_control, (*state).T_max_slider, SET_VALUE=project.procPramArray[CI].T_max
          end else begin
              project.procPramArray[CI].T_max = float(event.value)
          end

          if project.procPramArray[CI].image_fit_auto eq 1 then image_fit, use_ci=CI
        end

        (*state).threshold_slider :BEGIN
            project.procPramArray[CI].image_fit_threshold = float(float(event.value)/100.0)
            project.procPramArray[CI].image_fit_flag = 0

            if project.procPramArray[CI].image_fit_auto eq 1 then image_fit, use_ci=CI
        end
        
        (*state).image_fit_batch:begin

          project.procPramArray[CI].image_fit_batch = event.SELECT
          print, project.procPramArray[CI].image_fit_batch
          
        end
        
        (*state).image_fit_process: begin
         if project.procPramArray[CI].image_fit_batch eq 1 then begin
            print, "Bt2"
            batch_T2_image_fit
         endif else begin
            image_fit, use_ci=CI
         endelse
         
        end
        
        (*state).image_fit_done: begin
            ptr_free, state
            widget_control, event.top, /DESTROY
        end
        
        (*state).image_fit_auto:begin
            project.procPramArray[CI].image_fit_auto = event.SELECT
        end

        else: print, "Unknown Event!"
        
    endcase

end


; Subroutine name: image_fit_init
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro image_fit_init
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci

    S0_min = project.procPramArray[CI].S0_min
    S0_max = project.procPramArray[CI].S0_max
    T_min = project.procPramArray[CI].T_min
    T_max = project.procPramArray[CI].T_max
    image_fit_threshold = project.procPramArray[CI].image_fit_threshold
    image_fit_type = project.procPramArray[CI].image_fit_type

    if image_fit_type eq 0 then begin
       title = 'Image Fit - T1 - Scan #'+string(ci+1, format='(I02)')
       S0_min_name = 'Set S0 min'
       S0_max_name = 'Set S0 max'
       T_min_name = 'Set T1 min (sec)'
       T_max_name = 'Set T1 max (sec)'
    endif else if image_fit_type eq 1 then begin
       title = 'Image Fit - T1 - Scan #'+string(ci+1, format='(I02)')
       S0_min_name = 'Set S0 min'
       S0_max_name = 'Set S0 max'
       T_min_name = 'Set T1 min (sec)'
       T_max_name = 'Set T1 max (sec)'
    endif else if image_fit_type eq 2 then begin
       title = 'Image Fit - T2 - Scan #'+string(ci+1, format='(I02)')
       S0_min_name = 'Set S0 min'
       S0_max_name = 'Set S0 max'
       T_min_name = 'Set T2 min (sec)'
       T_max_name = 'Set T2 max (sec)'
    endif else if image_fit_type eq 3 then begin
       title = 'Image Fit - ADC (exp) - Scan #'+string(ci+1, format='(I02)')
       S0_min_name = 'Set S0 min'
       S0_max_name = 'Set S0 max'
       T_min_name = 'Set ADC min (um^2/ms)'
       T_max_name = 'Set ADC max (um^2/ms)'
    endif else if image_fit_type eq 4 then begin
       title = 'Image Fit - ADC (str exp) - Scan #'+string(ci+1, format='(I02)')
       S0_min_name = 'Set S0 min'
       S0_max_name = 'Set S0 max'
       T_min_name = 'Set ADC min (um^2/ms)'
       T_max_name = 'Set ADC max (um^2/ms)'
    endif else if image_fit_type eq 5 then begin
       title = 'Image Fit - T1 - Scan #'+string(ci+1, format='(I02)')
       S0_min_name = 'Set S0 min'
       S0_max_name = 'Set S0 max'
       T_min_name = 'Set T1 min (sec)'
       T_max_name = 'Set T1 max (sec)'
    endif

    image_fit_base = widget_base(/COLUMN, TITLE=title)
    
    disp_scl_base = widget_base(image_fit_base, /column, /align_center, frame=1)
    
    lbl = widget_label(disp_scl_base, value="Image Display Scaling", /align_center)
    ROW1 = widget_base(disp_scl_base, /ROW , FRAME=0)

    S0_min_slider = CW_FSLIDER(ROW1, MINIMUM=0, FRAME=0, MAXIMUM=50, $
                           TITLE = S0_min_name ,  /EDIT , $
                           UNAME='S0_min_slider', $
                           VALUE= S0_min )

    S0_max_slider = CW_FSLIDER(ROW1, MINIMUM=0, FRAME=0, MAXIMUM=500, $
                           TITLE = S0_max_name ,   /EDIT , $
                           UNAME='S0_max_slider', $
                           VALUE=S0_max  )
    ROW2 = widget_base( disp_scl_base, /ROW , FRAME=0)


    T_min_slider = CW_FSLIDER(ROW2, MINIMUM=0, FRAME=0, MAXIMUM=10, $
                           TITLE =T_min_name  ,    /EDIT , $
                           UNAME='T_min_slider', $
                           VALUE= T_min )

    T_max_slider = CW_FSLIDER (ROW2, MINIMUM=0, FRAME=0, MAXIMUM=10, $
                           TITLE = T_max_name ,    /EDIT , $
                           UNAME='T_max_slider', $
                           VALUE=T_max  )
    ROW3 = widget_base( image_fit_base, /ROW , FRAME=1)

    threshold_slider = widget_slider(ROW3, MINIMUM=0, FRAME=0, MAXIMUM=100, $
                           TITLE = 'Threshold % cutoff' ,  SCR_XSIZE=130 , $
                           UNAME='threshold_slider', $
                           VALUE= image_fit_threshold*100)

    COLUMN1 = widget_base( ROW3, /COLUMN , FRAME=0 )

    image_fit_process= widget_button(COLUMN1,value='Process'$
                 ,UNAME = 'image_fit_process')

    image_fit_done= widget_button(COLUMN1,value='Done'$
                 ,UNAME = 'image_fit_done')

    exb = widget_base( ROW3 ,/NONEXCLUSIVE)

    image_fit_auto= widget_button(exb ,value='Auto Update' ,UNAME = 'image_fit_auto')
    if project.procPramArray[CI].image_fit_auto eq 1 then WIDGET_CONTROL, /SET_BUTTON, image_fit_auto

    batch_sens = 0
    if image_fit_type eq 2 then batch_sens = 1


    image_fit_batch = widget_button(exb, value='Batch T2 Export', UNAME = 'image_fit_batch', SENSITIVE = batch_sens)

    state = { proj_id: ci, $
              tlb: image_fit_base, $
              S0_min_slider: S0_min_slider, $
              S0_max_slider: S0_max_slider, $
              T_min_slider: T_min_slider, $
              T_max_slider: T_max_slider, $
              threshold_slider: threshold_slider, $
              image_fit_process: image_fit_process, $
              image_fit_done: image_fit_done, $
              image_fit_auto: image_fit_auto, $
              image_fit_batch: image_fit_batch }
              
    widget_control, image_fit_base, /realize, set_uvalue=ptr_new(state, /no_copy)

    xmanager, 'image_fit_init',image_fit_base,/no_block

end


; Subroutine name: mas_curvefit
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; parameters:
;   type:
;     1 for T1 progressive saturation ROI
;     2 for T1 progressive saturation IMAGE
;     3 for T1 Inversion Recovery ROI
;     4 for T1 Inversion Recovery IMAGE
;     5 for T2 Multi Spin Echo ROI
;     6 for T2 Multi Spin Echo IMAGE
;     7 for ADC ROI
;     8 for ADT Image
; Editing Information:
    ;Edited by HS 2006/10/09
    ;Fix spelling mistakes and commenting.

pro mas_curvefit, type, event
    COMPILE_OPT IDL2

    ;print,'curve fit'

    COMMON scan_data
    CI = project.ci

    ;chk that there are mulitple TR's or TE's in the scan selected
    if (type eq 1 or type eq 2 or type eq 3 or type eq 4 or type eq 9 or type eq 10) then begin
       if not(PTR_VALID(project.imndArray[CI].rep_Time_ptr)) then begin
         update_status_bar, 'Scan selected does not have multiple TR`s'
         return
       end
    end else  if (type eq 5 or type eq 6 ) then begin
       if not(PTR_VALID(project.imndArray[CI].echo_Time_ptr)) then begin
         update_status_bar, 'Scan selected does not have multiple TE`s'
         return
       end
    end else if (type eq 7 or type eq 8 or type eq 13 or type eq 14 ) then begin
       if not(PTR_VALID(project.imndArray[CI].bval_Array)) then begin
         update_status_bar, 'Scan selected does not have multiple B-Values'
         return
       end
    end else if (type eq 11 or type eq 12) then begin
       if not(PTR_VALID(project.imndarray[CI].inversion_Time_ptr)) then begin
         update_status_bar, 'Scan selected does not have multiple inversion times.'
         return
       end
    endif else return


    if type eq 1 then begin
       update_status_bar, 'T1 Saturation Recovery GE ROI'
       project.procPramArray[CI].ROI_fit_type = 0
       roi_Curve_Fit, event

    end else if type eq 2 then begin
       update_status_bar, 'T1 Saturation Recovery GE IMAGE'
       project.procPramArray[CI].image_fit_type = 0
       image_fit_init
       
    end else if type eq 3 then begin
       update_status_bar, 'Unsupported Curvefit: T1 Inversion Recovery GE ROI'
       void = dialog_message('T1 Inversion Recovery GE ROI is not supported at this time.', $
                             /center, /error)
    end else if type eq 4 then begin
       update_status_bar, 'Unsupported Curvefit: T1 Inversion Recovery GE IMAGE'
       void = dialog_message('T1 Inversion Recovery GE IMAGE is not supported at this time.', $
                             /center, /error)

    end else if type eq 5 then begin
       update_status_bar, 'T2 Multi Spin Echo ROI'
       project.procPramArray[CI].ROI_fit_type = 2
       roi_Curve_Fit, event

    end else if type eq 6 then begin
       update_status_bar, 'T2 Multi Spin Echo IMAGE'
       project.procPramArray[CI].image_fit_type = 2
       image_fit_init
       
    end else if type eq 7 then begin
        update_status_bar, 'ADC ROI (STANDARD S=S0*exp(-bd))'
       project.procPramArray[CI].ROI_fit_type = 3
       roi_Curve_Fit, event

    end else if type eq 8 then begin
       update_status_bar, 'ADC IMAGE (STANDARD S=S0*exp(-bd))'
       project.procPramArray[CI].image_fit_type = 3
       image_fit_init
       
    end else if type eq 9 then begin
       update_status_bar, 'T1 Saturation Recovery SE ROI'
       project.procPramArray[CI].ROI_fit_type = 1
       roi_Curve_Fit, event
    
     end else if type eq 10 then begin
       update_status_bar, 'T1 Saturation Recovery SE Image'
       project.procPramArray[CI].image_fit_type = 1
       image_fit_init
       
    end else if type eq 11 then begin
       update_status_bar, 'Unsupported Curvefit: T1 Inversion Recovery SE ROI'
       ;return
       project.procPramArray[CI].ROI_fit_type = 5
       update_status_bar, 'T1 Inversion Recovery SE ROI'
       roi_Curve_Fit, event

    end else if type eq 12 then begin
       update_status_bar, 'Unsupported Curvefit: T1 Inversion Recovery SE IMAGE'
       ;return
       project.procpramarray[CI].image_fit_type = 5
       image_fit_init
    end else if type eq 13 then begin
       update_status_bar, 'ADC ROI (EXPERIMENTAL S=S0*exp(-(bd)^alpha))'
       project.procPramArray[CI].ROI_fit_type = 4
       roi_Curve_Fit, event

    end else if type eq 14 then begin
       update_status_bar, 'ADC IMAGE (EXPERIMENTAL S=S0*exp(-(bd)^alpha)'
       project.procPramArray[CI].image_fit_type = 4
       image_fit_init
    
    end

END

pro mas_curvefit_pd_align, state, best_tx=best_tx

    common scan_data
    common _mc_working_set, working_set, _mc_working_set
    ;; this is a common block to share current working
    ;; image data with the objective function. Pointers
    ;; should reduce memory usage and keep the footprint
    ;; of the common block to a mimimum

    T1_ind = (*state).selected_T1
    T2_ind = (*state).selected_T2
    if (T1_ind lt 0 or T2_ind lt 0) then return
    
    T1_S0_data = reform((*project.dataarray[T1_ind].imgfit1)[0,*,*])
    T2_S0_data = reform((*project.dataarray[T2_ind].imgfit1)[0,*,*])
    
    binsize = 4
    
    working_set = { ref_image: T1_S0_data, $
                    tar_image: T2_S0_data, $
                    bs:        binsize,   $
                    normalized: 0 }
    max_mi     = 0
    best_mi    = 0
    scls       = [0,-2,2,4,-4,-6,6]
    start_pt   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    best_ever  = start_pt
    best_shift = [ start_pt[0], start_pt[1] ]
    
    use_trl = 1
    use_rot = 0
    use_shr = 0
    use_scl = 0 ;float(project.procpramarray[0].mc_mi_use_dilation)
    
    tol = 1e-3
    
    optvec = [use_trl, use_trl, use_shr, use_shr, use_rot, use_scl, use_scl]
    
    print, "******** Using Cascade Method ********"
    dir_vec = fltarr(n_elements(optvec))
    iterations = total([use_trl, use_shr, use_rot])
    n_iter = 0
    
    for tx = 0, 3-1 do begin ;;iterations-1 do begin
    
        if (([use_trl, use_shr, use_rot])[tx] eq 0) then continue
        
        dir_vec_wk = dir_vec
        case tx of
            0: dir_vec_wk[0:1] = optvec[0:1]
            1: dir_vec_wk[2:3] = optvec[2:3]
            2: dir_vec_wk[4]   = optvec[4]
        endcase
        
        dir_vec_wk = diag_matrix(dir_vec_wk)
        
        powell, start_pt, dir_vec_wk, tol, max_mi,$
            'mc_objective_function', ITER=n_iter
            
        best_mi = -max_mi
        best_ever = start_pt
        
        print, start_pt
    endfor

    best_tx = best_ever
    
end


pro mas_curvefit_pd_init_event, event

    common scan_data, project
    forward_function mc_transform_image
    
    widget_control, event.top, get_uvalue=state
    
    if (not ptr_valid(state)) then return
    
    ui = (*state).ui
    
    case event.id of

        ui.btn_select_T1: begin
            ind = widget_info(ui.lst_scans, /list_select)
            if (ind lt 0) then return
            
            if (not ptr_valid(project.imndarray[ind].rep_time_ptr) and not ptr_valid(project.imndarray[ind].inversion_time_ptr)) then begin
                void = dialog_message(['This does not appear to be a variable TR dataset.'], $
                                      /error, /center)
                return
            endif else if (not ptr_valid(project.dataarray[ind].imgfit1)) then begin
                void = dialog_message(['The T1 fit for this scan has not been performed', $
                                        'Please perform the T1 fit and then try again'], $
                                        /error, /center)
                return
            endif
            
            (*state).selected_T1 = widget_info(ui.lst_scans, /list_select)
            
            slice = project.procpramarray[ind].sdim_start
            widget_control, ui.txt_slice, set_value=string(slice, format='(I0)')
            if (ptr_valid(project.imndarray[ind].inversion_time_ptr)) then begin
                TR = project.imndarray[ind].recov_time
            endif else begin
            TR = max(*project.imndarray[ind].rep_time_ptr, max_tr_ind) ;; in ms
            endelse
            widget_control, ui.txt_TR, set_value=string(TR, format='(G0)')
            TE = project.imndarray[ind].echo_time ;; in ms
            widget_control, ui.txt_TE, set_value=string(TE, format='(G0)')
            
        end
        
        ui.btn_select_T2: begin
            ind = widget_info(ui.lst_scans, /list_select)
            if (ind lt 0) then return
            
            if (not ptr_valid(project.imndarray[ind].echo_time_ptr)) then begin
                void = dialog_message(['This does not appear to be a variable TE dataset.'], $
                                      /error, /center)
                return
            endif else if (not ptr_valid(project.dataarray[ind].imgfit1)) then begin
                void = dialog_message(['The T2 fit for this scan has not been performed', $
                                        'Please perform the T2 fit and then try again'], $
                                        /error, /center)
                return
            endif
             
            (*state).selected_T2 = widget_info(ui.lst_scans, /list_select)
                     
        end
        ui.btn_cancel: begin
            ptr_free, state
            widget_control, event.top, /destroy
            return
            
        end
        
        ui.btn_compute: begin
        
            T1_ind = (*state).selected_T1
            T2_ind = (*state).selected_T2
            if (T1_ind lt 0 or T2_ind lt 0) then return
            
            slice = project.procpramarray[T1_ind].sdim_start

            T1_data = reform((*project.dataarray[T1_ind].imgfit1)[1,*,*])
            T2_data = reform((*project.dataarray[T2_ind].imgfit1)[1,*,*])

            if (widget_info(ui.btn_align, /button_set)) then begin
                mask = reform(long(T2_data ne 0), (size(T2_data, /dim))[0], (size(T2_data, /dim))[1]) * 1000
                mask *= 1000.0
             
                mas_curvefit_pd_align, state, best_tx=best_tx
                
                mask = mc_transform_image(float(mask), best_tx, scale_first=0)
                mask = round(mask/max(mask)) < 1.0
                T2_data = mc_transform_image(T2_data, best_tx, scale_first=0) * mask
            endif
            
            TE = project.imndarray[T1_ind].echo_time ;; in ms
            if (ptr_valid(project.imndarray[T1_ind].inversion_time_ptr)) then begin
                TI_max = max(*project.imndarray[T1_ind].inversion_time_ptr, max_ti_ind)
                TR = project.imndarray[T1_ind].recov_time
                SIG_index = max_ti_ind
            endif else begin 
                TR = max(*project.imndarray[T1_ind].rep_time_ptr, max_tr_ind) ;; in ms
            SIG_index = max_tr_ind
            endelse
            
            SIG = reform((*project.dataarray[T1_ind].state1)[*,*,slice,SIG_index])
            
            print, string(TR, TE, format='(%"TR = %0.4f ms; TE = %0.4f ms")')
            
            TR /= 1000.0
            TE /= 1000.0
            PD = SIG/( (1.0-exp(-TR/T1_data))*(exp(-TE/T2_data)) )
            
            infnan = where(1.0-finite(PD), num_infnan)
            if (num_infnan gt 0) then begin
                PD[infnan] = 0.0
            endif
        
            if (project.iImage_flag eq 0) then begin
            
                p_PD = ptr_new(PD)

                mas_rotate_flip, p_PD
                mas_zoom, p_PD
        
                PD_nonneg = where(PD gt 0, n_PD_nonneg)
                scl_mean  = mean(PD[PD_nonneg])
                scl_std   = stddev(PD[PD_nonneg])
                
                widget_control, ui.txt_scl_min, get_value=txt_scl_min
                widget_control, ui.txt_scl_max, get_value=txt_scl_max
                scl_min = float(txt_scl_min[0])
                scl_max = float(txt_scl_max[0])
                PD_scl = bytscl(*p_PD, min=scl_min, max=scl_max)
                        
                PD = *p_PD & ptr_free, p_PD
                
                mas_display_multi, PD_scl, $
                                   tab_title=title,  $
                                   /standalone,$
                                   fov_x=xfov, $
                                   fov_y=yfov, $
                                   fp_vals=PD, $
                                   show_axis_labels=0
            endif else begin
            
                iimage, PD, /block

            endelse

            ;ptr_free, state
            ;widget_control, event.top, /destroy
            
            return
            
        end
        
        else: return

    endcase
    
    if ((*state).selected_T1 ge 0 and (*state).selected_T2 ge 0) then begin
        widget_control, ui.btn_compute, sensitive=1
    endif else begin
        widget_control, ui.btn_compute, sensitive=0
    endelse
    
    temp = (*state).scan_titles
    for i = 0, n_elements(temp)-1 do begin
        if (i eq (*state).selected_T2) then begin   
            temp[i] = '[T2] '+temp[i]
        endif else if (i eq (*state).selected_T1) then begin
            temp[i] = '[T1] '+temp[i]
        endif else begin
            temp[i] = '     '+temp[i]
        endelse
    endfor
    
    widget_control, ui.lst_scans, set_value = temp
    
    
end

pro mas_curvefit_pd_init

    common scan_data
    common common_widgets
    
    num_opened_scans = project.ni
    scan_titles = strarr(num_opened_scans)
    
    for i = 0, num_opened_scans-1 do begin
    
        display_name = project.imndArray[i].display_Name
        if (display_name eq '') then begin
            display_name = project.imndArray[i].file_Path
        endif
        
        scan_titles[i] = display_name
    
    endfor
    
    base = widget_base(title='Curve Fit PD', /column, group_leader=WID_BASE_MAIN);, /modal)
    
    lbl = widget_label(base, value='Please select T1/T2 fit data sets:', /align_center)
    
    lst_scans = widget_list(base, value='     '+scan_titles, $
                            ysize=10, $
                            xsize=5+min([80,max(strlen(scan_titles))]))

    base_info = widget_base(base, /row, /align_center, /base_align_center)
    lbl = widget_label(base_info, value='TR (ms):')
    txt_TR = widget_text(base_info, value='', xsize=5)

    lbl = widget_label(base_info, value='TE (ms):')
    txt_TE = widget_text(base_info, value='', xsize=5)

    lbl = widget_label(base_info, value='Slice:')
    txt_slice = widget_text(base_info, value='', xsize=4)

    align_base = widget_base(base, /row, /nonexclusive, /align_center, /base_align_center)
    btn_align = widget_button(align_base, value='Attempt to align data sets')
    
    scl_base = widget_base(base, /row, /align_center, /base_align_center)
    lbl = widget_label(scl_base, value='Display scale min:')
    txt_scl_min = widget_text(scl_base, value='0', xsize=6, /editable)
    lbl = widget_label(scl_base, value=' max:')
    txt_scl_max = widget_text(scl_base, value='500', xsize=6, /editable)
    
    btn_base = widget_base(base, /row, /align_center, /base_align_center)
    btn_select_T1 = widget_button(btn_base, value="Select T1 Fit")
    btn_select_T2 = widget_button(btn_base, value="Select T2 Fit")
    btn_compute   = widget_button(btn_base, value="Compute", sensitive=0)
    btn_cancel   = widget_button(btn_base, value="Cancel")
    
    state = { ui: { lst_scans: lst_scans, $
                    btn_select_T1: btn_select_T1, $
                    btn_select_T2: btn_select_T2, $
                    btn_align: btn_align, $
                    btn_compute: btn_compute, $
                    btn_cancel: btn_cancel, $
                    txt_scl_min: txt_scl_min, $
                    txt_scl_max: txt_scl_max, $
                    txt_TR:txt_TR, $
                    txt_TE:txt_TE, $
                    txt_slice: txt_slice }, $
              scan_titles: scan_titles, $
              selected_T1: -1L, $
              selected_T2: -1L }
    
    widget_control, base, /realize, set_uvalue=ptr_new(state)
    
    xmanager, 'mas_curvefit_pd_init', base, /no_block
    
end
