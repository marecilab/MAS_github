;; $Id$
;;

;; Subroutine name: mas_nnls_QR_solve
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    A:  Input Matrix
;;    b:  Right-hand Side
;;    S:  Desired S value
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Uses QR factorization to solve a linear problem.
;;
;; Editing Information:

pro mas_nnls_QR_solve, A, b, S, double=double

    mas_nnls_QR, A, Q, R, /double

    QTb = transpose(Q) ## b

    S = invert(R) ## QTb

end

;; Subroutine name: mas_nnls_QR
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    A:  Input Matrix
;;    Q:  Orthogonal Matrix Q result
;;    R:  Upper-triangular R result
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Factorized a matrix into an orthogonal matrix Q and a upper
;;   triangular matrix R
;;
;; Editing Information:

pro mas_nnls_QR, A, Q, R, double=double

    if (keyword_set(double)) then begin
       A_mat = double(A)
    endif else begin
       A_mat = A
    endelse

    A_sz = size(A_mat, /dimensions) ; => A_sz[col,row]
    
    col_placeholder = fltarr(1,A_sz[1])
    I_n = diag_matrix(replicate(1D, A_sz[1]))
    P_cum = I_n
    for col = 0, A_sz[0]-1 do begin
       
       x = A_mat[col,col:*]
       
       u = x
       u[0] -= sqrt(total(x*x))
       P_sub = I_n[col:*,col:*] - (2D/total(u*u)) * (u ## transpose(u))
       P = I_n
       P[col:*,col:*] = P_sub
       P_cum = P_cum ## P

       A_mat[col:A_sz[0]-1,col:A_sz[1]-1] = (P_sub ## A_mat[col:A_sz[0]-1,col:A_sz[1]-1])
       
    endfor

    Q = P_cum
    R = A_mat

end

;; Subroutine name: mas_nnls_pinv
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    a_mat:   The input matrix
;;  tik_mat:   A diagonal tikhonov matrix for regularization 
;;  use_tik:   Set this keyword to instruct the procedure to look for
;;             a tikhonov matrix
;;check_inv:   Set this keyword to display some statistics about the
;;             quality of the pseudoinverse
;;
;; Bugs or Important Comments to Developers:
;;
;;    This would be better handled in another way. Maybe using SVD
;;
;; Purpose of subroutine:
;;
;;    Compute a "pseudoinverse" of a non-square matrix.
;;
;; Editing Information:

FUNCTION mas_nnls_pinv, a_mat, tik_mat=tik_mat, use_tik=use_tik, $
                        check_inv=check_inv, DOUBLE=DOUBLE

    if (keyword_set(use_tik)) then begin
       DLS_mat = invert( (matrix_multiply(a_mat, a_mat, /btranspose) + tik_mat) ) # a_mat
    endif else begin
       DLS_mat = invert( matrix_multiply(a_mat, a_mat, /btranspose) ) # a_mat
    endelse

    IF ( keyword_set(check_inv) NE 0 ) THEN BEGIN
        
        IF ( !D.NAME EQ 'X' ) THEN BEGIN
            
            WINDOW, /FREE, TITLE = 'SVD_MATRIX_INVERT check'
            SURFACE, transpose(dls_mat) ## a_mat, /LEGO, $
              TITLE = 'Check of SVD matrix inversion', $
              XTITLE = 'Column index', YTITLE = 'Row index', $
              CHARSIZE = 1.5
            
        ENDIF
        
    ENDIF
    
    return, transpose(DLS_mat)
    
end

;; Subroutine name: mas_nnls
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    Given, Ex = f. Compute the non-negative least squares solution
;;    x.
;;
;;    resid:  Set this to a named variable to receive the residual
;;    vector
;;
;;    resnorm: Set this to a named variable to receive the size of the
;;    residual vector
;;
;; Bugs or Important Comments to Developers:
;;
;;
;;
;; Purpose of subroutine:
;;
;;    Solve the linear equation Ex = f for x in a non-negative sense.
;;
;; Editing Information:

pro mas_nnls, E, f, x, resid=resid, resnorm=resnorm

    sz_E = size(E, /dimensions)
    m = sz_E[0]
    n = sz_E[1]

    PPP = lonarr(n)
    ZZZ = lonarr(n) + 1
    PP  = where(PPP eq 0)
    ZZ  = where(ZZZ ne 0)

    E_PP = dblarr(m,n)

    f = double(f)
    x = transpose(dblarr(n))

    resid  = f - matrix_multiply(E, x, /btranspose)
    w      = E ## resid
    tol    = 1e-5
    iter   = 0
    itmax  = fix(n)
    check_inv = 0

    lambda = 1e-7 ;;0.0000000001D;0.01
    tik_mat = diag_matrix(replicate(lambda,m))

    while (total(ZZZ) ne 0 && total(w[ZZ] gt tol) gt 0) do begin

        t = where(w eq max(w[ZZ]))
        PPP[t] = 1
        ZZZ[t] = 0
        PP = where(PPP ne 0, ctp)
        ZZ = where(ZZZ ne 0, ctz)

        if ctp eq 0 or ctz eq 0 then begin
            ;;print, 'mas_nnls: no indices left'
            break
        endif

        E_PP[*,PP] = E[*,PP]
        E_PP[*,ZZ] = 0
        
        ;z = LA_LEAST_SQUARES(transpose(E_PP), f, method=3, status=st, /double)
        ;if (st ne 0) then begin
            z = mas_nnls_pinv(E_PP, /double, tik_mat=tik_mat, /use_tik) # f
        ;endif
        z[ZZ] = 0

        while ( total(z[PP] lt tol) ) do begin
            
            iter++
            if iter gt itmax then begin
                ;;print, 'mas_nnls: iteration count is exceeded, exiting loop.'
                exitflag = 0
                resnorm = sqrt(total(resid * resid))
                x = z
                return
            endif
        
            QQ = where( logical_and((z lt tol), PPP) ne 0 )
            alpha = min(x[QQ]/(x[QQ] - z[QQ]))
            x = x + (alpha # (z - x))
            
            ij = ( logical_and((abs(x) lt tol), PPP) )
            ;ij = ((abs(x) lt tol) and PPP) 
            ij = where(ij ne 0, ctij)
            if (ctij eq 0) then begin
               ZZZ[*] = 0
               PPP[*] = 1
            endif else begin
               ZZZ[ij] = 1
               PPP[ij] = 0
            endelse

            ZZ = where(ZZZ ne 0, ctp)
            PP = where(PPP ne 0, ctz)

            if (ctz gt 0) then E_PP[*,PP] = E[*,PP]
            if (ctp gt 0) then E_PP[*,ZZ] = 0
            
            ;z = LA_LEAST_SQUARES(transpose(E_PP), f, method=3, /double, status=st)
            ;if (st ne 0) then begin
                z = mas_nnls_pinv(E_PP, /double, tik_mat=tik_mat, /use_tik) # f
            ;endif
            if (ctp gt 0) then z[ZZ] = 0
            
        endwhile

        x = z
        resid = f - (E # x)
        w = transpose(E) # resid
        
    endwhile

    resnorm = sqrt(total(resid * resid))

end

