#+nil ; for compiling this file in ecl into a standalone binary
(progn
  (compile-file "cholesky.lisp" :output-file "cholesky.o" :system-p t)
  (c::build-program "cholesky" :lisp-files '("cholesky.o")))



(declaim (optimize (safety 3) (debug 3) (speed 1)))
;(declaim (optimize (safety 0) (debug 0) (speed 3)))

(deftype fun1 (args retval) 
  `(function ,args (values ,retval &optional)))

(defparameter *s*
  '((4 4 4 0 0 0 0 0)
    (4 5 2 0 0 0 0 0)
    (4 2 9 1 2 1 0 0)
    (0 0 1 10 -7 1 0 0)
    (0 0 2 -7 14 0 2 1)
    (0 0 1 1 0 9 -8 -2)
    (0 0 0 0 2 -8 9 0)
    (0 0 0 0 1 -2 0 6)))

(deftype mat ()
  `(simple-array single-float 2))

(deftype vec ()
  `(simple-array single-float 1))


(defmacro make-vec (n &rest rest)
  `(make-array ,n :element-type 'single-float
	       ,@rest))

(defmacro make-mat (r &rest rest)
  `(make-array (list ,r ,r) :element-type 'single-float
	       :initial-element 0s0
	       ,@rest))

(defun afloat1 (list)
  (make-vec (length list)
	    :initial-contents
	    (mapcar #'(lambda (x) (* 1s0 x)) list)))

(defun afloat2 (lists)
  (make-array (list (length lists) (length (first lists)))
	      :element-type 'single-float
	      :initial-contents 
	      (mapcar #'(lambda (row)
			  (mapcar #'(lambda (x) (* 1s0 x)) row))
		      lists)))


(defmacro with-arrays (arrays &body body)
  `(macrolet ,(mapcar 
	       (lambda (array)
		 `(,array (&rest indices) 
			  `(aref ,',array ,@indices)))
	       arrays)
     ,@body))


(declaim (ftype (fun1 (mat) mat) cholesky))
(defun cholesky (s)
  (destructuring-bind (r c) (array-dimensions s)
    (assert (= r c))
    (let ((l (make-mat r)))
      (with-arrays (l s)
	(dotimes (i r)
	 (dotimes (k (+ i 1))
	   (setf (l i k)
		 (if (= i k)
		     (sqrt (- (s i i)
			      (loop for j below i
				 sum (expt (l i j) 2))))
		     (/ (- (s i k)
			   (loop for j below i
			      sum (* (l i j)
				     (l k j))))
			(l k k)))))))
      l)))


#+nil
(cholesky (afloat2 *s*))

(defparameter *yt*  '(4 1 10 -18 16 8 -9 1))

(declaim (ftype (fun1 (mat vec) vec) forward-eliminate))
(defun forward-eliminate (l b)
  "Solve l*y=b for y with lower triangular matrix L."
  (let* ((n (length b))
	 (y (make-vec n)))
    (with-arrays (l y b)
      (dotimes (i n)
	(setf (y i) (b i))
	(dotimes (j i)
	  (decf (y i) (* (y j) (l i j))))
	(setf (y i) (/ (y i)
		       (l i i)))))
    y))

#+nil
(forward-eliminate (cholesky (afloat2 *s*))
		   (afloat1 *yt*))

(declaim (ftype (fun1 (mat) mat)
		transpose))
(defun transpose (a)
  (destructuring-bind (r c) (array-dimensions a)
   (let ((b (make-array (reverse (array-dimensions a))
			:element-type 'single-float)))
     (with-arrays (a b)
       (dotimes (i r)
	 (dotimes (j c)
	   (setf (b j i) (a i j)))))
     b)))

#+nil
(transpose (afloat2 *s*))

(declaim (ftype (fun1 (mat vec) vec) back-substitute))
(defun back-substitute (u b)
  "Solve u*y=b for y with upper triangular matrix U."
  (let* ((n (length b))
	 (y (make-vec n)))
    (with-arrays (u y b)
      (loop for i from (1- n) downto 0 do
	   (setf (y i) (b i))
	   (loop for j from (1- n) above i do
		(decf (y i) (* (y j) (u i j))))
	   (setf (y i) (/ (y i)
			  (u i i)))))
    y))

#+nil
(let ((l (cholesky (afloat2 *s*))))
  (back-substitute (transpose l)
		   (forward-eliminate l
				      (afloat1 *yt*))))
(declaim (ftype (fun1 (mat vec) vec) solve))
(defun solve (m y)
  "Solve M x = y for x."
  (let ((l (cholesky m)))
    (back-substitute (transpose l)
		     (forward-eliminate l
					y))))

#+nil
(solve (afloat2 *s*) (afloat1 *yt*))

(defparameter *problem*
  '((:sending-end-voltage . 2s0)
    (:resistivity . 1s0) ; per unit-length
    (:conductivity . .5s0) ; per unit-length
    (:segments . ((3 2))))) ; another element with length .5
#+nil
(assoc :segments *problem*)

(defun def-problem (v r g segs)
  (setf *problem*
	`((:sending-end-voltage . ,v)
	  (:resistivity . ,r) ; per unit-length
	  (:conductivity . ,g) ; per unit-length
	  (:segments . ,segs))))

(defun get-param (name)
  (rest (assoc name *problem*)))
#+nil
(get-param :segments)

(defun setup ()
   (let* ((segs (get-param :segments))	; elements with length
	  (x 0)
	  (pos (let ((res ()))
		 (mapcar #'(lambda (nl)
			     (destructuring-bind (num len) nl
			       (dotimes (i num)
				 (push (incf x (/ len num))
				       res))))
			 segs)
		 (reverse res))))
     (push 0s0 pos)
     (afloat1 pos)))
#+nil 
(setup)


(defun ndis (ncon)
  (* 2 (- ncon 1)))
#+nil
(ndis (length (setup)))
(defun make-connection-matrix (ncon)
  (let* ((ndis (ndis ncon))
	 (c (make-array (list ndis ncon)
			:element-type 'single-float
			:initial-element 0s0)))
    (with-arrays (c)
      (dotimes (i ncon)
	(setf (c i i) 1s0))
      (dotimes (i (- ndis ncon))
	(setf (c (+ ncon i) (1+ i)) 1s0))
      c)))
#+nil
(print-binary-matrix
 (make-connection-matrix (length (setup))))

(defun make-connected-coefficient-matrix2 ()
  (let* ((s (setup))
	 (n (length s))
	 (md (make-disjoint-coefficient-matrix s))
	 (m (make-mat n)))
    (with-arrays (m md)
      (dotimes (i n)
	(dotimes (j n)
	  (cond ((= i j)
		 (incf (m i j) (md i j)))
		((= (+ n i -1) j)
		 (incf (m i j) (md i j)))))))
    m))



(defun print-binary-matrix (m)
  (destructuring-bind (r c) (array-dimensions m)
    (dotimes (i r)
      (dotimes (j c)
	(if (< (abs (aref m i j)) 0.1s0)
	    (format t "_")
	    (format t "*")))
      (terpri))))

(defun print-matrix (m)
  (destructuring-bind (r c) (array-dimensions m)
    (dotimes (i r)
      (dotimes (j c)
	(if (< -0.01s0 (aref m i j) 0.01s0)
	    (format t "______")
	    (format t "~5@f " (aref m i j))
	    ))
      (terpri))))

(declaim (ftype (fun1 (vec) mat) make-disjoint-coefficient-matrix))
(defun make-disjoint-coefficient-matrix (x)
  (let* ((ncon (length x))
	 (ndis (ndis ncon))
	 (c (make-mat ndis))
	 (co (get-param :conductivity))
	 (re (get-param :resistivity)))
    (with-arrays (x c)
     (loop for i below (- ndis 1) by 2 do
	(let* ((j (1+ i))
	       (p (floor i 2))
	       (l (- (x (1+ p)) 
		     (x p))))
	  (when (< l 0s0)
	    (error "Non-positive element length."))
	  (let ((diag (+ (/ (* re l))
			 (* 1/3 co l)))
		(off (+ (/ (* -1 re l))
			(* 1/6 co l))))
	    (setf (c i i) diag
		  (c i j) off
		  (c j i) off
		  (c j j) diag
		  )))))
    c))


#+nil
(print-binary-matrix
 (make-disjoint-coefficient-matrix (setup)))

#+nil
(print-matrix
 (make-disjoint-coefficient-matrix (setup)))




;; their code
;; M_kl = c_ki Md_ij c_lj
;;   oo     o*    **   o*
;; * .. many, o .. not many



(declaim (ftype (fun1 (mat mat) mat) 
		make-connected-coefficient-matrix3))
(defun make-connected-coefficient-matrix3 (c m)
  (destructuring-bind (mr mcol) (array-dimensions m)
   (destructuring-bind (r col) (array-dimensions c)
     (assert (= mr mcol col))
     (let ((mc (make-mat r)))
       (with-arrays (mc c m)
	 (dotimes (k r)
	   (dotimes (l r)
	     (let ((s 0s0))
	       (dotimes (j col)
		 (dotimes (i col)
		   (incf s (* (c k i) (m i j) (c l j)))))
	       (setf (mc k l) s))))
	 mc)))))
;; W=v^T_i M_ik v_k
;;  =w^T_i C_ji M_jk C_kn w_n
;; Mcon_in = C_ji C_kn M_jk
(declaim (ftype (fun1 (mat mat) mat)
		make-connected-coefficient-matrix))
(defun make-connected-coefficient-matrix (c m)
  (destructuring-bind (mz ms) (array-dimensions m)
   (destructuring-bind (z s) (array-dimensions c)
     (assert (= z mz ms))
     (let ((mc (make-mat s)))
       (with-arrays (mc c m)
	 (dotimes (i s)
	   (dotimes (n s)
	     (dotimes (j z)
	       (dotimes (k z)
		 (incf (mc i n) (* (c j i) (c k n) (m j k)))))))
	 mc)))))

;; 0 -- 4  1 -- 5  2 -- 3          ((3 2))
;; 0 -- 5  1 -- 6  2 -- 7  3 -- 4  ((4 2))

(declaim (ftype (fun1 (mat fixnum fixnum) mat)
		swap-rows))
(defun swap-rows (m r q)
  (let* ((n (array-dimension m 0))
	 (b (make-mat n)))
    (with-arrays (m b)
      (dotimes (i n)
	(cond 
	  ((= i r) (dotimes (j n)
		     (setf (b i j) (m q j))))
	  ((= i q) (dotimes (j n)
		     (setf (b i j) (m r j))))
	  (t (dotimes (j n)
	       (setf (b i j) (m i j)))))))
    b))

(declaim (ftype (fun1 (vec fixnum fixnum) vec)
		swap-cols))
(defun v-swap-cols (v r q)
  (let* ((n (length v))
	 (b (make-vec n)))
    (with-arrays (v b)
      (dotimes (i n)
	(cond ((= i r) (setf (b i) (v q)))
	      ((= i q) (setf (b i) (v r)))
	      (t (setf (b i) (v i))))))
    b))

#+nil
(#+nil print-binary-matrix
 #+nil print-matrix
 progn
 (progn (def-problem 1 1 1 '((6 2)))
  (let* ((x (setup))
	 (md (make-disjoint-coefficient-matrix x))
	 (c (make-connection-matrix (length x)))
	 (cmc (make-connected-coefficient-matrix
	       c
	       md)))
    (print-matrix md)
    (print-matrix c)
    (print-matrix cmc)
    (multiple-value-bind (m y) (transpose-rhs cmc)
      (print-matrix m)
      (format t "~A~%" y)
      (format t "~a~%" (solve m y))
      (format t "~A~%" (v-swap-cols y 1 0))
      (format t "~a~%" (solve (swap-rows m 1 0) 
			      (v-swap-cols y 1 0)))))))

(declaim (ftype (function (mat) (values mat vec &optional))
		transpose-rhs))
(defun transpose-rhs (mcon)
  (destructuring-bind (r c) (array-dimensions mcon)
    (assert (= r c))
    (let* ((v (get-param :sending-end-voltage))
	   (p (1- r))
	   (mm (make-mat p))
	   (rhs (make-vec p)))
      (with-arrays (mm mcon rhs)
       (dotimes (j p)
	 (dotimes (i p)
	   (setf (mm i j) (mcon i j))))
       (dotimes (i p)
	 (setf (rhs i) (* -1 (mcon i (1- c)) v))))
      (values mm rhs))))

(declaim (ftype (fun1 () vec) solve-lossy-line))
(defun solve-lossy-line ()
  (let* ((x (setup))
	 (mcon (make-connected-coefficient-matrix
		(make-connection-matrix (length x))
		(make-disjoint-coefficient-matrix x))))
    (multiple-value-bind (m y) (transpose-rhs mcon)
      (let* ((sol (solve m y))
	     (n (length sol))
	     (res (make-array (1+ n)
			      ;; append end voltage
			      :element-type 'single-float)))
	#+nil(with-arrays (res sol)
	 (dotimes (i n)
	   (setf (res i) (sol i)))
	 (setf (res n) (get-param :sending-end-voltage)))
	res))))

#+nil
(def-problem 1s0 1s0 1s0 '((4 2)))
#+nil
(solve-lossy-line)

(declaim (ftype (fun1 (vec single-float) fixnum)
		find-interval))
(defun find-interval (x z)
  (with-arrays (x v)
     (dotimes (i (1- (length x)))
       (when (<= (x i) z (x (1+ i)))
	 (return-from find-interval i)))
     (error "value z=~a not in range of x." z)))

(declaim (ftype (fun1 (vec vec single-float) single-float)
		interpolate-solution))
(defun interpolate-solution (x v z)
  (with-arrays (v x)
   (let* ((i (find-interval x z))
	  (l (x i))
	  (r (x (1+ i)))
	  (len (- r l)))
     (/ (+ (* (v (1+ i)) (- z l))
	   (* (v i) (- r z))) len))))

#+nil
(let ((x (setup))
      (v (solve-lossy-line)))
  (interpolate-solution x v .22))

(defun analytic-solution (z)
  (let* ((r (get-param :resistivity))
	 (g (get-param :conductivity))
	 (v (get-param :sending-end-voltage))
	 (l (let ((s (setup)))
	      (aref s (1- (length s)))))
	 (a (sqrt (* r g))))
    (* v (/ (+ (exp (* a z)) (exp (* -1 a z)))
	    (+ (exp (* a l)) (exp (* -1 a l)))))))

#+nil
(analytic-solution 0)

(defun compare (&optional (n 10))
  (let* ((x (setup))
	 (v (solve-lossy-line))
	 (l (aref x (1- (length x)))))
   (loop for i upto  n collect
	(let ((z (/ (* i l) n)))
	  (list 
	   z 
	   (analytic-solution z)
	   (interpolate-solution x v z))))))

#+nil
(compare)
#+nil ; ecl
(progn
 (format t "~{~{~6f ~}~%~}~%"
	 (compare))
 (ext:quit 0))

#+nil
(defun plot ()
  (with-open-file (s "/dev/shm/o.dat" :direction :output
		     :if-does-not-exist :create
		     :if-exists :supersede)
   (let ((dat (compare (length (setup)))))
     (dolist (e dat)
       (format s "~f ~f~%" (first e) (- (second e)
					  (third e))))
     #+nil (format s "~{~{~6f ~}~%~}~%" dat)))
  (with-open-file (s "/dev/shm/o.gp" :direction :output
		     :if-does-not-exist :create
		     :if-exists :supersede)
    (format s "set term postscript~%")
    (format s "set output \"/dev/shm/o.eps\"~%")
    #+nil (format s "plot \"/dev/shm/o.dat\" u 1:2 w l, \"/dev/shm/o.dat\" u 1:3 w l~%")
    (format s "plot \"/dev/shm/o.dat\" u 1:2 w l~%"))
  (sb-ext:run-program "/usr/bin/gnuplot" (list "/dev/shm/o.gp")))

#+nil
(progn
  (def-problem 1s0 1s0 1s0 '((85 2)))
  (time (plot)))