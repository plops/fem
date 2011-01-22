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

(defmacro make-mat (n &rest rest)
  `(make-array (list ,n ,n) :element-type 'single-float
	       ,@rest))

(defun afloat1 (list)
  (make-vec (length list)
	    :initial-contents
	    (mapcar #'(lambda (x) (* 1s0 x)) list)))

(defun afloat2 (lists)
  (make-mat (length lists)
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


(defun cholesky (s)
  (declare (mat s)
	   (values mat &optional))
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

(defun transpose (a)
  (declare (mat a)
	   (values mat &optional))
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

(defun solve (m y)
  "Solve M x = y for x."
  (declare (mat m)
	   (vec y)
	   (values vec &optional))
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
    (:segments . ((12 1s0)
		  (12 .5s0)
		  (12 .1s0))))) ; another element with length .5
#+nil
(assoc :segments *problem*)

(defun get-param (name)
  (rest (assoc name *problem*)))
#+nil
(get-param :segments)

(defun setup ()
   (let* ((segs (get-param :segments))	; elements with length
	  (x (destructuring-bind (num len) (first segs)
	       (- (/ len num))))
	  (pos (let ((res ()))
		 (mapcar #'(lambda (nl)
			     (destructuring-bind (num len) nl
			       (dotimes (i num)
				 (push (incf x (/ (* len (1+ i)) num))
				       res))))
			 segs)
		 (reverse res))))
     (make-array (length pos) :element-type 'single-float
		 :initial-contents pos)))
#+nil 
(setup)

(defun ndis (ncon)
  (* 2 (1- ncon)))

(defun make-connection-matrix (ncon)
  (let* ((ndis (ndis ncon))
	 (c (make-array (list ncon ndis)
			:element-type 'single-float
			:initial-element 0s0)))
    (with-arrays (c)
      (dotimes (i ncon)
	(setf (c i i) 1s0))
      (dotimes (i (- ndis ncon))
	(setf (c (1+ i) (+ ncon i)) 1s0))
      c)))
#+nil
(make-connection-matrix 6)

(defun make-disjoint-coefficient-matrix (x)
  (declare (vec x)
	   (values mat &optional))
  (let* ((ncon (length x))
	 (ndis (ndis ncon))
	 (c (make-mat ndis))
	 (co (get-param :conductivity))
	 (re (get-param :resistivity)))
    (with-arrays (x c)
     (loop for i below (- ndis 1) by 2 do
	(let* ((j (1+ i))
	       (p (floor j 2))
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
(make-disjoint-coefficient-matrix (setup))



;; their code
;; M_kl = c_ki Md_ij c_lj
;;   oo     o*    **   o*
;; * .. many, o .. not many


(defun make-connected-coefficient-matrix (c m)
  (declare (mat c m)
	   (values mat &optional))
  (destructuring-bind (mr mcol) (array-dimensions m)
   (destructuring-bind (r col) (array-dimensions c)
     ;(assert (= mr mcol col))
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
#+nil
(let ((x (setup)))
   (make-connected-coefficient-matrix
    (make-connection-matrix (length x))
    (make-disjoint-coefficient-matrix x)))


(defun transpose (a)
  (declare (mat a)
	   (values mat &optional))
  (let ((b (make-array (reverse (array-dimensions a))
		       :element-type 'single-float)))
    (destructuring-bind (r c) (array-dimensions a)
      (with-arrays (a b)
       (dotimes (i r)
	 (dotimes (j c)
	   (setf (b j i) (a i j)))))
      b)))

(defun transpose-rhs (mcon)
  (declare (mat mcon)
	   (values mat vec &option))
  (destructuring-bind (r c) (array-dimensions mcon)
    (assert (= r c))
    (let* ((v (get-param :sending-end-voltage))
	   (mm (make-array (list (1- r) (1- c))
			   :element-type 'single-float))
	   (rhs (make-vec (1- r))))
      (with-arrays (mm mcon rhs)
       (dotimes (j (1- c))
	 (dotimes (i (1- r))
	   (setf (mm i j) (mcon i j))))
       (dotimes (i (1- r))
	 (setf (rhs i) (* -1 (mcon i (1- c)) v))))
      (values mm rhs))))

(defun solve-lossy-line ()
  (let* ((x (setup))
	 (mcon  (make-connected-coefficient-matrix
		 (make-connection-matrix (length x))
		 (make-disjoint-coefficient-matrix x))))
    (multiple-value-bind (m y) (transpose-rhs mcon)
      (let* ((sol (solve m y))
	    (res (make-array (1+ (length sol)) ;; append end voltage
			     :element-type 'single-float)))
	(with-arrays (res sol)
	 (dotimes (i (length sol))
	   (setf (res i) (sol i)))
	 (setf (res (length sol)) (get-param :sending-end-voltage)))
	res))))

#+nil
(solve-lossy-line)

(defun find-interval (x z)
  (declare (vec x)
	   (single-float z)
	   (values fixnum &optional))
  (with-arrays (x v)
     (dotimes (i (1- (length x)))
       (when (<= (x i) z (x (1+ i)))
	 (return-from find-interval i)))
     (error "value z=~a not in range of x." z)))

(defun interpolate-solution (x v z)
  (declare (vec x v)
	   (single-float z)
	   (values single-float &optional))
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
  (interpolate-solution x v .4))

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
(analytic-solution .3)

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

(defun plot ()
  (with-open-file (s "/dev/shm/o.dat" :direction :output
		     :if-does-not-exist :create
		     :if-exists :supersede)
   (let ((dat (compare 100)))
     (format s "~{~{~6f ~}~%~}~%" dat)))
  (with-open-file (s "/dev/shm/o.gp" :direction :output
		     :if-does-not-exist :create
		     :if-exists :supersede)
    (format s "set term postscript~%")
    (format s "set output \"/dev/shm/o.eps\"~%")
    (format s "plot \"/dev/shm/o.dat\" u 1:2 w l, \"/dev/shm/o.dat\" u 1:3 w l~%"))
  (sb-ext:run-program "/usr/bin/gnuplot" (list "/dev/shm/o.gp")))

#+nil
(plot)