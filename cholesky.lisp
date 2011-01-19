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
    (:segments . ((2 2s0) ; 2 elements each with length 1 near sending end
		 (1 .5s0))))) ; another element with length .5
#+nil
(assoc :segments *problem*)
(defun setup (prob)
   (let* ((segs (rest (assoc :segments prob))) ; elements with length
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
     pos))
 
(setup *problem*)

#+nil
(defun make-connection-matrix ()
  (let* ((ncon *nodes*)
	 (ndis (* 2 (1- ncon)))
	 (c (make-array (list ndis ncon))))
    (with-arrays (c)
      (dotimes (i ncon)
	(dotimes (j ndis)
	  (setf (c j i)
		(if (= i (1+ (floor j 2)))
		    1s0 0s0))))
      c)))
#+nil
(make-connection-matrix)