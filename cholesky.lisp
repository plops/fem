(defparameter *s*
  '((4 4 4 0 0 0 0 0)
    (4 5 2 0 0 0 0 0)
    (4 2 9 1 2 1 0 0)
    (0 0 1 10 -7 1 0 0)
    (0 0 2 -7 14 0 2 1)
    (0 0 1 1 0 9 -8 -2)
    (0 0 0 0 2 -8 9 0)
    (0 0 0 0 1 -2 0 6)))

(defun afloat1 (list)
  (make-array (length list)
	      :element-type 'single-float
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


(defun cholesky (s)
  (declare ((simple-array single-float 2) s)
	   (values (simple-array single-float 2) &optional))
  (destructuring-bind (r c) (array-dimensions s)
    (assert (= r c))
    (let ((l (make-array (list r r)
			 :element-type 'single-float)))
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
	 (y (make-array n
			:element-type 'single-float)))
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
  (declare ((simple-array single-float 2) a)
	   (values (simple-array single-float 2) &optional))
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
	 (y (make-array n
			:element-type 'single-float)))
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