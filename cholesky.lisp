(defparameter *s*
  '((4 4 4 0 0 0 0 0)
    (4 5 2 0 0 0 0 0)
    (4 2 9 1 2 1 0 0)
    (0 0 1 10 -7 1 0 0)
    (0 0 2 -7 14 0 2 1)
    (0 0 1 1 0 9 -8 -2)
    (0 0 0 0 2 -8 9 0)
    (0 0 0 0 1 -2 0 6)))

(defparameter s 
  (make-array 
   (list 8 8)
   :initial-contents
   *s*))

(defmacro with-arrays (arrays &body body)
  `(macrolet ,(mapcar 
	       (lambda (array)
		 `(,array (&rest indices) 
			  `(aref ,',array ,@indices)))
	       arrays)
     ,@body))


(defun cholesky (s)
  (declare ((simple-array single-float 2) a l)
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



(cholesky (make-array '(8 8) 
		      :element-type 'single-float 
		      :initial-contents (mapcar #'(lambda (row) 
						    (mapcar #'(lambda (x) (* 1s0 x)) row))
						*s*)))

; fem for engineers p 475
(defun cholesky2 (a l)
  (declare ((simple-array single-float 2) a l)
	   (values (simple-array single-float 2) &optional))
  (let ((n (array-dimension a 0)))
    (with-arrays (a l)
      (dotimes (j n)
	(setf (l j j) (a j j))
	(dotimes (k (1- j))
	  (decf (l j j) (expt (l j k) 2)))
	(when (< 0 (l j j))
	  (setf (l j j) (sqrt (l j j)))
	  (loop for i from j below n do
	       (setf (l i j) (a i j))
	       (dotimes (k (1- j))
		 (decf (l i j) (* (l i k)
				  (l j k))))
	       (setf (l i j) (/ (l i j)
			   (l j j))
		(l j i) (l i j)))))))
  l)

#+nil
(cholesky2 s s)

; numerical recipes in C p 121
(defun chold (a)
  (let* ((n (array-dimension a 0))
	 (l (make-array (list n n)))
	 (s 0s0))
    (with-arrays (a l)
     (dotimes (i n)
       (dotimes (j n)
	 (setf s (a i j))
	 (dotimes (k (1- j))
	   (decf s (* (a i k)
		      (a j k))))
	 (if (= i j)
	     (setf (l i i) (sqrt s))
	     (setf (l j i) (/ s (l i i)))))))
    l))

(chold (make-array '(8 8) :initial-contents *s*))