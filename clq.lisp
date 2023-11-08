(defstruct machine
  quantum-state
  measurement-register)

(defun make-quantum-state (n)
  (let ((s (make-array (expt 2 n) :initial-element 0.0d0)))
    (setf (aref s 0) 1.0d0)
    s))

(defun dimension-qubits (d)
  (- (integer-length d) 1))

(defun apply-operator (unitary-operator ket)
  (let* ((oper-size (array-dimension mat 0))
	 (new-ket (make-array oper-size :initial-element 0.0d0)))
    (dotimes (i oper-size)
      (let ((element 0.0d0))
	(dotimes (j oper-size)
	  (incf element (* (aref unitary-operator i j) (aref ket j))))
	(setf (aref new-ket i) element)))
    (replace ket new-ket)))
