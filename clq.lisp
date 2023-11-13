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

(defun compose-operator (unitary-left unitary-right)
  (let* ((m (array-dimension unitary-left 0))
	 (n (array-dimension unitary-left 1))
	 (p (array-dimension unitary-right 1))
	 (unitary-new (make-array (list m p) :initial-element 0.0d0)))
    (dotimes (i m unitary-new)
      (dotimes (j p)
	(let ((dot-prod 0.0d0))
	  (dotimes (k n)
	    (setf dot-prod (+ dot-prod
			      (* (aref unitary-left i k)
				 (aref unitary-right k j)))))
	  (setf (aref unitary-new i j) dot-prod))))))

(defun observe (machine)
  (let ((b (sample (machine-quantum-state machine))))
    (collapse (machine-quantum-state machine) b)
    (setf (machine-measurement-register machine) b)
    machine))
