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
  (let* ((oper-size (array-dimension unitary-operator 0))
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

(defun sample (state)
  (let ((u (random 1.0d0)))
    (dotimes (i (length state))
      (decf u (expt (abs (aref state i)) 2))
      (when (minusp u) (return i)))))

(defun collapse (state qubit)
  (fill state 0.0d0)
  (setf (aref state qubit) 1.0d0))

(defparameter +I+ #2A((1 0)
		      (0 1))
	      "The identity gate.")

(defun apply-gate (state unitary-operator qubits)
  (assert (= (length qubits)
	     (dimension-qubits (array-dimension unitary-operator 0))))
  (if (= (length qubits) 1)
      (apply-singleq-gate state unitary-operator (first qubits))
      (apply-multiq-gate state unitary-operator qubits)))

(defun kronecker-product (unitary-A unitary-B)
  (destructuring-bind (m n) (array-dimensions unitary-A)
    (destructuring-bind (p q) (array-dimensions unitary-B)
      (let ((kron-dot-unitary (make-array (list (* m p) (* n q)))))
	(dotimes (i m kron-dot-unitary)
	  (dotimes (j n)
	    (let ((unitary-A-ij (aref unitary-A i j))
		  (pointer-row (* i p))
		  (pointer-col (* j q)))
	      (dotimes (k p)
		(dotimes (l q)
		  (setf (aref kron-dot-unitary
			      (+ pointer-row k) (+ pointer-col l))
			(* unitary-A-ij (aref unitary-B k l))))))))))))

(defun kronecker-expt (unitary n)
  (cond
    ((< n 1) #2A((1)))
    ((= n 1) unitary)
    (t (kronecker-product (kronecker-expt unitary (- n 1)) unitary))))

(defun lift (unitary i n)
  (let ((left-factors (kronecker-expt +I+ (- n i (dimension-qubits
						  (array-dimension unitary 0)))))
	(right-factors (kronecker-expt +I+ i)))
    (kronecker-product left-factors
		       (kronecker-product unitary right-factors))))

(defparameter +SWAP+ #2A((1 0 0 0)
			 (0 0 1 0)
			 (0 1 0 0)
			 (0 0 0 1)))

(defun perm-as-trans (permutation)
  (let ((transpositions nil))
    (dotimes (natural (length permutation) (nreverse transpositions))
      (let ((permuted (elt permutation natural)))
	(loop :while (< permuted natural) :do
	  (setf permuted (elt permutation permuted)))
	(when (> permuted natural)
	  (push (cons natural permuted) transpositions))))))

(defun trans-as-adjacent (transpositions)
  (mapcan (lambda (trans)
            (let ((start (car trans))
                  (end (cdr trans)))
              (if (= end (1+ start))
                  (list start)
                  (nconc (loop :for i :from start :below (1- end) :collect i)
			 (loop :for i :from (1- end) :downto start :collect i)))))
          transpositions))

(defun apply-singleq-gate (state unitary q)
  (apply-operator (lift unitary q (dimension-qubits (length state)))
                  state))

(defun apply-multiq-gate (state unitary qubits)
  (let ((n (dimension-qubits (length state))))
    (labels ((trans-to-op (trans)
               (reduce #'compose-operator trans
		       :key (lambda (i) (lift +swap+ i n)))))
      (let* ((unitary-init (lift unitary 0 n))
             (perm (append (reverse qubits)
			   (remove-if (lambda (i) (member i qubits))
                                      (loop for i :below n :collect i))))
             (trans (trans-as-adjacent (perm-as-trans perm)))
             (to->from (trans-to-op trans))
             (from->to (trans-to-op (reverse trans)))
             (unitary-conform (compose-operator to->from
						(compose-operator unitary-init
								  from->to))))
        (apply-operator unitary-conform state)))))

(defun clq (qprog machine)
  (loop :for (instruction . parameters) :in qprog
        :do (ecase instruction
              ((GATE)
               (destructuring-bind (gate &rest qubits) parameters
                 (apply-gate (machine-quantum-state machine) gate qubits)))
              ((MEASURE)
               (observe machine)))
        :finally (return machine)))

(defparameter +H+
  (make-array '(2 2)
	      :initial-contents
	      (let ((s (/ (sqrt 2))))
                (list (list s s)
		      (list s (- s))))))

(defparameter +CNOT+ #2A((1 0 0 0)
                         (0 1 0 0)
                         (0 0 0 1)
                         (0 0 1 0)))

(defun bell-state (p q)
  `((GATE ,+H+ ,p)
    (GATE ,+CNOT+ ,p ,q)))

(clq (bell-state 0 1)
     (make-machine :quantum-state (make-quantum-state 2)
                   :measurement-register 0))

(defun cphase (angle)
  (make-array '(4 4) :initial-contents `((1 0 0 0)
                                         (0 1 0 0)
                                         (0 0 1 0)
                                         (0 0 0 ,(cis angle)))))

(defun cphase-reorder (angle)
  (make-array '(4 4) :initial-contents `((1 0 0 0)
                                         (0 ,(cis angle) 0 0)
                                         (0 0 1 0)
                                         (0 0 0 1))))

(defun qft (qubits)
  (labels ((bit-reversal (qubits)
             (let ((n (length qubits)))
               (if (< n 2)
                   nil
                   (loop :repeat (floor n 2)
                         :for qs :in qubits
                         :for qe :in (reverse qubits)
                         :collect `(GATE ,+swap+ ,qs ,qe)))))
           (%qft (qubits)
             (destructuring-bind (q . qs) qubits
               (if (null qs)
                   (list `(GATE ,+H+ ,q))
                   (let ((cR (loop :with n := (1+ (length qs))
                                   :for i :from 1
                                   :for qi :in qs
                                   :for angle := (/ pi (expt 2 (- n i)))
                                   :collect `(GATE ,(cphase angle) ,q ,qi))))
                     (append
                      (qft qs)
                      cR
                      (list `(GATE ,+H+ ,q))))))))
    (append (%qft qubits) (bit-reversal qubits))))

(qft '(0 1 2))

(clq (qft '(0 1 2))
     (make-machine :quantum-state (make-quantum-state 3)
                   :measurement-register 0))

(defun shor-fifteen ()
  `((GATE ,+H+ 0)
    (GATE ,+H+ 1)
    (GATE ,+H+ 2)
    (GATE ,+CNOT+ 2 3)
    (GATE ,+CNOT+ 2 4)
    (GATE ,+H+ 1)
    (GATE ,(cphase-reorder (/ pi 2)) 0 1)
    (GATE ,+H+ 0)
    (GATE ,(cphase (/ pi 4)) 1 2)
    (GATE ,(cphase (/ pi 2)) 0 2)
    (GATE ,+H+ 2)
    (MEASURE)))

(clq (shor-fifteen)
     (make-machine :quantum-state (make-quantum-state 5)
                   :measurement-register 0))

(defun counts (realizations qubits)
  (let* ((hilbert (expt 2 qubits))
         (results (make-array hilbert :initial-element 0)))
    (dotimes (i realizations results)
      (let ((state (machine-quantum-state
                    (clq (shor-fifteen)
                         (make-machine :quantum-state (make-quantum-state qubits)
                                       :measurement-register 0)))))
	(dotimes (j hilbert)
	  (when (> (aref state j) 0) (incf (aref results j))))))))
(counts 1024 5)
