(module image-processing *
(import chicken scheme extras lolevel)
(use traversal define-structure scheme2c-compatibility
     linear-algebra format miscmacros imlib2 files srfi-1)

;;; Accessing colors

(define r x)
(define g y)
(define b z)

(define h x)
(define s y)
(define v z)

;;; Color space conversions

(define *max-red* 255)
(define *max-green* 255)
(define *max-blue* 255)
(define *max-grey* 255)
(define *max-hue* 360)
(define *max-saturation* 100)
(define *max-value* 100)

(define (rgb->hsv rgb)
 ;; The constants are hardwired to be inexact for efficiency.
 (let* ((red (/ (vector-ref rgb 0) *max-red*))
	(green (/ (vector-ref rgb 1) *max-green*))
	(blue (/ (vector-ref rgb 2) *max-blue*))
	(value (max red green blue))
	(m (min red green blue))
	(saturation (if (zero? value) 0.0 (/ (- value m) value)))
	(hue (if (zero? saturation)
		 0.0
		 (/ (let ((rc (/ (- value red) (- value m)))
			  (gc (/ (- value green) (- value m)))
			  (bc (/ (- value blue) (- value m))))
		     (cond ((= red value) (- bc gc))
			   ((= green value) (+ 2.0 (- rc bc)))
			   (else (+ 4.0 (- gc rc)))))
		    6.0)))
	(hue (if (negative? hue) (+ hue 1.0) hue)))
  (vector (inexact->exact (floor (* *max-hue* hue)))
	  (inexact->exact (floor (* *max-saturation* saturation)))
	  (inexact->exact (floor (* *max-value* value))))))
(define (hsv->rgb hsv)
 ;; The constants are hardwired to be inexact for efficiency.
 (let ((hue (/ (vector-ref hsv 0) *max-hue*))
       (saturation (/ (vector-ref hsv 1) *max-saturation*))
       (value (* (max *max-red* *max-green* *max-blue*)
		 (/ (vector-ref hsv 2) *max-value*))))
  (if (zero? saturation)
      (vector (inexact->exact (floor value))
	      (inexact->exact (floor value))
	      (inexact->exact (floor value)))
      (let* ((hue (if (negative? hue) (+ hue 1.0) hue))
	     (hue (* 6.0 hue))
	     (fract (- hue (floor hue)))
	     (new1 (inexact->exact (floor (* value (- 1.0 saturation)))))
	     (new2 (inexact->exact
		    (floor (* value (- 1.0 (* saturation fract))))))
	     (new3 (inexact->exact
		    (floor (* value (- 1.0 (* saturation (- 1.0 fract)))))))
	     (value (inexact->exact (floor value))))
       (case (inexact->exact (floor hue))
	((0) (vector value new3 new1))
	((1) (vector new2 value new1))
	((2) (vector new1 value new3))
	((3) (vector new1 new2 value))
	((4) (vector new3 new1 value))
	((5) (vector value new1 new2))
	((6) (vector value new3 new1))
	(else (error "Inappropriate hue angle")))))))

(define (rgb->cd rgb)
 (let* ((red (vector-ref rgb 0))
	(green (vector-ref rgb 1))
	(blue (vector-ref rgb 2))
	(intensity (max (+ red green blue) 1)))
  (vector (inexact->exact (floor (* *max-red* (/ red intensity))))
	  (inexact->exact (floor (* *max-green* (/ green intensity)))))))

(define (rgb->cmyk rgb)
 (let* ((r (vector-ref rgb 0))
	(g (vector-ref rgb 1))
	(b (vector-ref rgb 2))
	(k (max r g b)))
  (if (zero? k)
      '#(0 0 0 0)
      `#(,(exact-round (* 255 (- 1 (/ r k))))
	 ,(exact-round (* 255 (- 1 (/ g k))))
	 ,(exact-round (* 255 (- 1 (/ b k))))
	 ,k))))

(define (rgb->uv-hsv c)
 (let ((hsv (rgb->hsv c)))
  `#(,(exact-round (+ (- (* 0.169 (x c))) (* 0.331 (y c)) (* 0.5 (z c))))
     ,(exact-round (+ (* 0.5 (x c)) (- (* 0.418 (y c))) (- (* 0.082 (z c)))))
     ,(x hsv) ,(y hsv) ,(z hsv))))

(define (rgb->xyz c)
 (define (adjust c)
  (let ((c (/ c 255)))
   (* 100 (if (> c 0.04045)
	      (expt (/ (+ c 0.055) 1.055) 2.4)
	      (/ c 12.92)))))
 (let ((r (adjust (r c))) (g (adjust (g c))) (b (adjust (b c))))
  ;; Observer. = 2 degrees, Illuminant = D65
  `#(,(+ (* r 0.4124) (* g 0.3576) (* b 0.1805))
     ,(+ (* r 0.2126) (* g 0.7152) (* b 0.0722))
     ,(+ (* r 0.0193) (* g 0.1192) (* b 0.9505)))))
(define (xyz->rgb c)
 (let ((c (map-vector (lambda (c) (/ c 100)) c)))
  (let ((r (+ (* (x c)  3.2406) (* (y c) -1.5372) (* (z c) -0.4986)))
	(g (+ (* (x c) -0.9689) (* (y c)  1.8758) (* (z c)  0.0415)))
	(b (+ (* (x c)  0.0557) (* (y c) -0.2040) (* (z c)  1.0570))))
   (define (adjust c)
    (exact-round
     (* 255 (if (> c 0.0031308)
		(+ (* 1.055 (expt c (/ 2.4))) -0.055)
		(* 12.92 c)))))
   ;; Observer. = 2 degrees, Illuminant = D65
   (map-vector adjust `#(,r ,g ,b)))))
(define (xyz->l*ab c)
 (define (adjust c)
  (if (> c 0.008856)
      (expt c (/ 3))
      (+ (* 7.787 c) (/ 16 116))))
 (let ((x (adjust (/ (x c)  95.047)))
       (y (adjust (/ (y c) 100.000)))
       (z (adjust (/ (z c) 108.883))))
  ;; Observer. = 2 degrees, Illuminant = D65
  `#(,(- (* 116 y) 16)
     ,(* 500 (- x y))
     ,(* 200 (- y z)))))
(define (l*ab->xyz c)
 (define (adjust c)
  (if (> (expt c 3) 0.008856)
      (expt c 3)
      (/ (- c (/ 16 116)) 7.787)))
 (let* ((cy (/ (+ 16 (x c)) 116))
	(cx (+ (/ (y c) 500) cy))
	(cz (- cy (/ (z c) 200))))
  ;; Observer. = 2 degrees, Illuminant = D65
  `#(,(* (adjust cx)  95.047)
     ,(* (adjust cy) 100.000)
     ,(* (adjust cz) 108.883))))

(define (rgb->l*ab c) (xyz->l*ab (rgb->xyz c)))
(define (l*ab->rgb c) (xyz->rgb (l*ab->xyz c)))

(define (rgb->html c)
 (define (number-radix->padded-string-of-length number radix length)
  (when (negative? number) (fuck-up))
  (let ((string (number->string number radix)))
   (string-append (make-string (- length (string-length string)) #\0) string)))
 (qmap-reduce
  string-append
  "#"
  (lambda (a) (number-radix->padded-string-of-length a 16 2))
  (vector->list c)))

;;; Images

(define-structure pbm raw? bitmap)
(define-structure pgm raw? maxval grey)
(define-structure ppm raw? maxval red green blue)

(define (pnm? m) (or (pbm? m) (pgm? m) (ppm? m)))

(define (image-ref i p)
 (cond ((pbm? i) (matrix-ref (pbm-bitmap i) (y p) (x p)))
       ((pgm? i) (matrix-ref (pgm-grey i) (y p) (x p)))
       ((ppm? i)
	`#(,(matrix-ref (ppm-red i) (y p) (x p))
	   ,(matrix-ref (ppm-blue i) (y p) (x p))
	   ,(matrix-ref (ppm-green i) (y p) (x p))))
       ((matrix? i) (matrix-ref i (y p) (x p)))
       (else (fuck-up))))

(define (pbm-ascii pbm) (make-pbm #f (pbm-bitmap pbm)))

(define (pnm-copy pnm)
 (cond
  ((ppm? pnm)
   (make-ppm (ppm-raw? pnm)
	     (ppm-maxval pnm)
	     (matrix-copy (ppm-red pnm))
	     (matrix-copy (ppm-green pnm))
	     (matrix-copy (ppm-blue pnm))))
  ((pgm? pnm)
   (make-pgm (pgm-raw? pnm) (pgm-maxval pnm) (matrix-copy (pgm-grey pnm))))
  ((pbm? pnm) (make-pbm (pbm-raw? pnm) (matrix-copy (pbm-bitmap pnm))))
  (else (error "Argument is not a PNM"))))

;;; Creating images

(define (read-pnm pathname)
 (define (read-pnm port)
  (define (read-pbm raw?)
   (let* ((width (read port))
	  (height (read port))
	  (bitmap (make-matrix height width 0)))
    (call-with-current-continuation
     (lambda (return)
      (cond
       (raw? (error "Cannot (yet) read a raw pbm image"))
       (else
	(for-each-n (lambda (y)
		     (for-each-n (lambda (x)
				  (let ((v (read port)))
				   (when (eof-object? v) (return #f))
				   ;; Yes, it really is the case (at least
				   ;; according to xv) that 0 means white and
				   ;; 1 means black for ascii pbm images.
				   (matrix-set! bitmap y x (zero? v))))
                      width))
         height)))))
    (make-pbm raw? bitmap)))
  (define (read-pgm raw?)
   (let* ((width (read port))
	  (height (read port))
	  (maxval (read port))
	  (size (* width height))
	  (grey (make-matrix height width 0)))
    (call-with-current-continuation
     (lambda (return)
      (cond
       (raw?
	(read-char port)
	(for-each-n
          (lambda (y)
           (for-each-n (lambda (x)
                        (let ((c (read-char port)))
                         (when (eof-object? c) (return #f))
                         (let ((int (if (< 255 maxval)
                                        (+ (bit-lsh (char->integer c) 8)
                                           (char->integer (read-char port)))
                                        (char->integer c))))
                          (matrix-set! grey y x int))))
            width))
	 height))
       (else (for-each-n (lambda (y)
			  (for-each-n (lambda (x)
				       (let ((v (read port)))
					(when (eof-object? v) (return #f))
					(matrix-set! grey y x v)))
                           width))
              height)))))
    (make-pgm raw? maxval grey)))
  (define (read-ppm raw?)
   (let* ((width (read port))
	  (height (read port))
	  (maxval (read port))
	  (size (* width height))
	  (red (make-matrix height width 0))
	  (green (make-matrix height width 0))
	  (blue (make-matrix height width 0)))
    (call-with-current-continuation
     (lambda (return)
      (cond (raw? (read-char port)
		  (for-each-n
                    (lambda (y)
                     (for-each-n
                       (lambda (x)
                        (let* ((c1 (read-char port))
                               (c2 (read-char port))
                               (c3 (read-char port)))
                         (when (eof-object? c1) (return #f))
                         (matrix-set! red y x (char->integer c1))
                         (matrix-set! green y x (char->integer c2))
                         (matrix-set! blue y x (char->integer c3))))
                      width))
		   height))
	    (else (for-each-n
                    (lambda (y)
                     (for-each-n (lambda (x)
                                  (let* ((v1 (read port))
                                         (v2 (read port))
                                         (v3 (read port)))
                                   (when (eof-object? v1) (return #f))
                                   (matrix-set! red y x v1)
                                   (matrix-set! green y x v2)
                                   (matrix-set! blue y x v3)))
                      width))
		   height)))))
    (make-ppm raw? maxval red green blue)))
  (let ((format (read port)))
   (read-char port)
   (while (char=? (peek-char port) #\#) (read-line port))
   (case format
    ((p1) (read-pbm #f))
    ((p2) (read-pgm #f))
    ((p3) (read-ppm #f))
    ((p4) (read-pbm #t))
    ((p5) (read-pgm #t))
    ((p6) (read-ppm #t))
    (else (error "Incorrect format for a pnm image")))))
 (cond ((port? pathname) (read-pnm pathname))
       ((string=? pathname "-") (read-pnm (current-input-port)))
       (else (call-with-input-file pathname read-pnm))))

(define (write-pnm pnm pathname)
 (define (write-pnm port)
  (define (write-pbm pbm)
   (let ((width (pnm-width pbm))
	 (height (pnm-height pbm))
	 (bitmap (pbm-bitmap pbm)))
    (write (if (pbm-raw? pbm) 'p4 'p1) port)
    (newline port)
    (write width port)
    (write-char #\space port)
    (write height port)
    (newline port)
    (if (pbm-raw? pbm)
	(error "Cannot (yet) write a raw pbm image")
	(for-each-n (lambda (y)
		     (for-each-n (lambda (x)
				  ;; Yes, it really is the case (at least
				  ;; according to xv) that 0 means white and
				  ;; 1 means black for ascii pbm images.
				  (write (if (matrix-ref bitmap y x) 0 1) port)
				  (newline port))
                      width))
         height))))
  (define (write-pgm pgm)
   (let ((width (pnm-width pgm))
	 (height (pnm-height pgm))
	 (grey (pgm-grey pgm)))
    (when (pgm-raw? pgm)
     (for-each-n
       (lambda (y)
        (for-each-n (lambda (x)
                     (when (> (matrix-ref grey y x) 255)
                      (error "Grey value too large for raw pgm file format")))
         width))
      height))
    (write (if (pgm-raw? pgm) 'p5 'p2) port)
    (newline port)
    (write width port)
    (write-char #\space port)
    (write height port)
    (newline port)
    (write (pgm-maxval pgm) port)
    (newline port)
    (if (pgm-raw? pgm)
	(for-each-n
          (lambda (y)
           (for-each-n
             (lambda (x) (write-char (integer->char (matrix-ref grey y x)) port))
            width))
	 height)
	(for-each-n (lambda (y)
		     (for-each-n (lambda (x)
				  (write (matrix-ref grey y x) port)
				  (newline port))
                      width))
         height))))
  (define (write-ppm ppm)
   (let ((width (pnm-width ppm))
	 (height (pnm-height ppm))
	 (red (ppm-red ppm))
	 (green (ppm-green ppm))
	 (blue (ppm-blue ppm)))
    (when (ppm-raw? ppm)
     (for-each-n
       (lambda (y)
        (for-each-n (lambda (x)
                     (when (or (> (matrix-ref red y x) 255)
			      (> (matrix-ref green y x) 255)
			      (> (matrix-ref blue y x) 255))
                      (error "Color value too large for raw ppm file format")))
         width))
      height))
    (write (if (ppm-raw? ppm) 'p6 'p3) port)
    (newline port)
    (write width port)
    (write-char #\space port)
    (write height port)
    (newline port)
    (write (ppm-maxval ppm) port)
    (newline port)
    (if (ppm-raw? ppm)
	(for-each-n
          (lambda (y)
           (for-each-n (lambda (x)
                        (write-char (integer->char (matrix-ref red y x)) port)
                        (write-char (integer->char (matrix-ref green y x)) port)
                        (write-char (integer->char (matrix-ref blue y x)) port))
            width))
	 height)
	(for-each-n (lambda (y)
		     (for-each-n (lambda (x)
				  (write (matrix-ref red y x) port)
				  (newline port)
				  (write (matrix-ref green y x) port)
				  (newline port)
				  (write (matrix-ref blue y x) port)
				  (newline port))
                      width))
         height))))
  (cond ((pbm? pnm) (write-pbm pnm))
	((pgm? pnm) (write-pgm pnm))
	((ppm? pnm) (write-ppm pnm))
	(else (error "Non-PNM argument to WRITE-PNM"))))
 (cond ((port? pathname) (write-pnm pathname))
       ((string=? pathname "-") (write-pnm (current-output-port)))
       (else (call-with-output-file (default-extension
                                     pathname
                                     (cond ((pbm? pnm) "pbm")
                                           ((pgm? pnm) "pgm")
                                           ((ppm? pnm) "ppm")
                                           (else (fuck-up))))
              write-pnm))))

(define (pbm-constant width height bit)
 (make-pbm #t (make-matrix height width bit)))

(define (pgm-constant width height grey)
 (make-pgm #t *max-grey* (make-matrix height width grey)))

(define (ppm-constant width height red green blue)
 (make-ppm #t
	   (max *max-red* *max-green* *max-blue*)
	   (make-matrix height width red)
	   (make-matrix height width green)
	   (make-matrix height width blue)))

(define (pbm-left-vertical-stripe width height left)
 ;; Creates a black (#F) stripe on a white (#T) background.
 (let ((m (make-matrix height width #t)))
  (do ((y 0 (+ y 1))) ((= y height))
   (do ((x 0 (+ x 1))) ((= x left))
    (matrix-set! m y x #f)))
  (make-pbm #t m)))

(define (crop-image pnm x y width height)
 (cond ((pbm? pnm) (make-pbm #f (crop (pbm-bitmap pnm) x y width height)))
       ((pgm? pnm) (make-pgm (pgm-raw? pnm)
			     (pgm-maxval pnm)
			     (crop (pgm-grey pnm) x y width height)))
       ((ppm? pnm) (make-ppm (ppm-raw? pnm)
			     (ppm-maxval pnm)
			     (crop (ppm-red pnm) x y width height)
			     (crop (ppm-green pnm) x y width height)
			     (crop (ppm-blue pnm) x y width height)))
       (else (panic "Image must be one of PBM, PGM or PPM"))))

(define (pnm->ppm pnm)
 (cond ((ppm? pnm) pnm)
       ((pgm? pnm) (pgm->ppm pnm))
       ((pbm? pnm) (pbm->ppm pnm))
       (else (fuck-up))))

(define (pnm->pgm pnm)
 (cond ((ppm? pnm) (ppm->pgm pnm))
       ((pgm? pnm) pnm)
       ((pbm? pnm) (pbm->pgm pnm))
       (else (fuck-up))))

(define (pbm->pgm pbm)
 (unless (pbm? pbm) (error "Argument to PBM->PGM is not a PBM"))
 (make-pgm
  (pbm-raw? pbm)
  *max-grey*
  (map-vector
   (lambda (row) (map-vector (lambda (bit) (if bit *max-grey* 0)) row))
   (pbm-bitmap pbm))))

(define (pgm->ppm pgm)
 (unless (pgm? pgm) (error "Argument to PGM->PPM is not a PGM"))
 (make-ppm (pgm-raw? pgm)
	   (pgm-maxval pgm)
	   (pgm-grey pgm)
	   (pgm-grey pgm)
	   (pgm-grey pgm)))

(define (pbm->ppm pbm) (pgm->ppm (pbm->pgm pbm)))

(define (ppm->pgm ppm)
 (unless (ppm? ppm) (error "Argument to PPM->PGM is not a PPM"))
 (make-pgm
  (ppm-raw? ppm)
  (ppm-maxval ppm)
  (map-vector
   (lambda (red-row green-row blue-row)
    (map-vector (lambda (red green blue)
		 (inexact->exact
		  (floor (+ (* 0.299 red) (* 0.587 green) (* 0.114 blue)))))
		red-row green-row blue-row))
   (ppm-red ppm) (ppm-green ppm) (ppm-blue ppm))))

(define (pgm->pbm pgm threshold)
 (unless (pgm? pgm) (error "Argument to PGM->PBM is not a PGM"))
 (make-pbm (pgm-raw? pgm)
	   (map-vector
	    (lambda (row) (map-vector (lambda (grey) (>= grey threshold)) row))
	    (pgm-grey pgm))))

(define (ppm->pbm ppm threshold) (pgm->pbm (ppm->pgm ppm) threshold))

;;; Basic image information

(define (pnm-width pnm)
 (matrix-columns (cond ((pbm? pnm) (pbm-bitmap pnm))
		       ((pgm? pnm) (pgm-grey pnm))
		       ((ppm? pnm) (ppm-red pnm))
		       (else (error "Argument not PNM")))))

(define (pnm-height pnm)
 (matrix-rows (cond ((pbm? pnm) (pbm-bitmap pnm))
		    ((pgm? pnm) (pgm-grey pnm))
		    ((ppm? pnm) (ppm-red pnm))
		    (else (error "Argument not PNM")))))

(define (ppm-hue ppm)
 (unless (ppm? ppm) (error "Argument to PPM-HUE is not a PPM"))
 (make-pgm (ppm-raw? ppm)
	   *max-hue*
	   (map-vector
	    (lambda (red-row green-row blue-row)
	     (map-vector
	      (lambda (red green blue)
	       (vector-ref (rgb->hsv (vector red green blue)) 0))
	      red-row
	      green-row
	      blue-row))
	    (ppm-red ppm)
	    (ppm-green ppm)
	    (ppm-blue ppm))))

(define (ppm-saturation ppm)
 (unless (ppm? ppm) (error "Argument to PPM-SATURATION is not a PPM"))
 (make-pgm (ppm-raw? ppm)
	   *max-saturation*
	   (map-vector
	    (lambda (red-row green-row blue-row)
	     (map-vector
	      (lambda (red green blue)
	       (vector-ref (rgb->hsv (vector red green blue)) 1))
	      red-row
	      green-row
	      blue-row))
	    (ppm-red ppm)
	    (ppm-green ppm)
	    (ppm-blue ppm))))

(define (ppm-value ppm)
 (unless (ppm? ppm) (error "Argument to PPM-VALUE is not a PPM"))
 (make-pgm (ppm-raw? ppm)
	   *max-value*
	   (map-vector
	    (lambda (red-row green-row blue-row)
	     (map-vector
	      (lambda (red green blue)
	       (vector-ref (rgb->hsv (vector red green blue)) 2))
	      red-row
	      green-row
	      blue-row))
	    (ppm-red ppm)
	    (ppm-green ppm)
	    (ppm-blue ppm))))

(define (ppm-mean image #!optional (colour-transform identity))
 (let ((acc (colour-transform '#(0 0 0))))
  (for-each-vector
   (lambda (red-row green-row blue-row)
    (for-each-vector
     (lambda (r g b)
      (set! acc (v+ acc (colour-transform `#(,r ,g ,b)))))
     red-row green-row blue-row))
   (ppm-red image) (ppm-green image) (ppm-blue image))
  (k*v (/ 1 (* (pnm-height image) (pnm-width image))) acc)))

(define (ppm-covariance image #!optional (colour-transform identity))
 (let* ((mu (ppm-mean image colour-transform))
	(acc (make-matrix (vector-length mu) (vector-length mu) 0)))
  (for-each-vector
   (lambda (red-row green-row blue-row)
    (for-each-vector
     (lambda (r g b)
      (set!
       acc
       (m+ acc (self-outer-product * (v- (colour-transform `#(,r ,g ,b)) mu)))))
     red-row green-row blue-row))
   (ppm-red image) (ppm-green image) (ppm-blue image))
  (k*m (/ 1 (* (pnm-height image) (pnm-width image))) acc)))

;;; Logical and morphological operations

(define (pbm-and pbm1 pbm2)
 (unless (and (pbm? pbm1)
            (pbm? pbm2)
            (= (pnm-width pbm1) (pnm-width pbm2))
            (= (pnm-height pbm1) (pnm-height pbm2)))
  (error "Arguments to PBM-AND are not matching PBMs"))
 (make-pbm (pbm-raw? pbm1)
	   (map-vector
	    (lambda (row1 row2)
	     (map-vector (lambda (bit1 bit2) (and bit1 bit2)) row1 row2))
	    (pbm-bitmap pbm1)
	    (pbm-bitmap pbm2))))

(define (pbm-or pbm1 pbm2)
 (unless (and (pbm? pbm1)
            (pbm? pbm2)
            (= (pnm-width pbm1) (pnm-width pbm2))
            (= (pnm-height pbm1) (pnm-height pbm2)))
  (error "Arguments to PBM-OR are not matching PBMs"))
 (make-pbm (pbm-raw? pbm1)
	   (map-vector
	    (lambda (row1 row2)
	     (map-vector (lambda (bit1 bit2) (or bit1 bit2)) row1 row2))
	    (pbm-bitmap pbm1)
	    (pbm-bitmap pbm2))))

(define (pbm-not pbm)
 (unless (pbm? pbm) (error "Argument to PBM-NOT is not a PBM"))
 (make-pbm
  (pbm-raw? pbm)
  (map-vector (lambda (row) (map-vector not row)) (pbm-bitmap pbm))))

(define (pbm-xor pbm1 pbm2)
 (unless (and (pbm? pbm1)
            (pbm? pbm2)
            (= (pnm-width pbm1) (pnm-width pbm2))
            (= (pnm-height pbm1) (pnm-height pbm2)))
  (error "Arguments to PBM-XOR are not matching PBMs"))
 (make-pbm (pbm-raw? pbm1)
	   (map-vector
	    (lambda (row1 row2)
	     (map-vector (lambda (bit1 bit2) (xor bit1 bit2)) row1 row2))
	    (pbm-bitmap pbm1)
	    (pbm-bitmap pbm2))))

(define (pgm-absolute-difference pgm1 pgm2)
 (unless (and (pgm? pgm1)
            (pgm? pgm2)
            (= (pgm-maxval pgm1) (pgm-maxval pgm2))
            (= (pnm-width pgm1) (pnm-width pgm2))
            (= (pnm-height pgm1) (pnm-height pgm2)))
  (error "Arguments to PGM-ABSOLUTE-DIFFERENCE are not matching PGMs"))
 (make-pgm (pgm-raw? pgm1)
	   (pgm-maxval pgm1)
	   (map-vector
	    (lambda (row1 row2)
	     (map-vector (lambda (e1 e2) (abs (- e1 e2))) row1 row2))
	    (pgm-grey pgm1)
	    (pgm-grey pgm2))))

(define (empty-pnm? pnm)
 (cond
  ((ppm? pnm)
   (and (every-vector (lambda (row) (every-vector zero? row)) (ppm-red pnm))
      (every-vector (lambda (row) (every-vector zero? row)) (ppm-green pnm))
      (every-vector (lambda (row) (every-vector zero? row)) (ppm-blue pnm))))
  ((pgm? pnm)
   (every-vector (lambda (row) (every-vector zero? row)) (pgm-grey pnm)))
  ((pbm? pnm)
   (not
    (some-vector (lambda (row) (some-vector (lambda (a) a) row)) (pbm-bitmap pnm))))
  (else (error "Argument to EMPTY-PNM? is not a PNM"))))

(define (pbm-skeletonize pbm)
 (let ((height (pnm-height pbm)) (width (pnm-width pbm)))
  (let loop ((bitmap (pbm-bitmap pbm)))
   (let ((new-bitmap (map-vector (lambda (row) (map-vector (lambda (a) a) row)) bitmap)))
    (for-each-n
     (lambda (y)
      (for-each-n
       (lambda (x)
	(when (and
	       (matrix-ref bitmap y x)
	       (<= 4
		   (+ (if (and (not (zero? y))
			       (not (zero? x))
			       (matrix-ref bitmap (- y 1) (- x 1)))
			  1
			  0)
		      (if (and (not (zero? y))
			       (matrix-ref bitmap (- y 1) x))
			  1
			  0)
		      (if (and (not (zero? y))
			       (not (= x (- width 1)))
			       (matrix-ref bitmap (- y 1) (+ x 1)))
			  1
			  0)
		      (if (and (not (zero? x))
			       (matrix-ref bitmap y (- x 1)))
			  1
			  0)
		      (if (and (not (= x (- width 1)))
			       (matrix-ref bitmap y (+ x 1)))
			  1
			  0)
		      (if (and (not (= y (- height 1)))
			       (not (zero? x))
			       (matrix-ref bitmap (+ y 1) (- x 1)))
			  1
			  0)
		      (if (and (not (= y (- height 1)))
			       (matrix-ref bitmap (+ y 1) x))
			  1
			  0)
		      (if (and (not (= y (- height 1)))
			       (not (= x (- width 1)))
			       (matrix-ref bitmap (+ y 1) (+ x 1)))
			  1
			  0))
		   7))
	 (matrix-set! new-bitmap y x #f)))
       width))
     height)
    (if (equal? new-bitmap bitmap)
	(make-pbm (pbm-raw? pbm) bitmap)
	(loop new-bitmap))))))

(define (pbm-bloat pbm n)
 ;; Bloats white (#T) pixels.
 (let* ((height (pnm-height pbm))
	(width (pnm-width pbm))
	(bitmap (pbm-bitmap pbm))
	(new (make-matrix height width #f)))
  (do ((y 0 (+ y 1))) ((>= y height))
   (do ((x 0 (+ x 1))) ((>= x width))
    (do ((y0 (- y n) (+ y0 1))) ((> y0 (+ y n)))
     (when (and (>= y0 0) (< y0 height))
      (do ((x0 (- x n) (+ x0 1))) ((> x0 (+ x n)))
       (when (and (>= x0 0) (< x0 width) (matrix-ref bitmap y0 x0))
	(matrix-set! new y x #t)))))))
  (make-pbm (pbm-raw? pbm) new)))

(define (pbm-flood pbm point)
 (let ((new-bitmap
	(map-vector (lambda (row) (map-vector (lambda (a) a) row)) (pbm-bitmap pbm)))
       (height (pnm-height pbm))
       (width (pnm-width pbm)))
  (let loop ((point point))
   (when (and (<= 0 (y point) (- height 1)) (<= 0 (x point) (- width 1)))
    (unless (matrix-ref new-bitmap (y point) (x point))
     (matrix-set! new-bitmap (y point) (x point) #t)
     (loop (vector (- (x point) 1) (y point)))
     (loop (vector (+ (x point) 1) (y point)))
     (loop (vector (x point) (- (y point) 1)))
     (loop (vector (x point) (+ (y point) 1))))))
  (make-pbm (pbm-raw? pbm) new-bitmap)))

;;; PNM movies

(define (pnm-movie-frame-pathname pathname i)
 (when (string=? pathname "-") (error "Invalid pathname"))
 (let ((i (+ i 1)))
  (replace-extension (string-append (strip-extension pathname)
				    "_"
				    (cond ((< i 10) "0000")
					  ((< i 100) "000")
					  ((< i 1000) "00")
					  ((< i 10000) "0")
					  (else ""))
				    (number->string i))
		     (extension pathname))))

(define (pnm-movie-length pathname)
 (let loop ((i 0))
  (if (file-exists? (pnm-movie-frame-pathname pathname i))
      (loop (+ i 1))
      i)))

(define (read-pnm-movie pathname)
 (list->vector
  (map-n (lambda (i) (read-pnm (pnm-movie-frame-pathname pathname i)))
   (pnm-movie-length pathname))))

(define (write-pnm-movie pnm-movie pathname)
 (for-each-indexed
  (lambda (pnm i) (write-pnm pnm (pnm-movie-frame-pathname pathname i)))
  (vector->list pnm-movie)))

;;; Image operations

(define (overlay-pbm-on-pnm pbm pnm)
 ;; The white (#T) pixels in the pbm become white in the result.
 (unless (and (pbm? pbm)
            (or (pbm? pnm) (pgm? pnm) (ppm? pnm))
            (= (pnm-width pbm) (pnm-width pnm))
            (= (pnm-height pbm) (pnm-height pnm)))
  (error "Arguments to OVERLAY-PBM-ON-PNM are not a matching PBM and PNM"))
 (cond
  ((ppm? pnm)
   (make-ppm (pbm-raw? pbm)
	     (max *max-red* *max-green* *max-blue*)
	     (map-vector (lambda (bitmap-row red-row)
			  (map-vector (lambda (bit red)
				       (if bit *max-red* red))
				      bitmap-row
				      red-row))
			 (pbm-bitmap pbm)
			 (ppm-red pnm))
	     (map-vector (lambda (bitmap-row green-row)
			  (map-vector (lambda (bit green)
				       (if bit *max-green* green))
				      bitmap-row
				      green-row))
			 (pbm-bitmap pbm)
			 (ppm-green pnm))
	     (map-vector (lambda (bitmap-row blue-row)
			  (map-vector (lambda (bit blue)
				       (if bit *max-blue* blue))
				      bitmap-row
				      blue-row))
			 (pbm-bitmap pbm)
			 (ppm-blue pnm))))
  ((pgm? pnm)
   (make-pgm (pbm-raw? pbm)
	     *max-grey*
	     (map-vector (lambda (bitmap-row grey-row)
			  (map-vector (lambda (bit grey)
				       (if bit *max-grey* grey))
				      bitmap-row
				      grey-row))
			 (pbm-bitmap pbm)
			 (pgm-grey pnm))))
  ((pbm? pnm) (pbm-or pbm pnm))
  (else (fuck-up))))

(define (pnm-shift pnm delta)
 (cond
  ((ppm? pnm)
   (let* ((height (pnm-height pnm))
	  (width (pnm-width pnm))
	  (red (ppm-red pnm))
	  (green (ppm-green pnm))
	  (blue (ppm-blue pnm))
	  (new-red (make-matrix height width 0))
	  (new-green (make-matrix height width 0))
	  (new-blue (make-matrix height width 0))
	  (dy (y delta))
	  (dx (x delta)))
    (do ((y 0 (+ y 1))) ((>= y height))
     (when (and (>= (- y dy) 0) (< (- y dy) height))
      (do ((x 0 (+ x 1))) ((>= x width))
       (when (and (>= (- x dx) 0) (< (- x dx) width))
	(matrix-set! new-red y x (matrix-ref red (- y dy) (- x dx)))
	(matrix-set! new-green y x (matrix-ref green (- y dy) (- x dx)))
	(matrix-set! new-blue y x (matrix-ref blue (- y dy) (- x dx)))))))
    (make-ppm (ppm-raw? pnm) (ppm-maxval pnm) new-red new-green new-blue)))
  ((pgm? pnm)
   (let* ((height (pnm-height pnm))
	  (width (pnm-width pnm))
	  (grey (pgm-grey pnm))
	  (new-grey (make-matrix height width 0))
	  (dy (y delta))
	  (dx (x delta)))
    (do ((y 0 (+ y 1))) ((>= y height))
     (when (and (>= (- y dy) 0) (< (- y dy) height))
      (do ((x 0 (+ x 1))) ((>= x width))
       (when (and (>= (- x dx) 0) (< (- x dx) width))
	(matrix-set! new-grey y x (matrix-ref grey (- y dy) (- x dx)))))))
    (make-pgm (pgm-raw? pnm) (pgm-maxval pnm) new-grey)))
  ((pbm? pnm)
   (let* ((height (pnm-height pnm))
	  (width (pnm-width pnm))
	  (bitmap (pbm-bitmap pnm))
	  (new-bitmap (make-matrix height width #f))
	  (dy (y delta))
	  (dx (x delta)))
    (do ((y 0 (+ y 1))) ((>= y height))
     (when (and (>= (- y dy) 0) (< (- y dy) height))
      (do ((x 0 (+ x 1))) ((>= x width))
       (when (and (>= (- x dx) 0) (< (- x dx) width))
	(matrix-set! new-bitmap y x (matrix-ref bitmap (- y dy) (- x dx)))))))
    (make-pbm (pbm-raw? pnm) new-bitmap)))
  (else (error "Argument is not a PNM"))))

(define (pgm-smooth pgm sigma)
 (unless (pgm? pgm) (error "Argument to PGM-SMOOTH is not a PGM"))
 (let* ((height (pnm-height pgm))
	(width (pnm-width pgm))
	(grey1 (pgm-grey pgm))
	(grey2 (make-matrix height width 0)))
  (do ((y sigma (+ y 1))) ((= y (- height sigma)))
   (do ((x sigma (+ x 1))) ((= x (- width sigma)))
    (do ((i (- y sigma) (+ i 1))) ((= i (+ y sigma 1)))
     (do ((j (- x sigma) (+ j 1))) ((= j (+ x sigma 1)))
      (matrix-set!
       grey2 y x (+ (matrix-ref grey2 y x) (matrix-ref grey1 i j)))))
    (matrix-set!
     grey2 y x
     (inexact->exact
      (floor (/ (matrix-ref grey2 y x) (sqr (+ sigma sigma 1))))))))
  (make-pgm (pgm-raw? pgm) *max-grey* grey2)))

(define (pbm-proximity-clusterer pbm threshold)
 ;; Clusters white (#T) pixels.
 (unless (pbm? pbm) (fuck-up))
 (let* ((width (pnm-width pbm))
	(height (pnm-height pbm))
	(m (map-vector (lambda (row) (map-vector (lambda (a) a) row))
		       (pbm-bitmap pbm)))
	(i -1)
	(threshold-squared (sqr threshold)))
  (do ((y 0 (+ y 1))) ((= y height))
   (do ((x 0 (+ x 1))) ((= x width))
    (when (eq? (matrix-ref m y x) #t)
     (matrix-set! m y x i)
     (let loop ()
      (let ((again? #f))
       (do ((y1 0 (+ y1 1))) ((= y1 height))
	(do ((x1 0 (+ x1 1))) ((= x1 width))
	 (when (eqv? (matrix-ref m y1 x1) i)
	  (do ((y2 (max 0 (- y1 threshold)) (+ y2 1)))
	    ((= y2 (min height (+ y1 threshold 1))))
	   (do ((x2 (max 0 (- x1 threshold)) (+ x2 1)))
	     ((= x2 (min width (+ x1 threshold 1))))
	    (when (and (<= (+ (sqr (- x1 x2)) (sqr (- y1 y2)))
                        threshold-squared)
                     (eq? (matrix-ref m y2 x2) #t))
	     (matrix-set! m y2 x2 i)
	     (set! again? #t)))))))
       (when again? (loop))))
     (set! i (- i 1)))))
  (let ((clusters (make-vector (- (- i) 1) '())))
   (do ((y 0 (+ y 1))) ((= y height))
    (do ((x 0 (+ x 1))) ((= x width))
     (let ((i (matrix-ref m y x)))
      (when i
       (vector-set! clusters
		    (- (- i) 1)
		    (cons (vector x y) (vector-ref clusters (- (- i) 1))))))))
   (vector->list clusters))))

;;; Image pairs

(define (normal-flow-magnitude pgm1 pgm2 epsilon sigma sensitivity)
 (unless (and (pgm? pgm1)
            (pgm? pgm2)
            (= (pgm-maxval pgm1) (pgm-maxval pgm2))
            (eq? (pgm-raw? pgm1) (pgm-raw? pgm2))
            (= (pnm-width pgm1) (pnm-width pgm2))
            (= (pnm-height pgm1) (pnm-height pgm2)))
  (error "Arguments to NORMAL-FLOW-MAGNITUDE are not matching PGMs"))
 (let* ((width (pnm-width pgm1))
	(height (pnm-height pgm1))
	(e1 (pgm-grey (pgm-smooth pgm1 sigma)))
	(e2 (pgm-grey (pgm-smooth pgm2 sigma)))
	(m (make-matrix height width 0)))
  (do ((i 0 (+ i 1))) ((= i (- height 1)))
   (do ((j 0 (+ j 1))) ((= j (- width 1)))
    (let* ((ex (/ (- (+ (matrix-ref e1 (+ i 1) j)
			(matrix-ref e1 (+ i 1) (+ j 1))
			(matrix-ref e2 (+ i 1) j)
			(matrix-ref e2 (+ i 1) (+ j 1)))
		     (+ (matrix-ref e1 i j)
			(matrix-ref e1 i (+ j 1))
			(matrix-ref e2 i j)
			(matrix-ref e2 i (+ j 1))))
		  4.0))
	   (ey (/ (- (+ (matrix-ref e1 i (+ j 1))
			(matrix-ref e1 (+ i 1) (+ j 1))
			(matrix-ref e2 i (+ j 1))
			(matrix-ref e2 (+ i 1) (+ j 1)))
		     (+ (matrix-ref e1 i j)
			(matrix-ref e1 (+ i 1) j)
			(matrix-ref e2 i j)
			(matrix-ref e2 (+ i 1) j)))
		  4.0))
	   (et (/ (- (+ (matrix-ref e2 i j)
			(matrix-ref e2 i (+ j 1))
			(matrix-ref e2 (+ i 1) j)
			(matrix-ref e2 (+ i 1) (+ j 1)))
		     (+ (matrix-ref e1 i j)
			(matrix-ref e1 i (+ j 1))
			(matrix-ref e1 (+ i 1) j)
			(matrix-ref e1 (+ i 1) (+ j 1))))
		  4.0))
	   (l (sqrt (+ (sqr ex) (sqr ey)))))
     (matrix-set!
      m i j
      (min *max-grey*
	   (inexact->exact
	    (floor
	     (* *max-grey*
		(/ (if (< l epsilon) 0.0 (/ (abs et) l)) sensitivity)))))))))
  (pgm-smooth (make-pgm (pgm-raw? pgm1) *max-grey* m) sigma)))

(define (threshold-normal-flow-magnitude pgm1 pgm2 epsilon sigma threshold)
 ;; The moving pixels become white (#T) in the result.
 (unless (and (pgm? pgm1)
            (pgm? pgm2)
            (= (pgm-maxval pgm1) (pgm-maxval pgm2))
            (eq? (pgm-raw? pgm1) (pgm-raw? pgm2))
            (= (pnm-width pgm1) (pnm-width pgm2))
            (= (pnm-height pgm1) (pnm-height pgm2)))
  (error "Arguments to THRESHOLD-NORMAL-FLOW-MAGNITUDE are not matching PGMs"))
 (let* ((width (pnm-width pgm1))
	(height (pnm-height pgm1))
	(e1 (pgm-grey (pgm-smooth pgm1 sigma)))
	(e2 (pgm-grey (pgm-smooth pgm2 sigma)))
	(m (make-matrix height width #f)))
  (do ((i 0 (+ i 1))) ((= i (- height 1)))
   (do ((j 0 (+ j 1))) ((= j (- width 1)))
    (let* ((ex (/ (- (+ (matrix-ref e1 (+ i 1) j)
			(matrix-ref e1 (+ i 1) (+ j 1))
			(matrix-ref e2 (+ i 1) j)
			(matrix-ref e2 (+ i 1) (+ j 1)))
		     (+ (matrix-ref e1 i j)
			(matrix-ref e1 i (+ j 1))
			(matrix-ref e2 i j)
			(matrix-ref e2 i (+ j 1))))
		  4.0))
	   (ey (/ (- (+ (matrix-ref e1 i (+ j 1))
			(matrix-ref e1 (+ i 1) (+ j 1))
			(matrix-ref e2 i (+ j 1))
			(matrix-ref e2 (+ i 1) (+ j 1)))
		     (+ (matrix-ref e1 i j)
			(matrix-ref e1 (+ i 1) j)
			(matrix-ref e2 i j)
			(matrix-ref e2 (+ i 1) j)))
		  4.0))
	   (et (/ (- (+ (matrix-ref e2 i j)
			(matrix-ref e2 i (+ j 1))
			(matrix-ref e2 (+ i 1) j)
			(matrix-ref e2 (+ i 1) (+ j 1)))
		     (+ (matrix-ref e1 i j)
			(matrix-ref e1 i (+ j 1))
			(matrix-ref e1 (+ i 1) j)
			(matrix-ref e1 (+ i 1) (+ j 1))))
		  4.0))
	   (l (sqrt (+ (sqr ex) (sqr ey)))))
     (matrix-set!
      m i j (and (>= l epsilon) (>= (/ (abs et) l) threshold))))))
  (make-pbm (pgm-raw? pgm1) m)))

(define (pnm-black-window pnm upper-left lower-right)
 (cond ((ppm? pnm)
	(let* ((ppm (pnm-copy pnm))
	       (red (ppm-red ppm))
	       (green (ppm-green ppm))
	       (blue (ppm-blue ppm))
	       (yl (y upper-left))
	       (yh (y lower-right))
	       (xl (x upper-left))
	       (xh (x lower-right))
	       (height (pnm-height ppm))
	       (width (pnm-width ppm)))
	 (do ((y 0 (+ y 1))) ((>= y height))
	  (do ((x 0 (+ x 1))) ((>= x width))
	   (unless (and (>= y yl) (< y yh) (>= x xl) (< x xh))
	    (matrix-set! red y x 0)
	    (matrix-set! green y x 0)
	    (matrix-set! blue y x 0))))
	 ppm))
       ((pgm? pnm)
	(let* ((pgm (pnm-copy pnm))
	       (grey (pgm-grey pgm))
	       (yl (y upper-left))
	       (yh (y lower-right))
	       (xl (x upper-left))
	       (xh (x lower-right))
	       (height (pnm-height pgm))
	       (width (pnm-width pgm)))
	 (do ((y 0 (+ y 1))) ((>= y height))
	  (do ((x 0 (+ x 1))) ((>= x width))
	   (unless (and (>= y yl) (< y yh) (>= x xl) (< x xh))
	    (matrix-set! grey y x 0))))
	 pgm))
       ((pbm? pnm)
	(let* ((pbm (pnm-copy pnm))
	       (bitmap (pbm-bitmap pbm))
	       (yl (y upper-left))
	       (yh (y lower-right))
	       (xl (x upper-left))
	       (xh (x lower-right))
	       (height (pnm-height pbm))
	       (width (pnm-width pbm)))
	 (do ((y 0 (+ y 1))) ((>= y height))
	  (do ((x 0 (+ x 1))) ((>= x width))
	   (unless (and (>= y yl) (< y yh) (>= x xl) (< x xh))
	    (matrix-set! bitmap y x #f))))
	 pbm))
       (else (error "Argument is not a PNM"))))

(define (pnm-white-window pnm upper-left lower-right)
 (cond ((ppm? pnm)
	(let* ((ppm (pnm-copy pnm))
	       (maxval (ppm-maxval ppm))
	       (red (ppm-red ppm))
	       (green (ppm-green ppm))
	       (blue (ppm-blue ppm))
	       (yl (y upper-left))
	       (yh (y lower-right))
	       (xl (x upper-left))
	       (xh (x lower-right))
	       (height (pnm-height ppm))
	       (width (pnm-width ppm)))
	 (do ((y 0 (+ y 1))) ((>= y height))
	  (do ((x 0 (+ x 1))) ((>= x width))
	   (unless (and (>= y yl) (< y yh) (>= x xl) (< x xh))
	    (matrix-set! red y x maxval)
	    (matrix-set! green y x maxval)
	    (matrix-set! blue y x maxval))))
	 ppm))
       ((pgm? pnm)
	(let* ((pgm (pnm-copy pnm))
	       (maxval (pgm-maxval pgm))
	       (grey (pgm-grey pgm))
	       (yl (y upper-left))
	       (yh (y lower-right))
	       (xl (x upper-left))
	       (xh (x lower-right))
	       (height (pnm-height pgm))
	       (width (pnm-width pgm)))
	 (do ((y 0 (+ y 1))) ((>= y height))
	  (do ((x 0 (+ x 1))) ((>= x width))
	   (unless (and (>= y yl) (< y yh) (>= x xl) (< x xh))
	    (matrix-set! grey y x maxval))))
	 pgm))
       ((pbm? pnm)
	(let* ((pbm (pnm-copy pnm))
	       (bitmap (pbm-bitmap pbm))
	       (yl (y upper-left))
	       (yh (y lower-right))
	       (xl (x upper-left))
	       (xh (x lower-right))
	       (height (pnm-height pbm))
	       (width (pnm-width pbm)))
	 (do ((y 0 (+ y 1))) ((>= y height))
	  (do ((x 0 (+ x 1))) ((>= x width))
	   (unless (and (>= y yl) (< y yh) (>= x xl) (< x xh))
	    (matrix-set! bitmap y x #t))))
	 pbm))
       (else (error "Argument is not a PNM"))))

(define (pbm-ppm-and pbm ppm)
 ;; White (#T) pixels in the PBM are kept from the PPM. Black (#F) pixels in
 ;; the PBM become white in the result.
 (unless (and (pbm? pbm)
            (ppm? ppm)
            (= (pnm-width pbm) (pnm-width ppm))
            (= (pnm-height pbm) (pnm-height ppm)))
  (error "Arguments to PBM-PPM-AND are not matching PBMs"))
 (make-ppm (ppm-raw? ppm)
	   (ppm-maxval ppm)
	   (map-vector
	    (lambda (row red)
	     (map-vector (lambda (bit red) (if bit red *max-red*))
			 row red))
	    (pbm-bitmap pbm)
	    (ppm-red ppm))
	   (map-vector
	    (lambda (row green)
	     (map-vector (lambda (bit green) (if bit green *max-green*))
			 row green))
	    (pbm-bitmap pbm)
	    (ppm-green ppm))
	   (map-vector
	    (lambda (row blue)
	     (map-vector (lambda (bit blue) (if bit blue *max-blue*))
			 row blue))
	    (pbm-bitmap pbm)
	    (ppm-blue ppm))))

(define (pnm-rotate pnm)
 (cond ((pbm? pnm)
	(make-pbm (pbm-raw? pnm) (matrix-transpose (pbm-bitmap pnm))))
       ((pgm? pnm)
	(make-pgm (pgm-raw? pnm)
		  (pgm-maxval pnm)
		  (matrix-transpose (pgm-grey pnm))))
       ((ppm? pnm)
	(make-ppm (ppm-raw? pnm)
		  (ppm-maxval pnm)
		  (matrix-transpose (ppm-red pnm))
		  (matrix-transpose (ppm-green pnm))
		  (matrix-transpose (ppm-blue pnm))))
       (else (error "Argument is not a PNM"))))

(define (pnm-flip pnm)
 (cond ((pbm? pnm)
	(make-pbm (pbm-raw? pnm)
		  (list->vector (reverse (vector->list (pbm-bitmap pnm))))))
       ((pgm? pnm)
	(make-pgm (pgm-raw? pnm)
		  (pgm-maxval pnm)
		  (list->vector (reverse (vector->list (pgm-grey pnm))))))
       ((ppm? pnm)
	(make-ppm (ppm-raw? pnm)
		  (ppm-maxval pnm)
		  (list->vector (reverse (vector->list (ppm-red pnm))))
		  (list->vector (reverse (vector->list (ppm-green pnm))))
		  (list->vector (reverse (vector->list (ppm-blue pnm))))))
       (else (error "Argument is not a PNM"))))

(define (pgm-and-pbm pgm pbm)
 (map-vector (lambda (eg eb) (map-vector (lambda (eg eb) (if eb eg 0)) eg eb))
	     (pgm-grey pgm) (pbm-bitmap pbm)))

;;; Thresholding

(define (flatten-ppm ppm colour-transform)
 (map-vector
  (lambda (red-row green-row blue-row)
   (map-vector
    (lambda (r g b)
     (colour-transform `#(,r ,g ,b)))
    red-row green-row blue-row))
  (ppm-red ppm) (ppm-green ppm) (ppm-blue ppm)))

(define (binary-threshold pgm threshold)
 (unless (pgm? pgm) (panic "Argument is not a PGM"))
 (make-pbm #f (map-vector
	       (lambda (row)
		(map-vector
		 (lambda (col) (> col threshold)) row))
	       (pgm-grey pgm))))

;; Optimal Threshold - Otsu's Method
(define (find-threshold-otsu pgm)
 (let* ((normalised-histogram
	 (map-vector (lambda (v) (/ v (* (pnm-width pgm) (pnm-height pgm))))
		     (find-histogram (pgm-grey pgm) (pgm-maxval pgm))))
	(first-cumulative-moments
	 (cumulative-histogram
	  (map-indexed-vector (lambda (n i) (* n (+ i 1))) normalised-histogram)))
	(between-class-sigmas
	 (find-between-class-variances
	  (cumulative-histogram normalised-histogram)
	  first-cumulative-moments
	  (vector-ref first-cumulative-moments
		      (- (vector-length first-cumulative-moments) 1)))))
  (vector-positione (qreduce-vector max between-class-sigmas -inf.0)
                    between-class-sigmas)))

(define (binary-threshold-optimal pgm)
 (unless (pgm? pgm) (panic "Argument is not a PGM"))
 (binary-threshold pgm (find-threshold-otsu pgm)))

;;; Optimal Threshold - Method of Successive Means
(define (find-threshold-means pgm)
 (let* ((histogram (find-histogram (pgm-grey pgm) (pgm-maxval pgm)))
	(l (vector-length histogram)))
  (let loop ((t (/ l 2)) (oldt 0))
   (if (<= (abs (- t oldt)) 1)
       t
       (loop (inexact->exact
	      (round (/ (+ (histogram-mean (subvector histogram 0 t) 0)
			   (histogram-mean (subvector histogram t l) t))
			2)))
	     t)))))

(define (binary-threshold-means pgm)
 (unless (pgm? pgm) (panic "Argument is not a PGM"))
 (binary-threshold pgm (find-threshold-means pgm)))

;;; Colour Threshold
(define (sample-image ppm colour-tx window-centre window-size)
 (let* ((halfwin (inexact->exact (/ window-size 2)))
	(tlx (- (x window-centre) halfwin))
	(tly (- (y window-centre) halfwin)))
  (map-vector
   (lambda (r g b) (colour-tx `#(,r ,g ,b)))
   (unshape-matrix (crop (ppm-red ppm) tlx tly window-size window-size))
   (unshape-matrix (crop (ppm-green ppm) tlx tly window-size window-size))
   (unshape-matrix (crop (ppm-blue ppm) tlx tly window-size window-size)))))

(define (binary-threshold-colour ppm colour-tx point threshold)
 (let* ((window-size 21)
	(colour-values (sample-image ppm colour-tx point window-size))
	(mu (vectors-mean colour-values))
	(isigma (invert-matrix (vectors-variance mu colour-values))))
  (make-pbm #f
	    (map-vector
	     (lambda (red-row green-row blue-row)
	      (map-vector
	       (lambda (r g b)
		(< (mahalanobis-distance (colour-tx `#(,r ,g ,b)) mu isigma)
		   threshold))
	       red-row
	       green-row
	       blue-row))
	     (ppm-red ppm)
	     (ppm-green ppm)
	     (ppm-blue ppm)))))

;;; Histogram Equalization
(define (histogram-equalise pgm)
 (let* ((w (pnm-width pgm))
	(h (pnm-height pgm))
	(cdf (cumulative-histogram
	      (find-histogram (pgm-grey pgm) (pgm-maxval pgm))))
	(min-cdf (qreduce-vector min cdf +inf.0)))
  (make-pgm (pgm-raw? pgm)
	    (pgm-maxval pgm)
	    (map-vector
	     (lambda (row)
	      (map-vector
	       (lambda (val)
		(inexact->exact
		 (round (* (- (pgm-maxval pgm) 1)
			   (/ (- (vector-ref cdf val) min-cdf)
			      (- (* w h) 1))))))
	       row))
	     (pgm-grey pgm)))))

(define (colour-threshold ppm colour-tx mu isigma threshold)
 (make-pbm #f
	   (map-vector
	    (lambda (red-row green-row blue-row)
	     (map-vector
	      (lambda (r g b)
	       (< (mahalanobis-distance (colour-tx `#(,r ,g ,b)) mu isigma)
		  threshold))
	      red-row
	      green-row
	      blue-row))
	    (ppm-red ppm)
	    (ppm-green ppm)
	    (ppm-blue ppm))))

;;; Adaptive Thresholding
(define (make-integral-matrix matrix)
 (let ((integral-matrix
	(make-matrix (matrix-rows matrix) (matrix-columns matrix) 0)))
  (map-n
   (lambda (i)
    (map-n
     (lambda (j)
      (let ((current-val (matrix-ref matrix i j)))
       (matrix-set!
	integral-matrix i j
	(cond ((and (zero? i) (zero? j)) current-val)
	      ((zero? i) (+ current-val
			    (matrix-ref integral-matrix i (- j 1))))
	      ((zero? j) (+ current-val
			    (matrix-ref integral-matrix (- i 1) j)))
	      (else (+ current-val
		       (matrix-ref integral-matrix i (- j 1))
		       (matrix-ref integral-matrix (- i 1) j)
		       (- (matrix-ref integral-matrix (- i 1) (- j 1)))))))))
     (matrix-columns matrix)))
   (matrix-rows matrix))
  integral-matrix))

(define (compute-integral-matrix-mu integral-matrix x y w)
 (let ((r (- (matrix-rows integral-matrix) 1))
       (c (- (matrix-columns integral-matrix) 1))
       (del1 (inexact->exact (/ w 2)))
       (del2 (inexact->exact (ceiling (/ w 2)))))
  (/ (- (+ (matrix-ref integral-matrix (min (+ x del1) r) (min (+ y del1) c))
	   (matrix-ref integral-matrix (max (- x del2) 0) (max (- y del2) 0)))
	(+ (matrix-ref integral-matrix (min (+ x del1) r) (max (- y del2) 0))
	   (matrix-ref integral-matrix (max (- x del2) 0) (min (+ y del1) c))))
     (* w w))))

(define (compute-integral-matrix-sigma squared-integral-matrix mu x y w)
 (sqrt (- (compute-integral-matrix-mu squared-integral-matrix x y w)
	  (* mu mu))))

(define (compute-adaptive-threshold
	 integral-matrix squared-integral-matrix x y w)
 (let ((k 0.2) (R 128)
       (mu (compute-integral-matrix-mu integral-matrix x y w)))
  (* mu (+ 1 (* k (- (/ (compute-integral-matrix-sigma
			squared-integral-matrix mu x y w) R)
		    1))))))

(define (adaptive-threshold pgm winsize)
 (let* ((pixmap (pgm-grey pgm))
	(integral-pixmap (make-integral-matrix pixmap))
	(squared-integral-pixmap
	 (make-integral-matrix (map-matrix sqr pixmap))))
  (make-pbm (pgm-raw? pgm)
	    (map-n-vector
	     (lambda (i)
	      (map-n-vector
	       (lambda (j)
		(> (matrix-ref pixmap i j)
		   (compute-adaptive-threshold
		    integral-pixmap squared-integral-pixmap i j winsize)))
	       (matrix-columns pixmap)))
	     (matrix-rows pixmap)))))

(define (pgm-mean image)
 (histogram-mean (find-histogram (pgm-grey image) (pgm-maxval image)) 0))

(define (pgm-variance image)
 (let ((hist (find-histogram (pgm-grey image) (pgm-maxval image))))
  (histogram-variance hist (histogram-mean hist 0) 0)))

;; no wrapping, #F instead of the element at the edges
(define (slide-window m f size)
 (map-n-vector
  (lambda (x)
   (map-n-vector
    (lambda (y)
     (f (submatrix m x y size size)))
    (matrix-columns m)))
  (matrix-rows m)))

(define (count-pixels a)
 (qreduce-vector
  +
  (map-vector
   (lambda (e)
    (qreduce-vector (lambda (a b) (+ (if (and a (> a 0)) 1 0)
			       (if (and b (> b 0)) 1 0))) e 0)) a) 0))

(define (ppm->label-closest ppm mu1 sigma1 mu2 sigma2 #!optional (colour-transform identity))
 (make-pbm
  #f
  (map-vector
   (lambda (red-row green-row blue-row)
    (map-vector
     (lambda (r g b)
      (< (mahalanobis-distance (colour-transform `#(,r ,g ,b)) mu1 sigma1)
	 (mahalanobis-distance (colour-transform `#(,r ,g ,b)) mu2 sigma2)))
     red-row
     green-row
     blue-row))
   (ppm-red ppm)
   (ppm-green ppm)
   (ppm-blue ppm))))

;;; Histogram

(define (find-histogram pixmap maxval)
 (let ((bin (make-vector (+ maxval 1) 0)))
  (for-each-vector
   (lambda (row)
    (for-each-vector
     (lambda (col)
      (vector-set! bin col (+ (vector-ref bin col) 1)))
     row))
   pixmap)
  bin))

(define (histogram-mean histogram i)
 (let ((n (qreduce-vector + histogram 0)))
  (/ (qreduce-vector
      + (map-indexed-vector (lambda (n j) (* n (+ i j))) histogram) 0) n)))

(define (histogram-variance histogram mu i)
 (let ((n (qreduce-vector + histogram 0)))
  (/ (qreduce-vector
      + (map-indexed-vector (lambda (n j) (* (sqr (- (+ i j) mu)) n)) histogram) 0)
     n)))

(define (normalised-histogram histogram val)
 (map-vector (lambda (v) (/ v val)) histogram))

(define (weighted-histogram histogram)
 (map-indexed-vector (lambda (n i) (* n (+ i 1))) histogram))

(define (cumulative-histogram histogram)
 (let ((h (map-vector (lambda (v) v) histogram)))
  (let loop ((i 1))
   (when (< i (vector-length h))
    (vector-set! h i (+ (vector-ref h (- i 1)) (vector-ref h i)))
    (loop (+ i 1))))
  h))

(define (find-between-class-variances omegas mus mu-total)
 (map-vector
  (lambda (omega mu)
   (if (or (= omega 0) (= omega 1)) 0 (/ (sqr (- (* mu-total omega) mu)) (* omega (- 1 omega)))))
  omegas
  mus))

;;; Rendering Line Segments

(define (midpoint l) (k*v 0.5 (v+ (p l) (q l))))

;;; needs work: screen coordinates
(define (orientation v) (atan (y v) (x v)))

(define (line-segment-orientation l) (orientation (v- (q l) (p l))))

(define (line-segment->points l)
 (let ((n (ceiling (line-segment-length l))) (v (v- (q l) (p l))))
  (if (zero? n)
      (list (p l))
      (map-n (lambda (i) (v+ (p l) (k*v (/ i n) v))) (+ n 1)))))

(define (line-segments->points ls)
 (qmap-reduce append '() line-segment->points ls))

;;; Points and bounding boxes

(define (points->pbm-of-size points height width)
 ;; Takes a list of what will be the white (#T) points in the pbm image.
 (let ((bitmap (make-matrix height width #f)))
  (for-each (lambda (point) (matrix-set! bitmap (y point) (x point) #t))
   points)
  (make-pbm #t bitmap)))

(define (points->line-segments ps)
 (map make-line-segment (but-last ps) (rest ps)))

(define (points->bounding-box-points points)
 (let* ((xs (map x points))
	(ys (map y points))
	(min-x (minimume xs))
	(max-x (maximume xs))
	(min-y (minimume ys))
	(max-y (maximume ys)))
  (line-segments->points
   (list (make-line-segment (vector min-x min-y) (vector max-x min-y))
	 (make-line-segment (vector min-x min-y) (vector min-x max-y))
	 (make-line-segment (vector max-x max-y) (vector max-x min-y))
	 (make-line-segment (vector max-x max-y) (vector min-x max-y))))))

(define (bounding-box-size bb)
 (v- (vector (vector-ref bb 2) (vector-ref bb 3))
     (vector (vector-ref bb 0) (vector-ref bb 1))))

(define (points-bounding-box points)
 (vector (minimume (map x points)) (minimume (map y points))
	 (maximume (map x points)) (maximume (map y points))))

(define (points->points-bounding-box points bb)
 (map (lambda (p) (vector (- (x p) (x bb)) (- (y p) (y bb)))) points))

(define (normalize-to-bounding-box ps)
 (let* ((bb (points-bounding-box ps))
	(bottom-left (vector (vector-ref bb 0) (vector-ref bb 1)))
	(top-right (vector (vector-ref bb 2) (vector-ref bb 3)))
	(size (v- top-right bottom-left)))
  (map (lambda (p) (let ((p (v- p bottom-left)))
	       (vector (/ (x p) (x size)) (/ (y p) (y size)))))
       ps)))

(define (normalize-to-other-bounding-box points ps)
 (let* ((bb (points-bounding-box ps))
	(bottom-left (vector (vector-ref bb 0) (vector-ref bb 1)))
	(top-right (vector (vector-ref bb 2) (vector-ref bb 3)))
	(size (v- top-right bottom-left)))
  (map (lambda (p) (let ((p (v- p bottom-left)))
	       (vector (/ (x p) (x size)) (/ (y p) (y size)))))
       points)))

(define (points->target-bounding-box points target-bb)
 (let* ((size (bounding-box-size target-bb)))
  (map (lambda (v) (v+ (map-vector * v size)
		  (vector (x target-bb) (y target-bb))))
       (normalize-to-bounding-box points))))

(define (points->other-target-bounding-box points ps target-bb)
 (let* ((size (bounding-box-size target-bb)))
  (map (lambda (v) (v+ (map-vector * v size)
		  (vector (x target-bb) (y target-bb))))
       (normalize-to-other-bounding-box points ps))))

(define (bounding-box-bloat bb p)
 (let ((xl (vector-ref bb 0))
       (yl (vector-ref bb 1))
       (xu (vector-ref bb 2))
       (yu (vector-ref bb 3)))
  (vector (max 0 (quantize-coordinate (- xl (* p (- xu xl)))))
	  (max 0 (quantize-coordinate (- yl (* p (- yu yl)))))
	  (max 0 (quantize-coordinate (+ xu (* p (- xu xl)))))
	  (max 0 (quantize-coordinate (+ yu (* p (- yu yl))))))))

(define (bounding-box-crop bb image)
 (let ((x1 (max 0 (vector-ref bb 0)))
       (x2 (max 0 (min (vector-ref bb 2) (pnm-width image))))
       (y1 (max 0 (vector-ref bb 1)))
       (y2 (max 0 (min (vector-ref bb 3) (pnm-height image)))))
  (crop-image image x1 y1 (- x2 x1) (- y2 y1))))

(define (pbm->points pbm)
 (qmap-reduce-n
  append
  '()
  (lambda (y)
   (removeq
    #f
    (map-n (lambda (x) (if (matrix-ref (pbm-bitmap pbm) y x) (vector x y) #f))
	   (pnm-width pbm))))
  (pnm-height pbm)))

(define (points->pbm points height width)
 (let ((m (make-matrix height width #f)))
  (for-each (lambda (point)
	     (let ((y (quantize-coordinate (y point)))
		   (x (quantize-coordinate (x point))))
	      (when (and (>= y 0) (>= x 0) (< y height) (< x width))
	       (matrix-set! m y x #t))))
	    points)
  (make-pbm #f m)))

;;; Quantization

(define (quantize-coordinate x) (inexact->exact (round x)))

(define (quantize-point p)
 (vector (quantize-coordinate (x p)) (quantize-coordinate (y p))))

(define (quantize-points ps) (remove-duplicatese (map quantize-point ps)))

(define (quantize-line-segment l)
 (make-line-segment (quantize-point (p l)) (quantize-point (q l))))

(define (quantize-line-segments ls)
 (remove-duplicatese (map quantize-line-segment ls)))

;;; Ellipses

(define-structure ellipse x0 y0 t0 a b)

(define (ellipse-center ellipse)
 (vector (ellipse-x0 ellipse) (ellipse-y0 ellipse)))

(define (ellipse-area ellipse)
 (* pi (ellipse-a ellipse) (ellipse-b ellipse)))

(define (ellipse-eccentricity ellipse)
 (/ (ellipse-a ellipse) (ellipse-b ellipse)))

(define (radial-distance theta phi) (normalize-rotation (- phi theta)))

(define (point-on-ellipse? p ellipse tolerance)
 (let* ((p0 (vector (ellipse-x0 ellipse) (ellipse-y0 ellipse)))
	(r (rotation-matrix-2d (- (ellipse-t0 ellipse))))
	(a (ellipse-a ellipse))
	(b (ellipse-b ellipse))
	(q (unit (m*v r (v- p p0)))))
  (<= (abs (- (distance p p0) (magnitude (vector (* a (x q)) (* b (y q))))))
      tolerance)))

(define (line-segment->ellipse l)
 (let ((m (midpoint l)) (a (line-segment-length l)))
  ;; hardwired
  (make-ellipse (x m) (y m) (line-segment-orientation l) (* 0.5 a) (* 0.15 a))))

(define (ellipse->points e #!optional (n 360))
 (let* ((x0 (ellipse-x0 e))
	(y0 (ellipse-y0 e))
	(t0 (ellipse-t0 e))
	(a (ellipse-a e))
	(b (ellipse-b e))
	(rxx (cos t0))
	(rxy (- (sin t0)))
	(ryx (- rxy))
	(ryy rxx))
  (map-n (lambda (i)
	  (let ((ellipse-x (* a (sin (degrees->radians (* i (/ 360.0 n))))))
		(ellipse-y (* b (cos (degrees->radians (* i (/ 360.0 n)))))))
	   (vector (+ (* rxx ellipse-x) (* rxy ellipse-y) x0)
		   (+ (* ryx ellipse-x) (* ryy ellipse-y) y0))))
	 n)))

;;; Resizing

(define (resize-image w h i)
 (with-temporary-file
  "/tmp/resize.ppm"
  (lambda (filename)
   (write-pnm i filename)
   (system (format #f "mogrify -resize ~ax~a! ~a &> /dev/null" w h filename))
   (read-pnm filename))))

(define (subsample-pbm pbm)
 ;; hardwired to a factor of 2
 (make-pbm
  (pbm-raw? pbm)
  (map-n-vector
   (lambda (y)
    (map-n-vector
     (lambda (x)
      (or (matrix-ref (pbm-bitmap pbm) (* 2 y) (* 2 x))
	  (matrix-ref (pbm-bitmap pbm) (+ (* 2 y) 1) (* 2 x))
	  (matrix-ref (pbm-bitmap pbm) (* 2 y) (+ (* 2 x) 1))
	  (matrix-ref (pbm-bitmap pbm) (+ (* 2 y) 1) (+ (* 2 x) 1))))
     (quotient (pnm-width pbm) 2)))
   (quotient (pnm-height pbm) 2))))

(define (scale-ppm ppm scale)
 (cond ((= scale 1) ppm)
       ((= scale 2)
	(let* ((rows (matrix-rows (ppm-red ppm)))
	       (columns (matrix-columns (ppm-red ppm)))
	       (r (make-matrix (* rows scale) (* columns scale) 0))
	       (g (make-matrix (* rows scale) (* columns scale) 0))
	       (b (make-matrix (* rows scale) (* columns scale) 0))
	       (ir (ppm-red ppm))
	       (ig (ppm-green ppm))
	       (ib (ppm-blue ppm)))
	 (define (fill-4! m i j v)
	  (matrix-set! m i j v)
	  (matrix-set! m i (+ j 1) v)
	  (matrix-set! m (+ i 1) j v)
	  (matrix-set! m (+ i 1) (+ j 1) v))
	 (for-each-n
	  (lambda (i)
	   (let ((i2 (* i 2)))
	    (for-each-n
	     (lambda (j)
	      (let ((j2 (* j 2)))
	       (fill-4! r i2 j2 (matrix-ref ir i j))
	       (fill-4! g i2 j2 (matrix-ref ig i j))
	       (fill-4! b i2 j2 (matrix-ref ib i j))))
	     columns)))
	  rows)
	 (make-ppm (ppm-raw? ppm) (ppm-maxval ppm) r g b)))
       (else (panic "Scale unsupported"))))

;;; Connected Components

;; TODO defined private for now, should be renamed

(define-private-structure vertex pixels vertex edges)
(define-private-structure edge u v)
(define-private-structure graph vertices edges)

(define (pbm->graph pbm delta)
 ;; This uses max(delta-y, delta-x) distance metric.
 ;; needs work: a variant that uses Euclidean distance.
 (let* ((bitmap (pbm-bitmap pbm))
	(height (pnm-height pbm))
	(width (pnm-width pbm))
	(grid (map-n-vector
	       (lambda (i)
		(map-n-vector
		 (lambda (j)
		  (if (matrix-ref bitmap i j)
		      (make-vertex (list (vector j i)) #f '())
		      #f))
		 width))
	       height))
	(vertices '())
	(edges '()))
  (for-each-vector
   (lambda (row)
    (for-each-vector
     (lambda (element) (when element (set! vertices (cons element vertices))))
     row))
   grid)
  (for-each-n
   (lambda (i)
    (for-each-n
     (lambda (j)
      (when (matrix-ref grid i j)
       (for-each-m-n
	(lambda (i1)
	 (for-each-m-n
	  (lambda (j1)
	   (when (matrix-ref grid i1 j1)
	    (unless (and (= i i1) (= j j1))
	     (set! edges
		   (cons (make-edge (matrix-ref grid i j)
				    (matrix-ref grid i1 j1))
			 edges)))))
	  (max (- j delta) 0)
	  (min (+ j delta) (- width 1))))
	(max (- i delta) 0)
	(min (+ i delta) (- height 1)))))
     width))
   height)
  (make-graph vertices edges)))

(define (labeling->graph labeling delta)
 ;; This uses a Manhattan distance metric, but scans all matches within
 ;; a max(dx, dy) = delta distance
 ;; needs work: a variant that uses Euclidean distance.
 (let* ((height (matrix-rows labeling))
	(width (matrix-columns labeling))
	(grid (map-n-vector
	       (lambda (i)
		(map-n-vector
		 (lambda (j)
		  (make-vertex (list (vector j i)) #f '()))
		 width))
	       height))
	(vertices '())
	(edges '()))
  (for-each-vector
   (lambda (row)
    (for-each-vector
     (lambda (element) (set! vertices (cons element vertices)))
     row))
   grid)
  (for-each-n
   (lambda (i)
    (for-each-n
     (lambda (j)
      (when (matrix-ref grid i j)
       (for-each-m-n
	(lambda (i1)
	 (for-each-m-n
	  (lambda (j1)
	   (when (= (matrix-ref labeling i1 j1) (matrix-ref labeling i j))
	    (unless (or (and (= i i1) (= j j1))
			(> (+ (abs (- i i1)) (abs (- j j1))) delta))
	     (set! edges
		   (cons (make-edge (matrix-ref grid i j)
				    (matrix-ref grid i1 j1))
			 edges)))))
	  (max (- j delta) 0)
	  (min (+ j delta) (- width 1))))
	(max (- i delta) 0)
	(min (+ i delta) (- height 1)))))
     width))
   height)
  (make-graph vertices edges)))

(define (dereference-vertex u)
 (if (vertex-vertex u)
     (let ((v (dereference-vertex (vertex-vertex u)))) (set-vertex-vertex! u v))
     u))

(define (connected-components g)
 (for-each (lambda (u)
	    (set-vertex-vertex! u #f)
	    (set-vertex-edges! u '()))
	   (graph-vertices g))
 (for-each (lambda (e)
	    (let ((u (dereference-vertex (edge-u e)))
		  (v (dereference-vertex (edge-v e))))
	     (unless (eq? u v) (set-vertex-vertex! u v))))
	   (graph-edges g))
 (let ((classes
	(transitive-equivalence-classes
	 (lambda (u v) (eq? (dereference-vertex u) (dereference-vertex v)))
	 (graph-vertices g))))
  (for-each
   (lambda (e)
    (when (eq? (dereference-vertex (edge-u e)) (dereference-vertex (edge-v e)))
     (set-vertex-edges!
      (dereference-vertex (edge-u e))
      (cons e (vertex-edges (dereference-vertex (edge-u e)))))))
   (graph-edges g))
  (map (lambda (class)
	(make-graph class (vertex-edges (dereference-vertex (first class)))))
       classes)))

(define (vertices->pbm vertices height width)
 (let ((bitmap (make-matrix height width #f)))
  (for-each (lambda (u)
	     (for-each (lambda (pixel)
			(let ((y (quantize-coordinate (y pixel)))
			      (x (quantize-coordinate (x pixel))))
			 (when (and (>= y 0) (>= x 0) (< y height) (< x width))
			  (matrix-set! bitmap y x #t))))
		       (vertex-pixels u)))
	    vertices)
  (make-pbm #t bitmap)))

(define (graph->pbm g height width)
 (vertices->pbm (graph-vertices g) height width))

;;; Chains

(define (pbm->chains pbm)
 ;; Pixels at junctions can be in more than one chain.
 (let* ((m (pbm-bitmap pbm))
	(l (make-matrix (matrix-rows m) (matrix-columns m) '()))
	(chains '()))
  (define (trace-it i j)
   (let ((i-next (y (first (matrix-ref l i j))))
	 (j-next (x (first (matrix-ref l i j)))))
    (matrix-set! l i j (removee (vector j-next i-next) (matrix-ref l i j)))
    (cons
     (vector j i)
     (let loop ((i i-next) (j j-next) (i-previous i) (j-previous j))
      (unless (member (vector j-previous i-previous) (matrix-ref l i j))
       (fuck-up))
      (case (length (matrix-ref l i j))
       ((0) (fuck-up))
       ((2) (let ((entry (removee (vector j-previous i-previous)
				 (matrix-ref l i j))))
	     (matrix-set! l i j '())
	     (cons (vector j i)
		   (loop (y (first entry)) (x (first entry)) i j))))
       (else (matrix-set! l i j
			  (removee (vector j-previous i-previous)
				  (matrix-ref l i j)))
	     (list (vector j i))))))))
  (for-each-n
   (lambda (i)
    (for-each-n
     (lambda (j)
      (when (matrix-ref m i j)
       (when (and (> i 0) (> j 0) (matrix-ref m (- i 1) (- j 1)))
	(matrix-set! l i j (cons (vector (- j 1) (- i 1)) (matrix-ref l i j))))
       (when (and (> i 0)
		  (< j (- (matrix-columns m) 1))
		  (matrix-ref m (- i 1) (+ j 1)))
	(matrix-set! l i j (cons (vector (+ j 1) (- i 1)) (matrix-ref l i j))))
       (when (and (< i (- (matrix-rows m) 1))
		  (> j 0)
		  (matrix-ref m (+ i 1) (- j 1)))
	(matrix-set! l i j (cons (vector (- j 1) (+ i 1)) (matrix-ref l i j))))
       (when (and (< i (- (matrix-rows m) 1))
		  (< j (- (matrix-columns m) 1))
		  (matrix-ref m (+ i 1) (+ j 1)))
	(matrix-set! l i j (cons (vector (+ j 1) (+ i 1)) (matrix-ref l i j))))
       (when (and (> i 0) (matrix-ref m (- i 1) j))
	(matrix-set! l i j (cons (vector j (- i 1)) (matrix-ref l i j))))
       (when (and (> j 0) (matrix-ref m i (- j 1)))
	(matrix-set! l i j (cons (vector (- j 1) i) (matrix-ref l i j))))
       (when (and (< i (- (matrix-rows m) 1)) (matrix-ref m (+ i 1) j))
	(matrix-set! l i j (cons (vector j (+ i 1)) (matrix-ref l i j))))
       (when (and (< j (- (matrix-columns m) 1)) (matrix-ref m i (+ j 1)))
	(matrix-set! l i j (cons (vector (+ j 1) i) (matrix-ref l i j))))))
     (matrix-columns m)))
   (matrix-rows m))
  (let loop ()
   (cond
    ((some-n (lambda (i)
	      (some-n (lambda (j) (= (length (matrix-ref l i j)) 1))
		      (matrix-columns l)))
	     (matrix-rows l))
     (for-each-n (lambda (i)
		  (for-each-n (lambda (j)
			       (when (= (length (matrix-ref l i j)) 1)
				(set! chains (cons (trace-it i j) chains))))
			      (matrix-columns l)))
		 (matrix-rows l))
     (loop))
    ((some-n (lambda (i)
	      (some-n (lambda (j) (not (null? (matrix-ref l i j))))
		      (matrix-columns l)))
	     (matrix-rows l))
     (let* ((i (find-if (lambda (i)
			 (some-n (lambda (j) (not (null? (matrix-ref l i j))))
				 (matrix-columns l)))
			(enumerate (matrix-rows l))))
	    (j (find-if (lambda (j) (not (null? (matrix-ref l i j))))
			(enumerate (matrix-columns l)))))
      (set! chains (cons (trace-it i j) chains))
      (loop)))))
  (append
   ;; This is because with the above method, singleton pixels, either initial
   ;; ones or one left after chain removal, will not be included in chains.
   (let ((m (pbm-bitmap
	     (pbm-and
	      pbm
	      (pbm-not
	       (chains->pbm chains (pnm-height pbm) (pnm-width pbm)))))))
    (qreduce
     append
     (map-n (lambda (i)
	     (qreduce
	      append
	      (map-n (lambda (j)
		      (if (matrix-ref m i j) (list (list (vector j i))) '()))
		     (matrix-columns m))
	      '()))
	    (matrix-rows m))
     '()))
   chains)))

(define (chains->pbm chains height width)
 (let ((m (make-matrix height width #f)))
  (for-each (lambda (chain)
	     (for-each (lambda (point)
			(let ((y (quantize-coordinate (y point)))
			      (x (quantize-coordinate (x point))))
			 (when (and (>= y 0) (>= x 0) (< y height) (< x width))
			  (matrix-set! m y x #t))))
		       chain))
	    chains)
  (make-pbm #t m)))

(define (chain-filter pbm threshold)
 (chains->pbm (remove-if (lambda (chain) (< (length chain) threshold))
			 (pbm->chains pbm))
	      (pnm-height pbm)
	      (pnm-width pbm)))

(define (break-chain chain l)
 (let loop ((chain chain) (chains '()))
  (if (< (length chain) l)
      (if (null? chain) chains (cons chain chains))
      (loop (sublist chain l (length chain))
	    (cons (sublist chain 0 l) chains)))))

(define (break-chains chains l)
 (qmap-reduce append '() (lambda (chain) (break-chain chain l)) chains))


;;; Conjuring

(define (connected-component-filter pbm delta threshold)
 (qmap-reduce
  pbm-or
  (pbm-constant (pnm-width pbm) (pnm-height pbm) #f)
  (lambda (points) (points->pbm points (pnm-height pbm) (pnm-width pbm)))
  (remove-if
   (lambda (points) (< (length points) threshold))
   (map (lambda (g)
	 (pbm->points (graph->pbm g (pnm-height pbm) (pnm-width pbm))))
	(connected-components (pbm->graph pbm delta))))))

(define (conjure pbms delta span threshold1 threshold2)
 (let* ((new-pbms (list->vector pbms))
	(pbms (list->vector pbms))
	(height (pnm-height (vector-ref pbms 0)))
	(width (pnm-width (vector-ref pbms 0))))
  (for-each-n
   (lambda (i)
    (format #t "Frame ~s~%" i)
    (for-each
     (lambda (g)
      (let ((points (pbm->points (graph->pbm g height width))))
       (format #t "~s point~a~%"
	       (length points)
	       (if (= (length points) 1) "" "s"))
       (when (> (length points) threshold1)
	(let ((pbm (points->pbm points height width)))
	 (do ((j 1 (+ j 1)))
	   ((or (= (+ i j) (vector-length pbms))
		(= j span)
		(> (count-if
		    (lambda (point)
		     (matrix-ref
		      (pbm-bitmap (vector-ref pbms j)) (y point) (x point)))
		    points)
		   threshold2))
	    (if (or (= (+ i j) (vector-length pbms)) (= j span))
		(begin
		 (if (= (+ i j) (vector-length pbms))
		     (format #t "End of movie~%")
		     (format #t "Span too long~%"))
		 #f)
		(begin
		 (cond
		  ((positive? (- j 2))
		   (format #t "Conjuring ~s point~a from frame ~s ~s frame~a~%"
			   (length points)
			   (if (= (length points) 1) "" "s")
			   i
			   (- j 2)
			   (if (= (- j 2) 1) "" "s")))
		  ((zero? (- j 2)) (format #t "zero~%"))
		  (else (format #t "negative~%")))
		 (for-each-n
		  (lambda (k)
		   (vector-set! new-pbms
				(+ i k 1)
				(pbm-or (vector-ref new-pbms (+ i k 1)) pbm)))
		  (- j 2))))))))))
     (connected-components (pbm->graph (vector-ref pbms i) delta))))
   (vector-length pbms))
  (vector->list new-pbms)))

;;; Distance Transform

(define (distance-transform pbm)
 ;; This does Manhattan distance.
 (let ((height (pnm-height pbm))
       (width (pnm-width pbm))
       (bitmap (pbm-bitmap pbm)))
  (let loop ((distance
	      (map-vector
	       (lambda (row)
		(map-vector (lambda (pixel) (if pixel 0 +inf.0)) row))
	       bitmap))
	     (closest
	      (map-n-vector
	       (lambda (i)
		(map-n-vector
		 (lambda (j) (if (matrix-ref bitmap i j) (vector j i) #f))
		 width))
	       height)))
   (let ((new-distance
	  (map-n-vector
	   (lambda (i)
	    (map-n-vector
	     (lambda (j)
	      (min (if (> i 0)
		       (+ (matrix-ref distance (- i 1) j) 1)
		       +inf.0)
		   (if (> j 0)
		       (+ (matrix-ref distance i (- j 1)) 1)
		       +inf.0)
		   (matrix-ref distance i j)
		   (if (< j (- width 1))
		       (+ (matrix-ref distance i (+ j 1)) 1)
		       +inf.0)
		   (if (< i (- height 1))
		       (+ (matrix-ref distance (+ i 1) j) 1)
		       +inf.0)))
	     width))
	   height))
	 (new-closest
	  (map-n-vector
	   (lambda (i)
	    (map-n-vector
	     (lambda (j)
	      (list-ref
	       (list (if (> i 0)
			 (matrix-ref closest (- i 1) j)
			 +inf.0)
		     (if (> j 0)
			 (matrix-ref closest i (- j 1))
			 +inf.0)
		     (matrix-ref closest i j)
		     (if (< j (- width 1))
			 (matrix-ref closest i (+ j 1))
			 +inf.0)
		     (if (< i (- height 1))
			 (matrix-ref closest (+ i 1) j)
			 +inf.0))
	       (positione (min (if (> i 0)
                                   (+ (matrix-ref distance (- i 1) j) 1)
                                   +inf.0)
                               (if (> j 0)
                                   (+ (matrix-ref distance i (- j 1)) 1)
                                   +inf.0)
                               (matrix-ref distance i j)
                               (if (< j (- width 1))
                                   (+ (matrix-ref distance i (+ j 1)) 1)
                                   +inf.0)
                               (if (< i (- height 1))
                                   (+ (matrix-ref distance (+ i 1) j) 1)
                                   +inf.0))
                          (list (if (> i 0)
                                    (+ (matrix-ref distance (- i 1) j) 1)
                                    +inf.0)
                                (if (> j 0)
                                    (+ (matrix-ref distance i (- j 1)) 1)
                                    +inf.0)
                                (matrix-ref distance i j)
                                (if (< j (- width 1))
                                    (+ (matrix-ref distance i (+ j 1)) 1)
                                    +inf.0)
                                (if (< i (- height 1))
                                    (+ (matrix-ref distance (+ i 1) j) 1)
                                    +inf.0)))))
	     width))
	   height)))
    (if (every-vector v= distance new-distance)
	(list distance closest)
	(loop new-distance new-closest))))))

(define (closest-transform-ref closest-transform p)
 (matrix-ref closest-transform (y p) (x p)))

#>
double *euclidean_1d_dt(double *f, unsigned n) {
  double *d = (double *)malloc(2*n*sizeof(double));
  int *v = (int *)malloc(n*sizeof(int));
  double *z = (double *)malloc((n+1)*sizeof(double));
  int k = 0;
  v[0] = 0;
  z[0] = -HUGE_VAL;
  z[1] = +HUGE_VAL;
  for (int q = 1; q <= n-1; q++) {
    double s  = ((f[q]+(q*q))-(f[v[k]]+(v[k]*v[k])))/(2*q-2*v[k]);
    while (s <= z[k]) {
      k--;
      s  = ((f[q]+(q*q))-(f[v[k]]+(v[k]*v[k])))/(2*q-2*v[k]);
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k+1] = +HUGE_VAL;
  }
  k = 0;
  for (int q = 0; q <= n-1; q++) {
    while (z[k+1] < q) {k++;}
    d[q] = (q-v[k])*(q-v[k]) + f[v[k]];
    d[q+n] = v[k];
  }
  free(v);
  free(z);
  return d;
}
<#

(define (euclidean-1d-dt v)
 (let* ((n (vector-length v))
	(c-v (malloc (* c-sizeof-double n)))
	(c-v (vector->c-inexact-array c-v v c-sizeof-double #t))
	(c-dt ((c-function c-pointer ("euclidean_1d_dt" c-pointer unsigned-int)) c-v n))
	(dt (c-inexact-array->vector c-dt c-sizeof-double (* 2 n) #t)))
  (free c-dt)
  (free c-v)
  (shape-matrix dt n)))

(define (euclidean-2d-dt m)
 (let* ((row-dt (map-vector (lambda (row) (euclidean-1d-dt row)) m))
	(row-dt-val (map-vector x row-dt))
	(row-dt-pos (map-matrix inexact->exact (map-vector y row-dt)))
	(column-dt (map-vector (lambda (col) (euclidean-1d-dt col)) (matrix-transpose row-dt-val)))
	(column-dt-val (matrix-transpose (map-vector x column-dt)))
	(column-dt-pos (map-matrix inexact->exact (map-vector y column-dt))))
  `#(,column-dt-val
     ,(matrix-transpose
       (map-indexed-vector
	(lambda (r j) (map-vector (lambda (i) `#(,i ,(vector-ref (vector-ref row-dt-pos i) j))) r))
	column-dt-pos)))))

#>
double *euclidean_1d_dt_vals(double *f, unsigned n) {
  double *d = (double *)malloc(n*sizeof(double));
  int *v = (int *)malloc(n*sizeof(int));
  double *z = (double *)malloc((n+1)*sizeof(double));
  int k = 0;
  v[0] = 0;
  z[0] = -HUGE_VAL;
  z[1] = +HUGE_VAL;
  for (int q = 1; q <= n-1; q++) {
    double s  = ((f[q]+(q*q))-(f[v[k]]+(v[k]*v[k])))/(2*q-2*v[k]);
    while (s <= z[k]) {
      k--;
      s  = ((f[q]+(q*q))-(f[v[k]]+(v[k]*v[k])))/(2*q-2*v[k]);
    }
    k++;
    v[k] = q;
    z[k] = s;
    z[k+1] = +HUGE_VAL;
  }
  k = 0;
  for (int q = 0; q <= n-1; q++) {
    while (z[k+1] < q) {k++;}
    d[q] = (q-v[k])*(q-v[k]) + f[v[k]];
  }
  free(v);
  free(z);
  return d;
}
<#

(define (euclidean-1d-dt-vals v)
 (let* ((n (vector-length v))
	(c-v (malloc (* c-sizeof-double n)))
	(c-v (vector->c-inexact-array c-v v c-sizeof-double #t))
	(c-dt ((c-function c-pointer ("euclidean_1d_dt_vals" c-pointer unsigned-int)) c-v n))
	(dt (c-inexact-array->vector c-dt c-sizeof-double n #t)))
  (free c-dt)
  (free c-v)
  dt))

(define (euclidean-2d-dt-vals m)
 (let* ((row-dt (map-vector (lambda (row) (euclidean-1d-dt-vals row)) m))
	(column-dt (map-vector (lambda (col) (euclidean-1d-dt-vals col)) (matrix-transpose row-dt))))
  (matrix-transpose column-dt)))

;;; Buffers

(define-structure pnm-buffer buffer width height pixfmt storage type)

(define (pnm->pixfmt p p4?)
 (cond ((pbm? p) "I")
       ((pgm? p) "I")
       ((ppm? p) (if p4? "BGR0" "RGB"))
       (else (fuck-up))))

(define (pixfmt->stride pixfmt)
 (cond ((equal? pixfmt "I") 1)
       ((equal? pixfmt "RGB") 3)
       ((equal? pixfmt "BGR0") 4)
       (else (fuck-up))))

(define (image-type->pixfmt t p4?)
 (cond ((eq? 'bilevel t) "I")
       ((eq? 'grayscale t) "I")
       ((eq? 'color t) (if p4? "BGR0" "RGB"))
       (else (fuck-up))))

(define (image-type->storage-size t p4?)
 (cond ((eq? 'bilevel t) 1)
       ((eq? 'grayscale t) 1)
       ((eq? 'color t) (if p4? 4 3))
       (else (fuck-up))))

(define (pnm->image-type p)
 (cond ((pbm? p) 'bilevel)
       ((pgm? p) 'grayscale)
       ((ppm? p) 'color)
       (else (fuck-up))))

(define (pnm->storage p) 8)

(define (pixfmt->red pixfmt)
 (cond ((equal? pixfmt "RGB") 0)
       ((equal? pixfmt "BGR0") 2)
       (else (fuck-up))))
(define (pixfmt->green pixfmt)
 (cond ((equal? pixfmt "RGB") 1)
       ((equal? pixfmt "BGR0") 1)
       (else (fuck-up))))
(define (pixfmt->blue pixfmt)
 (cond ((equal? pixfmt "RGB") 2)
       ((equal? pixfmt "BGR0") 0)
       (else (fuck-up))))

(define (pnm-buffer-size p p4?)
 (* (pnm-width p) (pnm-height p) (cond ((pbm? p) 1)
				       ((pgm? p) 1)
				       ((ppm? p) (if p4? 4 3))
				       (else (fuck-up)))))

(define (pnm->pnm-buffer! p . p4?)
 (let* ((p4? (if (null? p4?) #f (car p4?)))
        (buffer (malloc (pnm-buffer-size p p4?))))
  (pnm-fill-buffer! p buffer (pnm->pixfmt p p4?))
  (make-pnm-buffer
   buffer
   (pnm-width p)
   (pnm-height p)
   (pnm->pixfmt p p4?)
   (pnm->storage p)
   (pnm->image-type p))))

(define (pnm-fill-buffer! p b pixfmt)
 (let* ((width (pnm-width p))
	(height (pnm-height p))
        (stride (pixfmt->stride pixfmt))
	(poke
	 (cond
	  ((pbm? p)
	   (lambda (r c)
	    (c-byte-set!
	     b
	     (+ (* width r) c)
	     (if  (matrix-ref (pbm-bitmap p) r c) 255 0))))
	  ((pgm? p)
	   (if (= (pgm-maxval p) 255)
	       (lambda (r c)
		(c-byte-set!
		 b
		 (+ (* width r) c)
		 (matrix-ref (pgm-grey p) r c)))
	       (fuck-up)))
	  ((ppm? p)
	   (let ((red-index (pixfmt->red pixfmt))
                 (green-index (pixfmt->green pixfmt))
                 (blue-index (pixfmt->blue pixfmt)))
            (if (= (ppm-maxval p) 255)
                (lambda (r c)
                 (c-byte-set! b (+ (* stride (+ (* width r) c)) red-index) (matrix-ref (ppm-red p) r c))
                 (c-byte-set! b (+ (* stride (+ (* width r) c)) blue-index) (matrix-ref (ppm-blue p) r c))
                 (c-byte-set! b (+ (* stride (+ (* width r) c)) green-index) (matrix-ref (ppm-green p) r c)))
                (fuck-up)))))))
  (for-each-n (lambda (r) (for-each-n (lambda (c) (poke r c)) width)) height)))

(define (free-pnm-buffer! p) (free (pnm-buffer-buffer p)) #f)

(define (pnm-buffer->pnm b)
 ((case (pnm-buffer-type b)
   ((bilevel) pnm-buffer->pbm)
   ((grayscale) pnm-buffer->pgm)
   ((color) pnm-buffer->ppm)
   (else (fuck-up))) b))

(define (pnm-buffer->pbm p)
 ;; TODO check maxval
 (make-pbm
  #f
  (map-n-vector
   (lambda (r)
    (map-n-vector
     (lambda (c)
      (> (c-byte-ref
	  (pnm-buffer-buffer p)
	  (+ (* r (pnm-buffer-width p)) c)) 1))
     (pnm-buffer-width p)))
   (pnm-buffer-height p))))

(define (pnm-buffer->pgm p)
 ;; TODO check maxval
 (make-pgm
  #f
  255
  (map-n-vector
   (lambda (r)
    (map-n-vector
     (lambda (c)
      (c-byte-ref
       (pnm-buffer-buffer p)
       (+ (* r (pnm-buffer-width p)) c)))
     (pnm-buffer-width p)))
   (pnm-buffer-height p))))

(define (pnm-buffer->ppm p)
 ;; TODO check maxval
 (let ((ppm (make-ppm
	     #f
	     255
	     (make-matrix (pnm-buffer-height p) (pnm-buffer-width p) 0)
	     (make-matrix (pnm-buffer-height p) (pnm-buffer-width p) 0)
	     (make-matrix (pnm-buffer-height p) (pnm-buffer-width p) 0)))
       (stride (pixfmt->stride (pnm-buffer-pixfmt p)))
       (red-index (pixfmt->red (pnm-buffer-pixfmt p)))
       (green-index (pixfmt->green (pnm-buffer-pixfmt p)))
       (blue-index (pixfmt->blue (pnm-buffer-pixfmt p))))
  (for-each-n
    (lambda (r)
     (for-each-n
       (lambda (c)
        (matrix-set! (ppm-red ppm)
                     r c
                     (c-byte-ref (pnm-buffer-buffer p) (+ (* stride (+ (* r (pnm-buffer-width p)) c)) red-index)))
        (matrix-set! (ppm-green ppm)
                     r c
                     (c-byte-ref (pnm-buffer-buffer p) (+ (* stride (+ (* r (pnm-buffer-width p)) c)) green-index)))
        (matrix-set! (ppm-blue ppm)
                     r c
                     (c-byte-ref (pnm-buffer-buffer p) (+ (* stride (+ (* r (pnm-buffer-width p)) c)) blue-index))))
      (pnm-buffer-width p)))
   (pnm-buffer-height p))
  ppm))

(define (imlib-image->pnm-buffer! image)
 ;; FIXME This may leak memory, I'm not sure
 (make-pnm-buffer (image-get-data-for-reading-only (image-clone image))
                  (image-width image)
                  (image-height image)
                  "BGR0"
                  8
                  'color))

(define (imlib-image->ppm image)
 (pnm-buffer->ppm (imlib-image->pnm-buffer! image)))

(define (ppm->imlib-image ppm)
 (let ((buf (pnm->pnm-buffer! ppm #t)))
  ;; (let ((i (image-create-using-copied-data (pnm-buffer-width buf)
  ;;                                          (pnm-buffer-height buf)
  ;;                                          (pnm-buffer-buffer buf))))
  ;;  (set-finalizer! i gc-collect-image)
  ;;  (free-pnm-buffer! buf)
  ;;  i)
  (error "This function is commented pending a patch to the imlib2 egg")))

(define (image->pnm-buffer! image)
 (cond ((pnm? image) (pnm->pnm-buffer! image #f))
       ((image? image) (imlib-image->pnm-buffer! image))
       (else (fuck-up))))

;;; Misc

(define (ppm-absolute-difference ppm1 ppm2)
 ;; maxval hardwired to 255
 (unless (and (ppm? ppm1)
	      (ppm? ppm2)
	      (= (ppm-maxval ppm1) 255)
	      (= (ppm-maxval ppm2) 255)
	      (eq? (ppm-raw? ppm1) (ppm-raw? ppm2))
	      (= (pnm-width ppm1) (pnm-width ppm2))
	      (= (pnm-height ppm1) (pnm-height ppm2)))
  (panic "Arguments to PPM-ABSOLUTE-DIFFERENCE are not matching PPMs"))
 (let ((scale
	(/ 255 (distance (vector 0.0 0.0 0.0) (vector 255.0 255.0 255.0)))))
  (make-pgm (ppm-raw? ppm1)
	    255
	    (map-vector
	     (lambda (red1 green1 blue1 red2 green2 blue2)
	      (map-vector (lambda (red1 green1 blue1 red2 green2 blue2)
			   (inexact->exact
			    (round
			     (* scale
				(distance (vector red1 green1 blue1)
					  (vector red2 green2 blue2))))))
			  red1 green1 blue1 red2 green2 blue2))
	     (ppm-red ppm1)
	     (ppm-green ppm1)
	     (ppm-blue ppm1)
	     (ppm-red ppm2)
	     (ppm-green ppm2)
	     (ppm-blue ppm2)))))

(define (set-ppm-pixel! ppm x y value)
 (matrix-set! (ppm-red ppm) y x (r value))
 (matrix-set! (ppm-green ppm) y x (g value))
 (matrix-set! (ppm-blue ppm) y x (b value))
 ppm)

(define (pnm-pixel? i x y)
 (and (>= x 0) (< x (pnm-width i))
      (>= y 0) (< y (pnm-height i))))

(define (map-ppm-values ppm f)
 (make-ppm
  (ppm-raw? ppm)
  (ppm-maxval ppm)
  (map-matrix (lambda (a) (inexact->exact (f a))) (ppm-red ppm))
  (map-matrix (lambda (a) (inexact->exact (f a))) (ppm-green ppm))
  (map-matrix (lambda (a) (inexact->exact (f a))) (ppm-blue ppm))))

(define (ppm-burn base mask colour)
 (let ((image (pnm-copy base)))
  (for-each-indexed-vector
   (lambda (row j)
    (for-each-indexed-vector
     (lambda (val i)
      (if val (set-ppm-pixel! image i j colour)))
     row))
   (pbm-bitmap mask))
  image))

;; Stacking

(define (pbm-stack-vertical pbm1 pbm2)
 (make-pbm (pbm-raw? pbm1)
	   (list->vector (append (vector->list (pbm-bitmap pbm1))
				 (vector->list (pbm-bitmap pbm2))))))

(define (pbm-stack-horizontal pbm1 pbm2)
 (make-pbm
  (pbm-raw? pbm1)
  (map-vector (lambda (row1 row2)
	       (list->vector (append (vector->list row1) (vector->list row2))))
	      (pbm-bitmap pbm1)
	      (pbm-bitmap pbm2))))

(define (ppm-stack-vertical ppm1 ppm2)
 (make-ppm (ppm-raw? ppm1)
	   (ppm-maxval ppm1)
	   (list->vector (append (vector->list (ppm-red ppm1))
				 (vector->list (ppm-red ppm2))))
	   (list->vector (append (vector->list (ppm-green ppm1))
				 (vector->list (ppm-green ppm2))))
	   (list->vector (append (vector->list (ppm-blue ppm1))
				 (vector->list (ppm-blue ppm2))))))

(define (ppm-stack-horizontal ppm1 ppm2)
 (make-ppm
  (ppm-raw? ppm1)
  (ppm-maxval ppm1)
  (map-vector (lambda (row1 row2)
	       (list->vector (append (vector->list row1) (vector->list row2))))
	      (ppm-red ppm1) (ppm-red ppm2))
  (map-vector (lambda (row1 row2)
	       (list->vector (append (vector->list row1) (vector->list row2))))
	      (ppm-green ppm1) (ppm-green ppm2))
  (map-vector (lambda (row1 row2)
	       (list->vector (append (vector->list row1) (vector->list row2))))
	      (ppm-blue ppm1) (ppm-blue ppm2))))

;;; Angles

(define (degrees->radians angle) (* two-pi/360 angle))

(define (radians->degrees angle) (* three-sixty/two-pi angle))

(define (normalize-rotation rotation)
 (cond ((> rotation pi) (normalize-rotation (- rotation two-pi)))
       ((<= rotation minus-pi) (normalize-rotation (+ rotation two-pi)))
       (else rotation)))

(define (rotation+ x y) (normalize-rotation (+ x y)))

(define (rotation- x y) (normalize-rotation (- x y)))

(define (angle-separation x y)
 (min (abs (rotation- x y)) (abs (rotation- y x))))

(define (rotation-matrix-2d theta)
 (let ((ct (cos theta))
       (st (sin theta)))
  (vector (vector ct (- st)) (vector st ct))))

(define (mean-angle angles)
 (atan (qmap-reduce + 0 sin angles) (qmap-reduce + 0 cos angles)))

;;; Misc

(define (show-image i) (show i))

(define (show i)
 (define (with-temporary-file filename f)
  (f (create-temporary-file (if (has-extension? filename)
                                (extension filename)
                                ""))))
 (cond ((or (ppm? i) (pgm? i) (pbm? i))
        (with-temporary-file
         (cond ((ppm? i) "/tmp/show.ppm")
               ((pgm? i) "/tmp/show.pgm")
               ((pbm? i) "/tmp/show.pbm")
               (else (fuck-up)))
         (lambda (filename)
          (write-pnm i filename)
          (system (string-append "(feh --force-aliasing " filename "; rm " filename ")&")))))
       ;; should be imlib?
       ((image? i)
	(with-temporary-file
         "/tmp/show.png"
         (lambda (filename)
          (image-save i filename)
          (system (string-append "(feh --force-aliasing " filename "; rm " filename ")&")))))
       (else (error "don't know how to show" i))))
)
