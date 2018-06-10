#lang typed/racket
(provide graph/linear
         graph/exponential
         graph/log
         graph/power

         linear-fit-params linear-fit
         exp-fit-params exp-fit
         log-fit-params log-fit
         power-fit-params power-fit

         Fit try-fits best-fit
         fit->string fit->fun)

(require plot/no-gui math/flonum
         "polyfit.rkt")

;; math from http://mathworld.wolfram.com/LeastSquaresFitting.html

(define-type Flonums (Listof Nonnegative-Flonum))

(define-type Grapher
  (->* (Flonums Flonums)
       ((Option (Listof Flonum)))
       (Values renderer2d renderer2d renderer2d)))
(define-type Grapher/poly
  (->* (Flonums Flonums Positive-Integer)
       ((Option (Listof Flonum)))
       (Values renderer2d renderer2d renderer2d)))

(define-type Fitter (-> Flonums Flonums (-> Real Real)))
(define-type Fitter/poly (-> Flonums Flonums Positive-Integer (-> Real Real)))
(define-type Parameterizer (-> Real Real (-> Real Real)))
(define-type Parameterizer* (-> (Listof Real) (-> Real Real)))


(: graph/gen : (-> (Listof Nonnegative-Flonum) (Listof Nonnegative-Flonum)
                   (Option (Listof Flonum))
                   Fitter
                   (Values renderer2d renderer2d renderer2d)))
(define (graph/gen pts-x pts-y error fit)
  (values (points (map (inst list Nonnegative-Flonum) pts-x pts-y)
                  #:x-min (min 0 (apply min pts-x))
                  #:y-min (min 0 (apply min pts-y)))
          (function (fit pts-x pts-y))
          (error-bars (map (lambda ([x : Real] [y : Real] [δ : Real])
                             (list x y (* y δ)))
                           pts-x pts-y (or error (build-list (length pts-x) (const 0)))))))

(: Σ : (-> (Listof Flonum) Flonum))
(define (Σ l) (foldl + 0.0 l))


;; --------------------
;; | linear fit

(: linear-fit-params : (-> Flonums Flonums (Values Real Real)))
(define (linear-fit-params pts-x pts-y)
  (define len (length pts-x))
  (define Σx (Σ pts-x))
  (define Σy (Σ pts-y))
  (define Σxy (Σ (map * pts-x pts-y)))
  (define ΣxΣy (* Σx Σy))
  (define Σx^2 (Σ (map sqr pts-x)))
  (define slope
    (/ (- (* len Σxy) ΣxΣy)
       (- (* len Σx^2) (sqr Σx))))
  (define offset
    (/ (- Σy (* slope Σx))
       len))
  (values offset slope))

(: linear-fit-params/list : (-> Flonums Flonums (List Real Real)))
(define (linear-fit-params/list pts-x pts-y)
  (define-values (offset slope) (linear-fit-params pts-x pts-y))
  (list offset slope))

(: linear-fit/with-params : Parameterizer)
(define (linear-fit/with-params offset slope)
  (lambda ([x : Real]) (+ offset (* slope (fl x)))))

(: linear-fit : Fitter)
(define (linear-fit pts-x pts-y)
  (define-values [a b]
    (linear-fit-params pts-x pts-y))
  (linear-fit/with-params a b))

(: graph/linear : Grapher)
(define (graph/linear pts-x pts-y [error #f])
  (graph/gen pts-x pts-y error linear-fit))


;; --------------------
;; | exp fit

;; see http://mathworld.wolfram.com/LeastSquaresFittingExponential.html
(: exp-fit-params : (-> Flonums Flonums (Values Real Real)))
(define (exp-fit-params pts-x pts-y)
  (define lny (map log pts-y))
  (define x^2 (map sqr pts-x))
  (define Σx^2y  (Σ (map * x^2 pts-y)))
  (define Σylny  (Σ (map * pts-y lny)))
  (define Σxy    (Σ (map * pts-x pts-y)))
  (define Σxylny (Σ (map * pts-x pts-y lny)))
  (define Σy     (Σ pts-y))
  (define ΣyΣx^2y-Σxy^2 (- (* Σy Σx^2y) (sqr Σxy)))
  (define a (/ (- (* Σx^2y Σylny) (* Σxy Σxylny))
               ΣyΣx^2y-Σxy^2))
  (define b (/ (- (* Σy Σxylny) (* Σxy Σylny))
               ΣyΣx^2y-Σxy^2))
  (values (exp a) b))

(: exp-fit-params/list : (-> Flonums Flonums (List Real Real)))
(define (exp-fit-params/list pts-x pts-y)
  (define-values (a b) (exp-fit-params pts-x pts-y))
  (list a b))

(: exp-fit/with-params : Parameterizer)
(define (exp-fit/with-params coeff exp-coeff)
  (lambda ([x : Real]) (* coeff (exp (* exp-coeff (fl x))))))

(: exp-fit : Fitter)
(define (exp-fit pts-x pts-y)
  (define-values [A B]
    (exp-fit-params pts-x pts-y))
  (exp-fit/with-params A B))

(: graph/exponential : Grapher)
(define (graph/exponential pts-x pts-y [error #f])
  (graph/gen pts-x pts-y error exp-fit))


;; --------------------
;; | log fit

;; see http://mathworld.wolfram.com/LeastSquaresFittingLogarithmic.html
(: log-fit-params : (-> Flonums Flonums (Values Real Real)))
(define (log-fit-params pts-x pts-y)
  (define n (length pts-x))
  (define lnx (map fllog pts-x))
  (define Σy*lnx (flsum (map fl* pts-y lnx)))
  (define Σy (flsum pts-y))
  (define Σlnx (flsum lnx))
  (define lnx^2 (map sqr lnx))
  (define Σlnx^2 (flsum lnx^2))

  (define b (/ (- (* n Σy*lnx)
                  (* Σy Σlnx))
               (- (* n Σlnx^2)
                  (sqr Σlnx))))

  (define a (/ (- Σy (* b Σlnx))
               n))
  (values a b))

(: log-fit-params/list : (-> Flonums Flonums (List Real Real)))
(define (log-fit-params/list pts-x pts-y)
  (define-values (a b) (log-fit-params pts-x pts-y))
  (list a b))

(: log-fit/with-params : Parameterizer)
(define (log-fit/with-params offset coeff)
  (lambda ([x : Real])
    (define fx (fl x))
    (if (nnn? fx)
        (+ offset (* coeff (log fx)))
        +nan.0)))

(: log-fit : Fitter)
(define-predicate nnn? Nonnegative-Flonum)
(define (log-fit pts-x pts-y)
  (define-values [a b]
    (log-fit-params pts-x pts-y))
  (log-fit/with-params a b))

(: graph/log : Grapher)
(define (graph/log pts-x pts-y [error #f])
  (graph/gen pts-x pts-y error log-fit))


;; --------------------
;; | power fit

;; see http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
(: power-fit-params : (-> Flonums Flonums (Values Real Real)))
(define (power-fit-params pts-x pts-y)
  (define lnx (map log pts-x))
  (define lny (map log pts-y))
  (define n (length pts-x))
  (define Σlnx (Σ lnx))
  (define Σlny (Σ lny))
  (define Σlnxlny (Σ (map * lnx lny)))
  (define Σlnx^2 (Σ (map sqr lnx)))
  (define b (/ (- (* n Σlnxlny) (* Σlnx Σlny))
               (- (* n Σlnx^2) (sqr Σlnx))))
  (define a (/ (- Σlny (* b Σlnx)) n))
  (values a b))

(: power-fit-params/list : (-> Flonums Flonums (List Real Real)))
(define (power-fit-params/list pts-x pts-y)
  (define-values (a b) (power-fit-params pts-x pts-y))
  (list a b))

(: power-fit/with-params : Parameterizer)
(define (power-fit/with-params coeff exp-coeff)
  (lambda ([x : Real])
    (define fx (fl x))
    (if (pn? fx)
        (* (exp coeff) (expt fx exp-coeff))
        +nan.0)))

(: power-fit : Fitter)
(define-predicate pn? Positive-Flonum)
(define (power-fit pts-x pts-y)
  (define-values [a b]
    (power-fit-params pts-x pts-y))
  (power-fit/with-params a b))

(: graph/power : Grapher)
(define (graph/power pts-x pts-y [error #f])
  (graph/gen pts-x pts-y error power-fit))


;; --------------------
;; | Polynomial fit

;; See polyfit.rkt
(: poly-fit-params : (-> Flonums Flonums Positive-Integer (Listof Real)))
(define (poly-fit-params pts-x pts-y degree)
  (gradient-descent/auto-tune pts-x pts-y degree))

(: poly-fit/with-params : Parameterizer*)
(define (poly-fit/with-params coeffs)
  (lambda ([x : Real])
    (define fx (fl x))
    (define vals : (Listof Real)
      (for/list ([c : Real (in-list coeffs)]
                 [pow : Natural (in-naturals)])
        (* c (flexpt fx (fl pow)))))
    (apply + vals)))

(: poly-fit : Fitter/poly)
(define (poly-fit pts-x pts-y degree)
  (define coeffs (poly-fit-params pts-x pts-y degree))
  (poly-fit/with-params coeffs))

(: graph/poly : Grapher/poly)
(define (graph/poly pts-x pts-y degree [error #f])
  (define (poly-fit/this-degree [pts-x : Flonums]
                                [pts-y : Flonums])
    (poly-fit pts-x pts-y degree))
  (graph/gen pts-x pts-y error poly-fit/this-degree))


;; --------------------
;; | Exploratory fit
;; Try every fit function, returning a summary or just the best

;; For now, just always use SSE for error
;; ;; f(x) = offset + slope*x
;; (define-struct LinearFit ([offset : Real] [slope : Real] [sse : Real]))
;; ;; f(x) = coeff*e^(exp-coeff*x)
;; (define-struct ExpFit ([coeff : Real] [exp-coeff : Real] [sse : Real]))
;; ;; f(x) = offset + coeff*log(x)
;; (define-struct LogFit ([offset : Real] [coeff : Real] [sse : Real]))
;; ;; f(x) = (e^coeff)*(x^exp-coeff)
;; (define-struct PowerFit ([coeff : Real] [exp-coeff : Real] [sse : Real]))

(: point-error : (-> Real Real Real))
(define (point-error actual expected)
  (expt (- expected actual) 2))

(: sse : (-> Flonums Flonums (-> Real Real) Real))
(define (sse pts-x pts-y f)
  (foldl + 0
         (map point-error
              (map f pts-x)
              pts-y)))

(struct Fit ([params : (Listof Real)] [sse : Real] [type : Symbol]))

(: try-fits : (-> Flonums Flonums (Listof Fit)))
(define (try-fits pts-x pts-y
                  #:poly-max-degree [poly-max-degree 0])
  (define fits-to-try `((linear ,linear-fit-params/list ,linear-fit/with-params)
                        (exp ,exp-fit-params/list ,exp-fit/with-params)
                        (log ,log-fit-params/list ,log-fit/with-params)
                        (power ,power-fit-params/list ,power-fit/with-params)))
  (define fits
    (for/fold ([fits : (Listof Fit) empty])
              ([fit-to-try fits-to-try])
      (let* ([type (first fit-to-try)]
             [get-params (second fit-to-try)]
             [fit/with-params (third fit-to-try)]
             [params (get-params pts-x pts-y)]
             [fun (apply fit/with-params params)]
             [sse (sse pts-x pts-y fun)])
        (cons (Fit params sse type) fits))))
  (define fits/with-poly
    (for/fold ([fits : (Listof Fit) fits])
              ([degree : Positive-Integer (in-range 1 poly-max-degree)])
      (let* ([params (poly-fit-params pts-x pts-y degree)]
             [fun (poly-fit/with-params params)]
             [sse (sse pts-x pts-y fun)])
        (cons (Fit params sse 'poly) fits))))
  fits/with-poly)

(: fit-format-str : (-> Symbol String))
(define (fit-format-str type)
  (match type
    ('linear "a + b*x")
    ('exp "a*e^(b*x)")
    ('log "a + b*log(x)")
    ('power "(e^a)*(x^b)")
    ('poly "a + b*x + c*x^2 + ...")
    (_ "unknown")))

(: fit->string : (-> Fit String))
(define (fit->string fit)
  (format (string-append (fit-format-str (Fit-type fit))
                         " with a, b, ... = ~a (SSE: ~a)")
          (Fit-params fit)
          (Fit-sse fit)))

(: fit->fun : (-> Fit (-> Real Real)))
(define (fit->fun fit)
  (define type (Fit-type fit))
  (define fit/with-params (match type
                            ['linear linear-fit/with-params]
                            ['exp exp-fit/with-params]
                            ['log log-fit/with-params]
                            ['power power-fit/with-params]))
  (match type
    ['poly (poly-fit/with-params (Fit-params fit))]
    [_
     (match-define (list a b) (Fit-params fit))
     (fit/with-params a b)]))

(: best-fit : (-> Flonums Flonums Fit))
(define (best-fit pts-x pts-y)
  (define all-fits (try-fits pts-x pts-y))
  (for/fold ([best-so-far (first all-fits)])
            ([fit all-fits])
    (if (< (Fit-sse fit) (Fit-sse best-so-far))
        fit
        best-so-far)))
