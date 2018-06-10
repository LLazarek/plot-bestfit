#lang typed/racket

(require "common.rkt")
(provide log-fit-params
         log-fit-params/list
         log-fit/with-params
         log-fit
         graph/log)

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


;; ---------- Tests ----------
(module+ test
  (define xs/log : Flonums (map fl (range 1 500)))
  (define ys/log : Flonums (cast (map (λ ([x : Flonum])
                                        (+ 5 (* 2 (fllog x)))) xs/log)
                                 Flonums))
 )
