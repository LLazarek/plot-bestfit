#lang typed/racket

(require "common.rkt")
(provide linear-fit-params
         linear-fit-params/list
         linear-fit/with-params
         linear-fit
         graph/linear)

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
