#lang typed/racket

(require plot/no-gui math/flonum)

(provide (all-from-out plot/no-gui)
         (all-from-out math/flonum)
         Flonums
         Grapher
         Grapher/poly
         Fitter
         Fitter/poly
         Parameterizer
         Parameterizer*
         graph/gen
         Σ)

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

