#lang typed/racket

(require "common.rkt")
(provide power-fit-params
         power-fit-params/list
         power-fit/with-params
         power-fit
         graph/power)

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


;; ---------- Tests ----------
(module+ test
  (require "test-common.rkt")
  (define xs/power : Flonums (map fl (range 1 500)))
  (define ys/power : Flonums (cast (map (λ ([x : Flonum])
                                          (* (exp 0.5) (expt x 0.2)))
                                        xs/power)
                                    Flonums))
  (define-values (a b) (power-fit-params xs/power ys/power))
  (check-≈ a 0.5 0.0001)
  (check-≈ b 0.2 0.0001))
