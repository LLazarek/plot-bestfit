#lang typed/racket

(require "common.rkt")
(provide exp-fit-params
         exp-fit-params/list
         exp-fit/with-params
         exp-fit
         graph/exponential)

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
