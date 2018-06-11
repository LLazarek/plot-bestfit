#lang typed/racket

(require typed/rackunit)
(define-syntax-rule (check-≈ x y threshold)
  (check-true (< (abs (- x y)) threshold)))

(provide (all-from-out typed/rackunit)
         check-≈)
