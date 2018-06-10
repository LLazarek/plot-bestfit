#lang typed/racket

(require math/matrix
         math/flonum)

(define-type Flonums (Listof Nonnegative-Flonum))

;; See http://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
(define (polyfit [xs : Flonums]
                 [ys : Flonums]
                 [degree : Positive-Integer])
  (define n (length xs))
  (define X (vandermonde-matrix xs (add1 degree)))
  (define Y (list->matrix n 1 ys))
  (define X^T (matrix-transpose X))
  (define X-bar (matrix* X^T X))
  (define Y-bar (matrix* X^T Y))
  ;; Solve X-bar * alpha = Y-bar to find coefficients

  (define X-bar^T (matrix-transpose X-bar))
  (define X-bar^T*X-bar (matrix* X-bar^T X-bar))
  (define X-bar^T*X-bar^-1 (matrix-inverse X-bar^T*X-bar))
  (matrix->list (matrix* X-bar^T*X-bar^-1
                         X-bar^T
                         Y-bar)))

;; todo: try
;; Using LU decomp: http://math.oit.edu/~watermang/math_341/341_ch7/F13_341_book_sec_7-2.pdf
;; Using Cramer's rule: https://en.wikipedia.org/wiki/Cramer%27s_rule

;; I could also consider making my gradient descent program use
;; matrices instead to be fast, since it seems to work pretty ok but
;; is way too slow

(module+ test
  (require typed/rackunit)
  (define-syntax-rule (check-~ x y threshold)
    (check-true (< (abs (- x y)) threshold)))
  (define ~-threshold 0.1)


  ;; (define xs : Flonums (drop (build-list 1000 fl) 1))
  (define xs : Flonums (drop (build-list 11 fl) 1))

  (define xs^2 : Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x x))
         xs))
  (match-define (list sq/zero-mult sq/one-mult sq/two-mult)
    (polyfit xs xs^2 2))
  (check-~ sq/zero-mult 0 ~-threshold)
  (check-~ sq/one-mult 0 ~-threshold)
  (check-~ sq/two-mult 1 ~-threshold)


  (define xs*2 : Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x 2))
         xs))
  (match-define (list double/zero-mult
                      double/one-mult
                      double/two-mult)
    (polyfit xs xs*2 2))
  (check-~ double/zero-mult 0 ~-threshold)
  (check-~ double/one-mult 2 ~-threshold)
  (check-~ double/two-mult 0 ~-threshold)

  (match-define (list double^1/zero-mult
                      double^1/one-mult)
    (polyfit xs xs*2 1))
  (check-~ double^1/zero-mult 0 ~-threshold)
  (check-~ double^1/one-mult 2 ~-threshold)


  (define xs^2*2 : Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x 2))
         xs^2))
  (match-define (list double-sq/zero-mult
                        double-sq/one-mult
                        double-sq/two-mult)
    (polyfit xs xs^2*2 2))
  (check-~ double-sq/zero-mult 0 ~-threshold)
  (check-~ double-sq/one-mult 0 ~-threshold)
  (check-~ double-sq/two-mult 2 ~-threshold)


  (define ~-threshold/higher 0.2)
  (define xs^3 :  Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x x x))
         xs))
  (match-define (list cube/zero-mult
                      cube/one-mult
                      cube/two-mult
                      cube/three-mult)
    (polyfit xs xs^3 3))
  (check-~ cube/zero-mult 0 ~-threshold/higher)
  (check-~ cube/one-mult 0 ~-threshold/higher)
  (check-~ cube/two-mult 0 ~-threshold/higher)
  (check-~ cube/three-mult 1 ~-threshold/higher)

  (define xs^3*17 :  Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* 17 x))
         xs^3))
  (match-define (list cube17/zero-mult
                      cube17/one-mult
                      cube17/two-mult
                      cube17/three-mult)
    (polyfit xs xs^3*17 3))
  (check-~ cube17/zero-mult 0 ~-threshold/higher)
  (check-~ cube17/one-mult 0 ~-threshold/higher)
  (check-~ cube17/two-mult 0 ~-threshold/higher)
  (check-~ cube17/three-mult 17 ~-threshold/higher)

  (define xs^4 :  Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x x x x))
         xs))
  (match-define (list quad/zero-mult
                      quad/one-mult
                      quad/two-mult
                      quad/three-mult
                      quad/four-mult)
    (polyfit xs xs^4 4))
  (check-~ quad/zero-mult 0 ~-threshold/higher)
  (check-~ quad/one-mult 0 ~-threshold/higher)
  (check-~ quad/two-mult 0 ~-threshold/higher)
  (check-~ quad/three-mult 0 ~-threshold/higher)
  (check-~ quad/four-mult 1 ~-threshold/higher))
