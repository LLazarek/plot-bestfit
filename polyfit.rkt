#lang typed/racket

(require math/flonum
         racket/match)

(define-type Flonums (Listof Nonnegative-Flonum))


(define #:forall (A B)
  (zip [list1 : (Listof A)]
       [list2 : (Listof B)]) : (List (Listof (List A B))
                                     (Listof A)
                                     (Listof B))
  (let zip-rest ([l1 list1]
                 [l2 list2]
                 [result : (Listof (List A B)) empty])
    (match (list l1 l2)
      [(list '() remainder)
       (list (reverse result) empty remainder)]
      [(list remainder '())
       (list (reverse result) remainder empty)]

      [(list (cons head1 tail1) (cons head2 tail2))
       (zip-rest tail1 tail2
                 (cons (list head1 head2) result))])))

(define (sse . [pairs : (List Flonum Flonum) *]) : Flonum
  (define (sse/pair [pair : (List Flonum Flonum)])
    (expt (- (first pair) (second pair)) 2))
  (define sses (map sse/pair pairs))
  (define sses/summed (flsum sses))
  (flsqrt sses/summed))

(define (gradient-descent [xs : Flonums]
                          [ys : Flonums]) : (values Flonum Flonum Flonum)
  (define alpha 0.0005)
  (define m (length xs))

  (let loop ([theta-0 0.5]
             [theta-1 0.5]
             [theta-2 0.5]
             [last-err +inf.0])
    (define (h-theta [x : Flonum]) : Flonum
      (+ (* theta-0 (expt x 0))
         (* theta-1 (expt x 1))
         (* theta-2 (expt x 2))))

    (define theta-0/new
      (- theta-0
         (* (/ alpha m)
            (flsum (for/list ([x (in-list xs)]
                              [y (in-list ys)])
                     (* (- (h-theta x) y)
                        (flexpt x (fl 0))))))))
    (define theta-1/new
      (- theta-1
         (* (/ alpha m)
            (flsum (for/list ([x (in-list xs)]
                              [y (in-list ys)])
                     (* (- (h-theta x) y)
                        (flexpt x (fl 1))))))))
    (define theta-2/new
      (- theta-2
         (* (/ alpha m)
            (flsum (for/list ([x (in-list xs)]
                              [y (in-list ys)])
                     (* (- (h-theta x) y)
                        (flexpt x (fl 2))))))))


    (match-define (list zipped _ _) (zip (map h-theta xs) ys))
    (define err (apply sse zipped))
    (assert (< err last-err))
    (if (< (- last-err err) 0.0000001)
        (values theta-0 theta-1 theta-2)
        (loop theta-0/new theta-1/new theta-2/new err))))

(module+ test
  (require typed/rackunit)
  (define-syntax-rule (check-~ x y threshold)
    (check-true (< (abs (- x y)) threshold)))
  (define ~-threshold 0.01)
  (define xs (map fl '(1 2 3 4 5 6 7 8 9 10)))
  (define ys (map fl '(1 4 9 16 25 36 49 64 81 100)))
  (match-define-values (sq/zero-mult sq/one-mult sq/two-mult)
                       (gradient-descent xs ys))
  (check-~ sq/zero-mult 0 ~-threshold)
  (check-~ sq/one-mult 0 ~-threshold)
  (check-~ sq/two-mult 1 ~-threshold)

  (define xs*2 : Flonums (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
                                 (* x 2)) xs))
  (match-define-values (double/zero-mult double/one-mult double/two-mult)
                       (gradient-descent xs xs*2))
  (check-~ double/zero-mult 0 ~-threshold)
  (check-~ double/one-mult 2 ~-threshold)
  (check-~ double/two-mult 0 ~-threshold)

  (define ys*2 : Flonums (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
                                 (* x 2)) ys))
  (match-define-values (double-sq/zero-mult
                        double-sq/one-mult
                        double-sq/two-mult)
                       (gradient-descent xs ys*2))
  (check-~ double-sq/zero-mult 0 ~-threshold)
  (check-~ double-sq/one-mult 0 ~-threshold)
  (check-~ double-sq/two-mult 2 ~-threshold)
  )
