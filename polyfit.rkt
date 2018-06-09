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

(define (gradient-descent/degree-2 [xs : Flonums]
                                   [ys : Flonums]
                                   [maxiters : Positive-Integer 1000000]
                                   #:fix-threshold
                                   [fix-threshold : Positive-Real 0.0000001]
                                   #:alpha [alpha : Positive-Real 0.0005])
  : (values Flonum Flonum Flonum)
  (define m (length xs))

  (let loop ([iter 0]
             [last-err +inf.0]
             [theta-0 0.5]
             [theta-1 0.5]
             [theta-2 0.5])
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
    (if (or (< (- last-err err) fix-threshold) (>= iter maxiters))
        (values theta-0 theta-1 theta-2)
        (loop (+ iter 1) err theta-0/new theta-1/new theta-2/new))))

(define (random-list [len : Positive-Integer]) : (Listof Flonum)
  (define (random/-1->1 [_ : Any]) : Flonum
    (fl (/ (- (random 200) 100) 100.)))
  (build-list len random/-1->1))

;; If no theta-guess provided then a set initial values is randomly chosen
(define (gradient-descent [xs : Flonums]
                          [ys : Flonums]
                          [degree : Nonnegative-Integer]
                          [maxiters : Positive-Integer 300000]
                          #:fix-threshold
                          [fix-threshold : Positive-Real 0.0000001]
                          #:alpha [alpha : Positive-Real 0.0005]
                          #:regularize [reg-λ : Nonnegative-Real 0.0]
                          #:guess [theta-guess : (Listof Flonum) empty])
  : (Listof Flonum)
  (define m (length xs))
  (define theta-init (if (empty? theta-guess)
                         (random-list (add1 degree))
                         theta-guess))
  (unless (= (length theta-init) (add1 degree))
    (error "Coefficient guess must provide a value for each coefficient,
i.e. one more than the degree to fit."))

  (let loop ([iter 1]
             [last-err +inf.0]
             [theta : (Listof Flonum) theta-init])
    (define (h-theta [x : Flonum]) : Flonum
      (flsum
       (for/list ([theta-i (in-list theta)]
                  [power (in-naturals)])
         (* theta-i (flexpt x (fl power))))))

    (define theta/new : (Listof Flonum)
      (for/list ([theta-i (in-list theta)]
                 [power (in-naturals)])
        (fl (- (* theta-i
                  (- 1 (/ (* alpha reg-λ) m)))
               (* (/ alpha m)
                  (flsum (for/list ([x (in-list xs)]
                                    [y (in-list ys)])
                           (* (- (h-theta x) y)
                              (flexpt x (fl power))))))))))

    (match-define (list zipped _ _) (zip (map h-theta xs) ys))
    (define err (apply sse zipped))
    (if (or (< (abs (- last-err err)) fix-threshold) (>= iter maxiters))
        (begin (printf "Last iter: ~v; err = ~v (last = ~v)\n"
                       iter err last-err)
               theta/new)
        (loop (+ iter 1) err theta/new))))


(module+ test
  (require typed/rackunit)
  (define-syntax-rule (check-~ x y threshold)
    (check-true (< (abs (- x y)) threshold)))
  (define ~-threshold 0.01)
  (define xs (map fl '(1 2 3 4 5 6 7 8 9 10)))
  (define ys (map fl '(1 4 9 16 25 36 49 64 81 100)))
  (match-define (list sq/zero-mult sq/one-mult sq/two-mult)
                       (gradient-descent xs ys 2))
  (check-~ sq/zero-mult 0 ~-threshold)
  (check-~ sq/one-mult 0 ~-threshold)
  (check-~ sq/two-mult 1 ~-threshold)

  (define xs*2 : Flonums (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
                                 (* x 2)) xs))
  (match-define (list double/zero-mult double/one-mult double/two-mult)
                       (gradient-descent xs xs*2 2))
  (check-~ double/zero-mult 0 ~-threshold)
  (check-~ double/one-mult 2 ~-threshold)
  (check-~ double/two-mult 0 ~-threshold)

  (define ys*2 : Flonums (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
                                 (* x 2)) ys))
  (match-define (list double-sq/zero-mult
                        double-sq/one-mult
                        double-sq/two-mult)
                       (gradient-descent xs ys*2 2))
  (check-~ double-sq/zero-mult 0 ~-threshold)
  (check-~ double-sq/one-mult 0 ~-threshold)
  (check-~ double-sq/two-mult 2 ~-threshold))
