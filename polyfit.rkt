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

(define (gradient-descent/auto-tune [xs : Flonums]
                                    [ys : Flonums]
                                    [degree : Nonnegative-Integer]
                                    #:fix-threshold
                                    [fix-threshold : Positive-Real 0.0000001]
                                    #:guess
                                    [theta-guess : (Listof Flonum) empty])
  : (Listof Flonum)
  ;; Hack: Higher-degree polynomials seem to do better with smaller alphas
  (gradient-descent xs ys degree
                    #:alpha (* 5 (expt 100 (- degree)))
                    #:fix-threshold fix-threshold
                    #:guess theta-guess))

(provide gradient-descent
         gradient-descent/auto-tune)

(module+ test
  (require typed/rackunit)
  (define-syntax-rule (check-~ x y threshold)
    (check-true (< (abs (- x y)) threshold)))
  (define ~-threshold 0.01)


  (define xs : Flonums
    (drop (build-list 11 fl) 1))

  (define xs^2 : Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x x))
         xs))
  (match-define (list sq/zero-mult sq/one-mult sq/two-mult)
    (gradient-descent xs xs^2 2 #:guess '(0.2 0.2 0.7)))
  (check-~ sq/zero-mult 0 ~-threshold)
  (check-~ sq/one-mult 0 ~-threshold)
  (check-~ sq/two-mult 1 ~-threshold)


  (define xs*2 : Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x 2))
         xs))
  (match-define (list double/zero-mult double/one-mult double/two-mult)
    (gradient-descent xs xs*2 2 #:guess '(0.2 1.5 0.2)))
  (check-~ double/zero-mult 0 ~-threshold)
  (check-~ double/one-mult 2 ~-threshold)
  (check-~ double/two-mult 0 ~-threshold)


  (define xs^2*2 : Flonums
    (map (λ ([x : Nonnegative-Flonum]) : Nonnegative-Flonum
            (* x 2))
         xs^2))
  (match-define (list double-sq/zero-mult
                        double-sq/one-mult
                        double-sq/two-mult)
    (gradient-descent xs xs^2*2 2 #:guess '(0.2 0.2 1.5)))
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
    (gradient-descent/auto-tune xs xs^3 3
                                #:guess '(0.2 0.2 0.3 0.7)))
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
    (gradient-descent/auto-tune xs xs^3*17 3
                                #:guess '(0.2 0.2 0.3 5.0)))
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
    (gradient-descent/auto-tune xs xs^4 4
                                #:guess '(0.2 0.2 0.1 -0.2 0.7)))
  (check-~ quad/zero-mult 0 ~-threshold/higher)
  (check-~ quad/one-mult 0 ~-threshold/higher)
  (check-~ quad/two-mult 0 ~-threshold/higher)
  (check-~ quad/three-mult 0 ~-threshold/higher)
  (check-~ quad/four-mult 1 ~-threshold/higher))
