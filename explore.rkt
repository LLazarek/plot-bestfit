#lang typed/racket

(require "common.rkt"
         "linear-fit.rkt"
         "exp-fit.rkt"
         "log-fit.rkt"
         "power-fit.rkt"
         "poly-fit.rkt")
(provide Fit try-fits best-fit
         fit->string fit->fun)

;; --------------------
;; | Exploratory fit
;; Try every fit function, returning a summary or just the best

;; For now, just always use SSE for error
;; ;; f(x) = offset + slope*x
;; (define-struct LinearFit ([offset : Real] [slope : Real] [sse : Real]))
;; ;; f(x) = coeff*e^(exp-coeff*x)
;; (define-struct ExpFit ([coeff : Real] [exp-coeff : Real] [sse : Real]))
;; ;; f(x) = offset + coeff*log(x)
;; (define-struct LogFit ([offset : Real] [coeff : Real] [sse : Real]))
;; ;; f(x) = (e^coeff)*(x^exp-coeff)
;; (define-struct PowerFit ([coeff : Real] [exp-coeff : Real] [sse : Real]))

(: point-error : (-> Real Real Real))
(define (point-error actual expected)
  (expt (- expected actual) 2))

(: sse : (-> Flonums Flonums (-> Real Real) Real))
(define (sse pts-x pts-y f)
  (foldl + 0
         (map point-error
              (map f pts-x)
              pts-y)))

(struct Fit ([params : (Listof Real)] [sse : Real] [type : Symbol]))

(: try-fits : (->* (Flonums Flonums)
                   (#:poly-max-degree Natural)
                   (Listof Fit)))
(define (try-fits pts-x pts-y
                  #:poly-max-degree [poly-max-degree 0])
  (define fits-to-try `((linear ,linear-fit-params/list ,linear-fit/with-params)
                        (exp ,exp-fit-params/list ,exp-fit/with-params)
                        (log ,log-fit-params/list ,log-fit/with-params)
                        (power ,power-fit-params/list ,power-fit/with-params)))
  (define fits
    (for/fold ([fits : (Listof Fit) empty])
              ([fit-to-try fits-to-try])
      (let* ([type (first fit-to-try)]
             [get-params (second fit-to-try)]
             [fit/with-params (third fit-to-try)]
             [params (get-params pts-x pts-y)]
             [fun (apply fit/with-params params)]
             [sse (sse pts-x pts-y fun)])
        (cons (Fit params sse type) fits))))
  (define fits/with-poly
    (for/fold ([fits : (Listof Fit) fits])
              ([degree : Positive-Integer (in-range 1 poly-max-degree)])
      (let* ([params (poly-fit-params pts-x pts-y degree)]
             [fun (poly-fit/with-params params)]
             [sse (sse pts-x pts-y fun)])
        (cons (Fit params sse 'poly) fits))))
  fits/with-poly)

(: fit-format-str : (-> Symbol String))
(define (fit-format-str type)
  (match type
    ('linear "a + b*x")
    ('exp "a*e^(b*x)")
    ('log "a + b*log(x)")
    ('power "(e^a)*(x^b)")
    ('poly "a + b*x + c*x^2 + ...")
    (_ "unknown")))

(: fit->string : (-> Fit String))
(define (fit->string fit)
  (format (string-append (fit-format-str (Fit-type fit))
                         " with a, b, ... = ~a (SSE: ~a)")
          (Fit-params fit)
          (Fit-sse fit)))

(: fit->fun : (-> Fit (-> Real Real)))
(define (fit->fun fit)
  (define type (Fit-type fit))
  (define fit/with-params (match type
                            ['linear linear-fit/with-params]
                            ['exp exp-fit/with-params]
                            ['log log-fit/with-params]
                            ['power power-fit/with-params]))
  (match type
    ['poly (poly-fit/with-params (Fit-params fit))]
    [_
     (match-define (list a b) (Fit-params fit))
     (fit/with-params a b)]))

(: best-fit : (->* (Flonums Flonums)
                   (#:poly-max-degree Natural)
                   Fit))
(define (best-fit pts-x pts-y #:poly-max-degree [poly-max-degree 0])
  (define all-fits (try-fits pts-x pts-y #:poly-max-degree poly-max-degree))
  (define (valid-fit [f : Fit])
    (define sse (Fit-sse f))
    (not (or (nan? sse) (infinite? sse))))
  (define valid-fits (filter valid-fit all-fits))
  (for/fold ([best-so-far (first valid-fits)])
            ([fit (rest valid-fits)])
    (if (< (Fit-sse fit) (Fit-sse best-so-far))
        (begin
          (printf "Decided that ~v\nis better than ~v\n\n"
                  (fit->string fit) (fit->string best-so-far))
          fit)
        (begin
          (printf "Decided that ~v\nis better than ~v\n\n"
                  (fit->string best-so-far) (fit->string fit))
          best-so-far))))
