#lang scribble/doc
@(require (except-in "base.rkt" ->)
          plot/utils plot/no-gui)
@(require (for-label
           (only-in math/flonum fl)
           (only-in typed/racket/base
                    Listof Nonnegative-Flonum Values -> Flonum Real)
           (only-in plot/utils renderer2d?)))
@title{Bestfit: Lines of Best Fit}
@defmodule[bestfit #:use-sources (bestfit)]

@author[(author+email "Spencer Florence" "spencer@florence.io")]

@(define eval
   (let ()
     (define e (make-base-eval))
     (e '(require racket))
     (e '(require bestfit))
     (e '(require plot/pict))
     (e '(require math/flonum))
     e))

@(define-syntax-rule (inter . a) (interaction #:eval eval . a))

Bestfit is a library for calculating lines of best fit using
@hyperlink["http://mathworld.wolfram.com/LeastSquaresFitting.html"]{Least Squares Fitting}.

@table-of-contents[]

@deftogether[(@defproc[(graph/linear [xs (Listof Nonnegative-Flonum)]
                                     [ys (Listof Nonnegative-Flonum)]
                                     [errors (U #f (Listof Flonum)) #f])
                       (Values renderer2d? renderer2d? renderer2d?)]
              @defproc[(graph/exponential [xs (Listof Nonnegative-Flonum)]
                                          [ys (Listof Nonnegative-Flonum)]
                                          [errors (U #f (Listof Flonum)) #f])
                       (Values renderer2d? renderer2d? renderer2d?)]
              @defproc[(graph/log [xs (Listof Nonnegative-Flonum)]
                                  [ys (Listof Nonnegative-Flonum)]
                                  [errors (U #f (Listof Flonum)) #f])
                       (Values renderer2d? renderer2d? renderer2d?)]
              @defproc[(graph/power [xs (Listof Nonnegative-Flonum)]
                                    [ys (Listof Nonnegative-Flonum)]
                                    [errors (U #f (Listof Flonum)) #f])
                       (Values renderer2d? renderer2d? renderer2d?)]
	      @defproc[(graph/poly [xs (Listof Nonnegative-Flonum)]
                                   [ys (Listof Nonnegative-Flonum)]
				   [degree Positive-Integer]
                                   [errors (U #f (Listof Flonum)) #f])
                       (Values renderer2d? renderer2d? renderer2d?)])]{

Uses @racket[linear-fit], @racket[exp-fit], @racket[log-fit], @racket[power-fit], to generate three
@racket[renderer2d?]s: A plot of the points given by @racket[xs]
and @racket[ys], a plot of the function of best fit, and error bars generated. The error bars are generated from @racket[error], which
is the percentage error on each y coordinate.

@inter[;(: 3x^2 : Nonnegative-Flonum -> Nonnegative-Flonum)
       (define (3x^2 x) (* 3.0 (expt x 2.0)))
       ;(: apply-error : Nonnegative-Flonum -> Nonnegative-Flonum)
       (define (add-error y) (+ y (* y (/ (- (random 4) 2) 10.0))))
       (define exact (function 3x^2 #:label "exact" #:color "blue"))
       (define-values (pts fit _)
         (graph/power (build-list 10 (compose fl add1))
                      (build-list 10 (compose 3x^2 fl add1))))
       (plot (list exact fit pts))
       (define-values (pts fit err)
         (graph/power (build-list 10 (compose fl add1))
                      (build-list 10 (compose add-error 3x^2 fl add1))
                      (build-list 10 (const 0.2))))
       (plot (list exact fit pts err))]

}

@deftogether[(@defproc[(linear-fit [xs (Listof Nonnegative-Flonum)]
                                   [ys (Listof Nonnegative-Flonum)])
                       (-> Nonnegative-Flonum Real)]
              @defproc[(exp-fit [xs (Listof Nonnegative-Flonum)]
                                [ys (Listof Nonnegative-Flonum)])
                      (-> Nonnegative-Flonum Real)]
              @defproc[(log-fit [xs (Listof Nonnegative-Flonum)]
                                [ys (Listof Nonnegative-Flonum)])
                      (-> Nonnegative-Flonum Real)]
              @defproc[(power-fit [xs (Listof Nonnegative-Flonum)]
                                  [ys (Listof Nonnegative-Flonum)])
                       (-> Nonnegative-Flonum Real)]
	      @defproc[(poly-fit [xs (Listof Nonnegative-Flonum)]
                                 [ys (Listof Nonnegative-Flonum)]
				 [degree Positive-Integer])
                       (-> Nonnegative-Flonum Real)])]{


Uses @hyperlink["http://mathworld.wolfram.com/LeastSquaresFitting.html"]{Least Squares Fitting} to
generate a best fit function of the given type.

@inter[(define line (linear-fit '(1.0 2.0 3.0) '(1.0 2.0 3.0)))
       (line 10.0)
       (line 12.0)]

@deftogether[(@defproc[(linear-fit-params [xs (Listof Nonnegative-Flonum)]
                                          [ys (Listof Nonnegative-Flonum)])
                       (Values Real Real)]
	      @defproc[(linear-fit-params/list [xs (Listof Nonnegative-Flonum)]
                                               [ys (Listof Nonnegative-Flonum)])
                       (List Real Real)])]{

Returns values @racket[_a] and @racket[_b] used to generate a best fit
function of the form @tt{y = a + b x}.

@inter[(linear-fit-params '(1.0 2.0 3.0) '(2.0 4.0 6.0))]

}


@deftogether[(@defproc[(exp-fit-params [xs (Listof Nonnegative-Flonum)]
                                       [ys (Listof Nonnegative-Flonum)])
                       (Values Real Real)]
	      @defproc[(exp-fit-params/list [xs (Listof Nonnegative-Flonum)]
                                            [ys (Listof Nonnegative-Flonum)])
                       (List Real Real)])]{

Returns values @racket[_A] and @racket[_B] used to generate a best fit
function of the form @tt{y = A e^(B x)}.

@inter[(exp-fit-params '(1.0 2.0 3.0) '(2.0 4.0 8.0))]

}


@deftogether[(@defproc[(log-fit-params [xs (Listof Nonnegative-Flonum)]
                                       [ys (Listof Nonnegative-Flonum)])
                       (Values Real Real)]
	      @defproc[(log-fit-params/list [xs (Listof Nonnegative-Flonum)]
                                            [ys (Listof Nonnegative-Flonum)])
                       (List Real Real)])]{

Returns values @racket[_a] and @racket[_b] used to generate a best fit
function of the form @tt{y = a + b ln(x)}.

@inter[(log-fit-params '(2.0 4.0 8.0) '(1.0 2.0 3.0))]

}


@deftogether[(@defproc[(power-fit-params [xs (Listof Nonnegative-Flonum)]
                                         [ys (Listof Nonnegative-Flonum)])
                       (Values Real Real)]
	      @defproc[(power-fit-params/list [xs (Listof Nonnegative-Flonum)]
                                              [ys (Listof Nonnegative-Flonum)])
                       (List Real Real)])]{

Returns values @racket[_A] and @racket[_B] used to generate a best fit
function of the form @tt{y = A x^B}.

@inter[(power-fit-params '(1.0 2.0 3.0) '(2.0 8.0 18.0))]

}


@defproc[(poly-fit-params [xs (Listof Nonnegative-Flonum)]
                          [ys (Listof Nonnegative-Flonum)]
			  [degree Positive-Integer])
         (Listof Real)]{

Returns values @racket[_A], @racket[_B], @racket[_C], ... used to
generate a best fit function of the form @tt{y = A + Bx + Cx^2 + ...}.

@inter[(poly-fit-params '(1.0 2.0 3.0) '(1.0 8.0 27.0) 3)]

}


}



@section{Fit Exploration}
You may not know which of the above fits applies best to your data;
these functions allow you to automatically find the best fit function
for a given data set.

@defproc[(try-fits [xs (Listof Nonnegative-Flonum)]
                   [ys (Listof Nonnegative-Flonum)]
		   [#:poly-max-degree poly-max-degree Natural 0])
	 (Listof Fit)]{

Try all of the above described fits on the given data, returning a
@tt{Fit} object for each attempted fit.

If @racket[poly-max-degree] is 0, do not try any polynomial fits. If
it is some other number @tt{N}, try all polynomials from degree 2
through @tt{N}.


@inter[(try-fits '(1.0 2.0 3.0) '(2.0 8.0 18.0))]
@inter[(try-fits '(1.0 2.0 3.0) '(1.0 8.0 27.0) 3)]

}


@defproc[(best-fit [xs (Listof Nonnegative-Flonum)]
                   [ys (Listof Nonnegative-Flonum)]
		   [#:poly-max-degree poly-max-degree Natural 0])
	 Fit]{

Try all of the above described fits on the given data, returning the
fit with the lowest error between the fit function and the data.

If @racket[poly-max-degree] is 0, do not try any polynomial fits. If
it is some other number @tt{N}, try all polynomials from degree 2
through @tt{N}.


@inter[(best-fit '(1.0 2.0 3.0) '(2.0 8.0 18.0))]
@inter[(best-fit '(1.0 2.0 3.0) '(1.0 8.0 27.0) #:poly-max-degree 3)]

}

@defproc[(fit->string [fit Fit])
         String]{

Returns a human-readable string representation of the given fit, along
with its parameter values and error.


@inter[(fit->string (best-fit '(1.0 2.0 3.0) '(2.0 8.0 18.0)))]

}

@defproc[(fit->fun [fit Fit])
         String]{

Return a procedure implimenting the given fit function.


@inter[(fit->string (best-fit '(1.0 2.0 3.0) '(2.0 8.0 18.0)))]

}

@deftogether[(
@defproc[(Fit-params [fit Fit]) (Listof Real)]
@defproc[(Fit-sse [fit Fit]) Real]
@defproc[(Fit-type [fit Fit]) Symbol]
)]{

Access relevant information (parameter values, total fit error, and
type, respectively) about a @tt{Fit} as returned by @racket[try-fits]
and @racket[best-fit].

}

