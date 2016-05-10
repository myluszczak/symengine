/**
 *  \file series_generic.h
 *  Class for univariate series.
 *
 **/
#ifndef SYMENGINE_SERIES_GENERIC_H
#define SYMENGINE_SERIES_GENERIC_H

#include <list>
#include <symengine/polynomial.h>
#include <symengine/rational.h>
#include <symengine/series.h>

namespace SymEngine
{
//! UnivariateSeries Class
class UnivariateSeries
    : public SeriesBase<UnivariateExprPolynomial, Expression, UnivariateSeries>
{
    // UnivariateSeries 1 + 2*x + x**2 + O(x**5) has dict_ = {{0, 1}, {1, 2},
    // {2, 1}} with var_ = "x" and prec_ = 5
public:
    IMPLEMENT_TYPEID(UNIVARIATESERIES)
    UnivariateSeries(const UnivariateExprPolynomial &sp,
                     const std::string varname, const unsigned degree)
        : SeriesBase(std::move(sp), varname, degree)
    {
    }

    static RCP<const UnivariateSeries> create(const RCP<const Symbol> &var,
                                              const unsigned int &prec,
                                              const UnivariateExprPolynomial &s)
    {
        return make_rcp<const UnivariateSeries>(std::move(s), var->get_name(),
                                                prec);
    }

    static RCP<const UnivariateSeries>
    series(const RCP<const Basic> &t, const std::string &x, unsigned int prec);
    virtual std::size_t __hash__() const;
    virtual int compare(const Basic &o) const;
    bool operator==(const UnivariateSeries &u) const;
    virtual RCP<const Basic> as_basic() const;
    virtual umap_int_basic as_dict() const;
    virtual RCP<const Basic> get_coeff(int) const;
    static UnivariateExprPolynomial var(const std::string &s);

    static Expression convert(const Basic &x);

    static int ldegree(const UnivariateExprPolynomial &s);
    static UnivariateExprPolynomial mul(const UnivariateExprPolynomial &s,
                                        const UnivariateExprPolynomial &r,
                                        unsigned prec);
    static UnivariateExprPolynomial pow(const UnivariateExprPolynomial &s,
                                        int n, unsigned prec);
    static Expression find_cf(const UnivariateExprPolynomial &s,
                              const UnivariateExprPolynomial &var, int deg);
    static Expression root(Expression &c, unsigned n);
    static UnivariateExprPolynomial diff(const UnivariateExprPolynomial &s,
                                         const UnivariateExprPolynomial &var);
    static UnivariateExprPolynomial
    integrate(const UnivariateExprPolynomial &s,
              const UnivariateExprPolynomial &var);
    static UnivariateExprPolynomial subs(const UnivariateExprPolynomial &s,
                                         const UnivariateExprPolynomial &var,
                                         const UnivariateExprPolynomial &r,
                                         unsigned prec);

    static Expression sin(const Expression &c);
    static Expression cos(const Expression &c);
    static Expression tan(const Expression &c);
    static Expression asin(const Expression &c);
    static Expression acos(const Expression &c);
    static Expression atan(const Expression &c);
    static Expression sinh(const Expression &c);
    static Expression cosh(const Expression &c);
    static Expression tanh(const Expression &c);
    static Expression asinh(const Expression &c);
    static Expression acosh(const Expression &c);
    static Expression atanh(const Expression &c);
    static Expression exp(const Expression &c);
    static Expression log(const Expression &c);
};

inline RCP<const UnivariateSeries>
univariate_series(RCP<const Symbol> i, unsigned int prec,
                  const UnivariateExprPolynomial &s)
{
    return make_rcp<const UnivariateSeries>(std::move(s), i->get_name(), prec);
}

//! MultivariateSeries Class
class MultivariateSeries : public SeriesBase<MultivariateExprPolynomial,
                                             Expression, MultivariateSeries>
{
public:
    map_sym_uint precs_;

public:
    IMPLEMENT_TYPEID(MULTIVARIATESERIES)
    MultivariateSeries(const MultivariateExprPolynomial &sp,
                       const std::string varname, const unsigned degree);

    MultivariateSeries(const MultivariateExprPolynomial &sp,
                       const std::string varname, const unsigned degree,
                       const map_sym_uint &precs);

    bool is_canonical(const MultivariateExprPolynomial p, const std::string var,
                      const long degree, const map_sym_uint precs_);

    static RCP<const MultivariateSeries>
    create(const RCP<const Symbol> &var, const unsigned int &prec,
           const MultivariateExprPolynomial &s)
    {
        return make_rcp<const MultivariateSeries>(std::move(s), var->get_name(),
                                                  prec);
    }

    static RCP<const MultivariateSeries>
    series(const RCP<const Basic> &t, const std::string &x, unsigned int prec);
    virtual std::size_t __hash__() const;
    virtual int compare(const Basic &o) const;
    virtual RCP<const Basic> as_basic() const;
    umap_int_basic as_dict() const;

    bool __eq__(const Basic &o) const;
    RCP<const Number> add(const Number &other) const;
    RCP<const Number> mul(const Number &other) const;
    // RCP<const Number> pow(const Number &other) const;

    RCP<const Basic> get_coeff(int) const;

    static MultivariateExprPolynomial add(const MultivariateExprPolynomial &s,
                                          const MultivariateExprPolynomial &r,
                                          map_sym_uint prec);
    static MultivariateExprPolynomial mul(const MultivariateExprPolynomial &s,
                                          const MultivariateExprPolynomial &r,
                                          map_sym_uint prec);
    // static MultivariateExprPolynomial pow (const MultivariateExprPolynomial
    // &s, int n, map_sym_uint prec);

    static MultivariateExprPolynomial var(const std::string &s);

    static Expression convert(const Basic &x);

    static int ldegree(const MultivariateExprPolynomial &s);
    static MultivariateExprPolynomial mul(const MultivariateExprPolynomial &s,
                                          const MultivariateExprPolynomial &r,
                                          unsigned prec);
    static MultivariateExprPolynomial pow(const MultivariateExprPolynomial &s,
                                          int n, unsigned prec);
    static Expression find_cf(const MultivariateExprPolynomial &s,
                              const MultivariateExprPolynomial &var, int deg);
    static Expression root(Expression &c, unsigned n);
    static MultivariateExprPolynomial
    diff(const MultivariateExprPolynomial &s,
         const MultivariateExprPolynomial &var);
    static MultivariateExprPolynomial
    integrate(const MultivariateExprPolynomial &s,
              const MultivariateExprPolynomial &var);
    static MultivariateExprPolynomial
    subs(const MultivariateExprPolynomial &s,
         const MultivariateExprPolynomial &var,
         const MultivariateExprPolynomial &r, unsigned prec);

    static Expression sin(const Expression &c);
    static Expression cos(const Expression &c);
    static Expression tan(const Expression &c);
    static Expression asin(const Expression &c);
    static Expression acos(const Expression &c);
    static Expression atan(const Expression &c);
    static Expression sinh(const Expression &c);
    static Expression cosh(const Expression &c);
    static Expression tanh(const Expression &c);
    static Expression asinh(const Expression &c);
    static Expression acosh(const Expression &c);
    static Expression atanh(const Expression &c);
    static Expression exp(const Expression &c);
    static Expression log(const Expression &c);
};

inline RCP<const MultivariateSeries>
multivariate_series(RCP<const Symbol> i, unsigned int prec,
                    const MultivariateExprPolynomial &s)
{
    return make_rcp<const MultivariateSeries>(std::move(s), i->get_name(),
                                              prec);
}

} // SymEngine
#endif
